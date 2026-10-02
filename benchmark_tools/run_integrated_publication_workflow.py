"""Install, infer, read back and score with separate frozen runtime environments.

Run this stdlib-only controller by absolute path with isolated Python. Assets
and acquired data must be supplied locally; it does not download or retry.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys

LOCKS = {
    "inference": "813751262c9405155cc2f1974a2c32187b352f88307e2332e0cec9df35708ae3",
    "reader": "53df1cabd91c2fb18179951ebc95166ae49ba843f4015e0f417468d2a8265063",
}
HARNESS = {
    "run_publication_pipeline.py": "acab3b802b90ea2c4bc4ebc679f5f64431b1672ad85ca38f7ced212965edec87",
    "candidate_hit_order_policy.py": "cc70b884e40fe73a3c25ef9ae60a2133508127323a11382a8c610e0d192ec811",
    "probe_ob_canonical_candidates.py": "c6f55261ce34e17681f793afce5bfeebcbc615c6df585e4bce4f62628e48060e",
}


def record(path):
    path = Path(path).absolute()
    data = path.read_bytes()
    return dict(path=str(path), bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def check(item):
    if record(item["path"]) != item:
        raise ValueError("Changed pinned file: " + item["path"])


def save(path, value):
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")


def validate_assets(assets, readers, assembly_manifest_sha256=None):
    if assembly_manifest_sha256 is not None:
        if not __package__:
            sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
        from benchmark_tools.assemble_publication_runtime_assets import validate_for_executor
        return validate_for_executor(assets, readers, assembly_manifest_sha256)
    path = assets / "copied_assets.json"
    if record(path)["sha256"] != "8f8be7f1609d549da79a4ec4231e937a08833a7b8acc415d03e66b16274d5067":
        raise ValueError("Changed validated native-asset manifest")
    rows = json.loads(path.read_text())
    original_root = Path(rows[0]["copied"]["path"]).parents[1]
    for row in rows:
        item = row["copied"]
        path = Path(item["path"])
        if path.is_relative_to(original_root):
            relative = path.relative_to(original_root)
        elif path.name in {"mafft-distance", "mafft-profile"} and path.parent.parts[-2:] == ("libexec", "mafft"):
            # The frozen manifest resolved two historical absolute convenience links.
            relative = Path("mafft/bin") / path.name
        else:
            raise ValueError("Unrecognized historical asset path")
        actual_path = assets / relative
        if not actual_path.resolve().is_relative_to(assets.resolve()):
            raise ValueError("Native asset link escapes portable tree")
        actual = record(actual_path)
        if any(actual[k] != item[k] for k in ("bytes", "sha256")):
            raise ValueError("Changed validated native asset")
    for name, digest in HARNESS.items():
        if record(assets / "benchmark_tools" / name)["sha256"] != digest:
            raise ValueError("Changed frozen inference harness")
    path = readers / "manifest.json"
    if record(path)["sha256"] != "c2a77c7a2635b9193f4321f96db65028a4d5158facd8e52bb4285fe68804f334":
        raise ValueError("Changed reader export manifest")
    expected = json.loads(path.read_text())["files"]
    actual = {str(p.relative_to(readers)) for p in readers.rglob("*") if p.is_file()}
    if actual != set(expected) | {"manifest.json"}:
        raise ValueError("Changed reader export inventory")
    for name, item in expected.items():
        actual = record(readers / name)
        if any(actual[k] != item[k] for k in ("bytes", "sha256")):
            raise ValueError("Changed exported reader")


def relocate_assets(source, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    if record(source / "copied_assets.json")["sha256"] != "8f8be7f1609d549da79a4ec4231e937a08833a7b8acc415d03e66b16274d5067":
        raise ValueError("Wrong historical assets")
    output.mkdir(parents=True)
    for name in ("wheels", "mafft", "benchmark_tools"):
        shutil.copytree(source / name, output / name, symlinks=True)
    for name in ("FastTree", "copied_assets.json", "manifest.json"):
        shutil.copy2(source / name, output / name)
    changes = []
    for name in ("mafft-distance", "mafft-profile"):
        path = output / "mafft/bin" / name
        if not path.is_symlink():
            raise ValueError("Expected historical convenience link")
        previous = os.readlink(path)
        target = "../libexec/mafft/" + name
        path.unlink()
        path.symlink_to(target)
        changes.append(dict(relative_path=str(path.relative_to(output)), previous=previous, target=target,
                            payload=record(path)))
    for path in output.rglob("*"):
        if path.is_symlink() and not path.resolve().is_relative_to(output.resolve()):
            raise ValueError("Remaining escaping asset link")
    return changes


def validate_data(data):
    dimensions = {"orthobench": (251378, 12, 70, 11), "installation_fixture": (16, 4, 3, 0)}
    observed = (data["genes"], len(data["fasta"]), len(data["references"]), len(data["uncertain"]))
    if dimensions.get(data["dataset"]) != observed:
        raise ValueError("Wrong dataset scope or dimensions")
    for role in ("fasta", "references", "uncertain"):
        names = [Path(r["path"]).name for r in data[role]]
        if len(names) != len(set(names)):
            raise ValueError("Duplicate input basename")
        for row in data[role]:
            check(row)


def stage(directory, name, command, environment, timeout):
    save(directory / (name + "_started.json"), dict(command=command, environment=environment))
    with (directory / (name + ".log")).open("x") as log:
        process = subprocess.Popen(command, cwd=directory, env=environment, stdout=log,
                                   stderr=subprocess.STDOUT, start_new_session=True)
        try:
            returncode = process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait()
            save(directory / (name + "_failed.json"), dict(timeout=True, retry=False))
            raise
    result = dict(command=command, returncode=returncode, log=record(directory / (name + ".log")))
    save(directory / (name + "_finished.json"), result)
    if returncode:
        raise RuntimeError("Stage failed without retry: " + name)
    return result


def install_commands(base, installer, wheels, lock, destination):
    python = str(destination / "bin/python")
    return [
        [str(base), "-I", "-S", "-m", "venv", "--without-pip", str(destination)],
        [str(installer), "-I", "-m", "pip", "--isolated", "--disable-pip-version-check", "--python", python,
         "install", "--no-index", "--require-hashes", "--only-binary=:all:", "--no-cache-dir",
         "--find-links", str(wheels), "--report", str(destination.parent / (destination.name + "_install.json")),
         "-r", str(lock)],
        [python, "-I", "-m", "pip", "check"],
    ]


def score_worker(readers, directory):
    sys.path.insert(0, str(readers))
    from benchmark_tools.audit_publication_pipeline import audit
    from benchmark_tools.audit_installed_orthobench import read_root_hogs
    from benchmark_tools.run_installed_orthobench import fasta_ids
    from benchmark_tools.score_orthobench_partition import score_partition
    data = json.loads((directory / "data.json").read_text())
    validate_data(data)
    scientific = audit(directory / "native", directory / "readback")
    universe = fasta_ids(sorted((directory / "input").iterdir()))
    if len(universe) != data["genes"]:
        raise ValueError("Wrong input gene universe")
    groups = read_root_hogs(directory / "native/inference/orthohmm_phylogeny/orthohmm_root_hogs.tsv", universe)
    refs, uncertain = [{Path(r["path"]).name: set(Path(r["path"]).read_text().splitlines())
                       for r in data[role]} for role in ("references", "uncertain")]
    if not (set().union(*refs.values(), *uncertain.values()) - {""}) <= universe:
        raise ValueError("Unknown reference gene")
    score = score_partition(groups, refs, uncertain)
    validate_data(data)
    save(directory / "score.json", dict(dataset=data["dataset"], genes=len(universe), groups=len(groups),
        score=score, scientific_readback=record(directory / "readback/result.json"),
        summary=scientific["summary"], biological_validation=False, publication_ready=False))


def assembly_base(args, output):
    base = args.base_python.absolute()
    if (not base.is_file() or not os.access(base, os.X_OK)
            or record(base)["sha256"] != getattr(args, "base_python_sha256", None)
            or args.installer_python.absolute() != base or base.parent.name != "bin"
            or output.resolve() != output or output.is_relative_to(base.parent.parent)
            or output.is_relative_to(args.assets.resolve().parent)):
        raise ValueError("Assembly requires an externally pinned private base/installer and external canonical output")
    return record(base)


def validate_assembly_host(args):
    if not __package__:
        sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from benchmark_tools.check_publication_glibc_floor import check_host
    return check_host(args.assets.absolute().parent, args.abi_inventory,
                      args.abi_inventory_sha256, args.assembly_manifest_sha256)


def validate_base_probe(value, base):
    prefix = base.absolute().parent.parent
    site = Path(value["site"])
    if (value["version"] != "3.10.13" or value["implementation"] != "CPython"
            or value["machine"] != "x86_64" or value["prefix"] != str(prefix)
            or site.resolve() != site or not site.is_relative_to(prefix)
            or value["distributions"] != [["pip", "26.2.1"]]):
        raise ValueError("Assembly base runtime/site/distributions differ")


def run(args):
    output = args.output.absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    if args.cpu < 1 or args.timeout < 1:
        raise ValueError("CPU and timeout must be positive")
    data_record = record(args.data)
    if data_record["sha256"] != args.data_sha256:
        raise ValueError("Data manifest checksum mismatch")
    data = json.loads(args.data.read_text())
    validate_data(data)
    assembly_digest = getattr(args, "assembly_manifest_sha256", None)
    preflight_only = getattr(args, "preflight_only", False)
    abi_arguments = (getattr(args, "abi_inventory", None), getattr(args, "abi_inventory_sha256", None))
    if assembly_digest is not None and not all(abi_arguments):
        raise ValueError("Assembled execution requires externally anchored ABI inventory and glibc preflight")
    if assembly_digest is None and any(abi_arguments):
        raise ValueError("ABI inventory arguments require explicit assembled execution")
    if preflight_only and assembly_digest is None:
        raise ValueError("Preflight-only mode requires explicit assembled execution")
    validate_assets(args.assets, args.readers, assembly_digest)
    if assembly_digest is not None and (
            args.reader_wheels.resolve() != args.assets.resolve().parent / "reader_wheels"
            or args.reader_lock.resolve() != args.assets.resolve().parent / "reader_requirements.txt"):
        raise ValueError("Reader wheels/lock must belong to the externally anchored assembly")
    base_record = assembly_base(args, output) if assembly_digest is not None else None
    host_preflight = validate_assembly_host(args) if assembly_digest is not None else None
    locks = dict(inference=args.assets / "benchmark_tools/results/publication_recovery_requirements_20260926.txt",
                 reader=args.reader_lock)
    if any(record(p)["sha256"] != LOCKS[n] for n, p in locks.items()):
        raise ValueError("Changed environment lock")
    watched = [data_record, record(__file__), *[record(p) for p in locks.values()]]
    if base_record is not None:
        watched.append(base_record)
    if host_preflight is not None:
        watched.extend([host_preflight["inventory"], host_preflight["source"],
                        *[item["identity"] for item in host_preflight["loaders"]]])
    for root in (args.assets / "wheels", args.assets / "mafft", args.assets / "benchmark_tools",
                 args.reader_wheels, args.readers):
        watched.extend(record(p) for p in sorted(root.rglob("*")) if p.is_file())
    watched.append(record(args.assets / "FastTree"))
    if preflight_only:
        for row in watched:
            check(row)
        if validate_assembly_host(args) != host_preflight:
            raise ValueError("Assembly glibc preflight changed during preflight-only checks")
        validate_data(data)
        result = dict(status="integrated_assembled_preflight_complete", dataset=data["dataset"],
            inputs=watched, host_preflight=host_preflight, source=record(__file__),
            installation_executed=False, private_base_runtime_probe_executed=False,
            scientific_inference_executed=False, reviewed_native_code_executed=False,
            inference_execution_permitted=False, full_host_compatibility_verified=False,
            controlled_timing=False, publication_ready=False,
            limitations=["Pinned input/asset/base-file and glibc floor checks only; no native stage is launched.",
                "Private-base package/runtime probing and installed-payload/scientific admission remain separate.",
                "Not CPU/loader/library/OS/security closure, resource readiness or timing authorization."])
        output.mkdir(parents=True)
        save(output / "preflight.json", result)
        return result
    output.mkdir(parents=True)
    (output / "home").mkdir()
    (output / "input").mkdir()
    private_data = {k: v for k, v in data.items() if k not in ("fasta", "references", "uncertain")}
    for role in ("fasta", "references", "uncertain"):
        target_dir = output / ("input" if role == "fasta" else "scoring_inputs/" + role)
        target_dir.mkdir(parents=True, exist_ok=True)
        private_data[role] = []
        for row in data[role]:
            target = target_dir / Path(row["path"]).name
            shutil.copyfile(row["path"], target)
            if any(record(target)[k] != row[k] for k in ("bytes", "sha256")):
                raise ValueError("Input changed while copying")
            private_data[role].append(record(target))
    save(output / "data.json", private_data)
    environment = dict(HOME=str(output / "home"), PATH="/usr/bin:/bin", LANG="C.UTF-8",
                       OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", PYTHONHASHSEED="0")
    save(output / "started.json", dict(source=record(__file__), inputs=watched, cpu=args.cpu,
        dataset=data["dataset"], attempts=1, native_checkpoint_reuse=False, publication_ready=False))
    if host_preflight is not None:
        save(output / "glibc_preflight.json", host_preflight)
    outcomes = []
    try:
        before_base = None
        if assembly_digest is not None:
            from benchmark_tools.build_publication_project_wheel import RUNTIME_PROBE
            probe = [str(args.base_python), "-I", "-B", "-c", RUNTIME_PROBE]
            outcomes.append(stage(output, "base_before", probe, environment, 120))
            before_base = json.loads((output / "base_before.log").read_bytes())
            validate_base_probe(before_base, args.base_python)
        for name, wheels in (("inference", args.assets / "wheels"), ("reader", args.reader_wheels)):
            commands = install_commands(args.base_python, args.installer_python, wheels, locks[name], output / name)
            if assembly_digest is not None:
                commands = [[command[0], "-B", *command[1:]] for command in commands]
            for index, command in enumerate(commands):
                outcomes.append(stage(output, f"{name}_install_{index}", command, environment, 600))
        native_environment = dict(environment, MAFFT_BINARIES=str(args.assets / "mafft/libexec/mafft"))
        native = [str(output / "inference/bin/python"), "-I", "-B",
            str(args.assets / "benchmark_tools/run_publication_pipeline.py"), "--input", str(output / "input"),
            "--output", str(output / "native"), "--cpu", str(args.cpu),
            "--aligner", str(args.assets / "mafft/bin/mafft"), "--tree-builder", str(args.assets / "FastTree")]
        outcomes.append(stage(output, "native", native, native_environment, args.timeout))
        command = [str(output / "reader/bin/python"), "-I", "-B", str(Path(__file__).absolute()),
                   "--score-worker", str(args.readers), "--output", str(output)]
        outcomes.append(stage(output, "readback_score", command, environment, args.timeout))
        if assembly_digest is not None:
            outcomes.append(stage(output, "base_after", probe, environment, 120))
            after_base = json.loads((output / "base_after.log").read_bytes())
            validate_base_probe(after_base, args.base_python)
            if before_base != after_base:
                raise ValueError("Supplied private base site changed during assembled execution")
        for row in watched:
            check(row)
        validate_assets(args.assets, args.readers, assembly_digest)
        if host_preflight is not None and validate_assembly_host(args) != host_preflight:
            raise ValueError("Assembly glibc preflight changed during execution")
    except BaseException as error:
        save(output / "failure.json", dict(type=type(error).__name__, error=str(error), outcomes=outcomes, retry=False))
        raise
    result = dict(status="integrated_install_inference_readback_scoring_complete", dataset=data["dataset"],
        started=record(output / "started.json"), score=record(output / "score.json"), outcomes=outcomes,
        native_checkpoint_reuse=False, publication_ready=False, controlled_timing=False,
        limitations=["Local supplied assets and a trusted installer are prerequisites; acquisition is not automated.",
            "Installation uses pinned wheels but this controller does not independently audit all installed payload bytes.",
            "Same-host workflow integration, not new biological validation or controlled comparative timing."])
    if assembly_digest is not None:
        result.update(assembly_manifest_sha256=assembly_digest, base_unchanged=True,
                      glibc_preflight=record(output / "glibc_preflight.json"),
                      base_site_payload_files=len(before_base["files"]),
                      limitations=result["limitations"] + [
                          "New explicitly anchored assembly; historical asset manifests/admissions remain untouched.",
                          "Private base-site snapshot checked, not every base/OS path or transitive dependency."])
    save(output / "complete.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--score-worker", type=Path)
    for name in ("assets", "readers", "reader-wheels", "reader-lock", "data", "base-python", "installer-python"):
        parser.add_argument("--" + name, type=lambda p: Path(p).absolute())
    parser.add_argument("--data-sha256")
    parser.add_argument("--assembly-manifest-sha256")
    parser.add_argument("--base-python-sha256")
    parser.add_argument("--abi-inventory", type=lambda p: Path(p).absolute())
    parser.add_argument("--abi-inventory-sha256")
    parser.add_argument("--preflight-only", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cpu", type=int, default=2)
    parser.add_argument("--timeout", type=int, default=86400)
    args = parser.parse_args()
    if args.score_worker:
        score_worker(args.score_worker.resolve(), args.output.resolve())
    else:
        if any(getattr(args, n) is None for n in
               ("assets", "readers", "reader_wheels", "reader_lock", "data", "data_sha256", "base_python", "installer_python")):
            parser.error("Require all asset, data, digest and interpreter arguments")
        print(run(args)["status"])
