"""Privately build the pinned double-precision FastTree source for baseline x86-64."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess

from benchmark_tools.acquire_publication_fasttree import FILES, identity, verify
from benchmark_tools.build_publication_mafft import execute


def compile_command():
    return ["/usr/bin/gcc", "-Wall", "-O3", "-finline-functions", "-funroll-loops",
            "-march=x86-64", "-mtune=generic", "-o", "FastTree", "FastTree.c", "-lm"]


def run(repo, source, mafft, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    originals = [verify(source / name, FILES[name]) for name in ("FastTree.c", "LICENSE")]
    output.mkdir(parents=True)
    env = dict(PATH="/usr/bin:/bin", LC_ALL="C", OMP_NUM_THREADS="1",
               OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    report = dict(status="building", sources=originals, script=identity(__file__),
        compile_command=compile_command(), compile_environment=env, builds=[], commands=[],
        historical_binary_modified=False, scientific_defaults_changed=False,
        cross_platform_verified=False, publication_ready=False)
    try:
        report["compiler"] = dict(binary=identity(Path("/usr/bin/gcc").resolve()),
            version=subprocess.check_output(["/usr/bin/gcc", "--version"], text=True, env=env))
        report["linker"] = identity(Path("/usr/bin/ld").resolve())
        for name in ("build_a", "build_b"):
            build = output / name
            build.mkdir()
            for item in originals:
                path = Path(item["path"])
                shutil.copyfile(path, build / path.name)
                verify(build / path.name, item["sha256"])
            command = compile_command()
            report["commands"].append(dict(argv=command, cwd=str(build), timeout_seconds=240))
            execute(command, build, env, output / f"{name}.log", timeout=240)
            report["builds"].append(identity(build / "FastTree"))
        report["repeat_build_byte_equal"] = report["builds"][0]["sha256"] == report["builds"][1]["sha256"]
        binary = output / "build_a/FastTree"
        probe = subprocess.run([str(binary), "-help"], env=env, capture_output=True, text=True, timeout=30)
        report["help_probe"] = dict(returncode=probe.returncode, stdout=probe.stdout, stderr=probe.stderr)
        if probe.returncode != 0 or "FastTree Version 2.2.0 Double precision" not in probe.stdout + probe.stderr:
            raise ValueError("Unexpected built FastTree identity/precision")
        for label, command in (("elf", ["/usr/bin/readelf", "-h", "-n", "-d", str(binary)]),
                               ("runtime_libraries", ["/usr/bin/ldd", str(binary)])):
            report[label] = subprocess.check_output(command, text=True, env=env, timeout=30)
        fixture = output / "phylogeny"
        command = ["/usr/bin/python3", "-S", "-m", "benchmark_tools.verify_frozen_phylogeny_install",
            "--repo", str(repo), "--artifact", str(repo / "benchmarks/work/publication_frozen_overlay_20260926"),
            "--mafft", str(mafft), "--fasttree", str(binary), "--output", str(fixture)]
        fixture_env = dict(env, MAFFT_BINARIES=str(mafft.parents[1] / "libexec/mafft"))
        report["fixture_environment"] = fixture_env
        report["commands"].append(dict(argv=command, cwd=str(repo), timeout_seconds=240))
        execute(command, repo, fixture_env, output / "phylogeny_driver.log", timeout=240)
        report["fixture_report"] = identity(fixture / "report.json")
        baseline = repo / "benchmarks/work/publication_frozen_phylogeny_smoke_20260926/inference"
        current = fixture / "inference"
        names = ["orthohmm_orthogroups.txt", "orthohmm_phylogeny/orthohmm_root_hogs.tsv",
                 "orthohmm_phylogeny/orthohmm_pairwise_orthologs.tsv",
                 "orthohmm_phylogeny/species_tree.rooted.nwk",
                 "orthohmm_phylogeny/gene_trees/Family0000002.reconciled.nwk"]
        report["fixture_comparisons"] = [dict(relative=name, baseline=identity(baseline / name),
            current=identity(current / name), byte_equal=(baseline / name).read_bytes() == (current / name).read_bytes())
            for name in names]
        for item in originals:
            verify(Path(item["path"]), item["sha256"])
        for item in report["builds"]:
            verify(Path(item["path"]), item["sha256"])
        report["status"] = "baseline_x86_64_build_and_installed_fixture_completed"
    except Exception as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        report["logs"] = [identity(p) for p in sorted(output.glob("*.log"))]
        with (output / "report.json").open("x") as stream:
            json.dump(report, stream, indent=2, sort_keys=True)
            stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--mafft", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.repo.resolve(), args.source.resolve(), args.mafft.resolve(), args.output.absolute())
