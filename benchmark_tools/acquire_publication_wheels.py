"""Acquire the exact historical inference/reader wheel sets without installation."""

import argparse
import json
from pathlib import Path
import re
import shutil
import sys
from urllib.parse import quote, unquote, urlsplit

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmark_tools import acquire_publication_base as network
from benchmark_tools.run_integrated_publication_workflow import record, save

INVENTORY_SHA = "97eaa8cca56db98b337ad161e619006ef7a91433d2c704e9dcc16ceacc5e7dab"
LOCK_SHAS = dict(inference="813751262c9405155cc2f1974a2c32187b352f88307e2332e0cec9df35708ae3",
                 reader="53df1cabd91c2fb18179951ebc95166ae49ba843f4015e0f417468d2a8265063")
ROLES = dict(inference={"dendropy", "igraph", "leidenalg", "llvmlite", "numba", "numpy",
                       "orthohmm", "pip", "python-igraph", "setuptools", "texttable"},
             reader={"biopython", "dendropy", "numpy", "pip", "setuptools"})


def normalized(name):
    return re.sub(r"[-_.]+", "-", name).lower()


def pinned(path, digest, maximum):
    path = Path(path).absolute()
    if path.is_symlink() or not path.is_file() or path.stat().st_size > maximum:
        raise ValueError("Require a bounded regular supplied file")
    ref = record(path)
    if ref["sha256"] != digest:
        raise ValueError("Supplied input differs from frozen identity")
    return ref


def wheel_rows(inventory):
    if inventory["status"] != "selected_wheel_elf_inventory":
        raise ValueError("Wrong admitted wheel inventory")
    rows = {}
    for item in inventory["wheels"]:
        name, wheel = normalized(item["name"]), item["wheel"]
        filename = Path(wheel["path"]).name
        if (name in rows or not re.fullmatch(r"[a-z0-9][a-z0-9-]*", name)
                or not re.fullmatch(r"[0-9]+(?:\.[0-9]+)+", item["version"])
                or not re.fullmatch(r"[A-Za-z0-9_.+-]+\.whl", filename)
                or type(wheel["bytes"]) is not int or not 0 < wheel["bytes"] <= 100_000_000
                or not re.fullmatch(r"[0-9a-f]{64}", wheel["sha256"])):
            raise ValueError("Invalid/duplicate retained wheel identity")
        rows[name] = dict(name=name, version=item["version"], filename=filename,
                          expected={k: wheel[k] for k in ("bytes", "sha256")})
    if set(rows) != set.union(*ROLES.values()):
        raise ValueError("Require exactly the 12 retained inference/reader artifacts")
    return rows


def public_selection(metadata, wheel):
    if (normalized(metadata["info"]["name"]) != wheel["name"]
            or metadata["info"]["version"] != wheel["version"]):
        raise ValueError("Wrong release-specific wheel metadata")
    matches = [row for row in metadata["urls"] if row["filename"] == wheel["filename"]]
    if len(matches) != 1:
        raise ValueError("Require exactly the retained platform wheel")
    row = matches[0]
    if (row["packagetype"] != "bdist_wheel" or row["yanked"] is not False
            or type(row["size"]) is not int or row["size"] != wheel["expected"]["bytes"]
            or row["digests"]["sha256"] != wheel["expected"]["sha256"]):
        raise ValueError("Published artifact differs from retained wheel")
    url = network.provider_url(row["url"], {"files.pythonhosted.org"})
    if unquote(urlsplit(url).path.rsplit("/", 1)[-1]) != wheel["filename"]:
        raise ValueError("Published wheel basename differs")
    return url


def copy_wheel(source, target, expected):
    if target.exists() or target.is_symlink():
        raise FileExistsError(target)
    with Path(source).open("rb") as src, target.open("xb") as dst:
        shutil.copyfileobj(src, dst, 1024 * 1024)
    target.chmod(0o644)
    ref = record(target)
    if any(ref[key] != value for key, value in expected.items()):
        raise ValueError("Copied wheel differs from retained identity")
    return ref


def acquire(args):
    output = args.output.absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    if (output.resolve() != output or type(args.timeout) is not int or args.timeout < 1
            or args.acknowledge_historical_runtime is not True):
        raise ValueError("Require fresh canonical output, positive timeout and historical acknowledgement")
    inputs = [pinned(args.inventory, INVENTORY_SHA, 1024 ** 2),
              *[pinned(getattr(args, role + "_lock"), digest, 8192) for role, digest in LOCK_SHAS.items()]]
    wheels = wheel_rows(json.loads(args.inventory.read_bytes()))
    supplied = dict(orthohmm=args.orthohmm_wheel, pip=args.pip_wheel)
    for name, path in supplied.items():
        ref = pinned(path, wheels[name]["expected"]["sha256"], 100_000_000)
        if path.name != wheels[name]["filename"] or ref["bytes"] != wheels[name]["expected"]["bytes"]:
            raise ValueError("Wrong supplied local wheel size/name")
        inputs.append(ref)
    inputs.extend(record(path) for path in (Path(__file__), Path(network.__file__)))
    output.mkdir(parents=True, exist_ok=False)
    for name in ("wheels", "metadata", "inference_wheels", "reader_wheels"):
        (output / name).mkdir()
    save(output / "started.json", dict(inputs=inputs, wheels=wheels, roles={k: sorted(v) for k, v in ROLES.items()},
        attempts=1, timeout=args.timeout, historical_runtime=True))
    acquisitions = []
    try:
        opener = network.network_opener()
        for name, wheel in sorted(wheels.items()):
            target = output / "wheels" / wheel["filename"]
            if name in supplied:
                ref = copy_wheel(supplied[name], target, wheel["expected"])
                acquisitions.append(dict(name=name, origin="supplied_exact_artifact", wheel=ref))
                continue
            url = "https://pypi.org/pypi/{}/{}/json".format(quote(name, safe=""), quote(wheel["version"], safe=""))
            metadata_path = output / "metadata" / (name + ".json")
            metadata = network.download(opener, url, metadata_path, network.METADATA_LIMIT, args.timeout)
            selected_url = public_selection(json.loads(metadata_path.read_bytes()), wheel)
            downloaded = network.download(opener, selected_url, target, wheel["expected"]["bytes"],
                                          args.timeout, wheel["expected"])
            acquisitions.append(dict(name=name, origin="public_provider", wheel=downloaded["file"],
                                     metadata=metadata, download=downloaded))
        environments = {}
        for role, names in ROLES.items():
            members = [copy_wheel(output / "wheels" / wheels[name]["filename"],
                output / (role + "_wheels") / wheels[name]["filename"], wheels[name]["expected"])
                for name in sorted(names)]
            lock = output / (role + "_requirements.txt")
            with lock.open("xb") as stream:
                stream.write(getattr(args, role + "_lock").read_bytes())
            lock.chmod(0o644)
            if record(lock)["sha256"] != LOCK_SHAS[role]:
                raise ValueError("Copied exact lock differs")
            environments[role] = dict(wheels=members, lock=record(lock))
        for item in inputs:
            if record(item["path"]) != item:
                raise ValueError("Supplied input/source changed during acquisition")
        for row in acquisitions:
            if record(row["wheel"]["path"]) != row["wheel"]:
                raise ValueError("Acquired wheel changed")
            if "metadata" in row and record(row["metadata"]["file"]["path"]) != row["metadata"]["file"]:
                raise ValueError("Provider metadata changed")
        for environment in environments.values():
            for ref in [*environment["wheels"], environment["lock"]]:
                if record(ref["path"]) != ref:
                    raise ValueError("Prepared environment artifact changed")
        result = dict(status="historical_inference_reader_wheels_acquired", inputs=inputs, acquisitions=acquisitions,
            environments=environments, union_wheels=12, public_wheels=10, supplied_wheels=2,
            union_artifact_bytes=sum(w["expected"]["bytes"] for w in wheels.values()),
            attempts=1, retry=False, installation_performed=False, native_inference_executed=False,
            publication_ready=False, security_clearance=False, redistribution_clearance=False,
            limitations=["The unpublished OrthoHMM wheel remains supplied, not publicly acquired or rebuilt here.",
                "Pip is a supplied exact artifact; public provider queries never substitute versions, wheels or platforms.",
                "Opaque exact historical lock bytes and admitted artifact role sets, not new dependency resolution or wheel-tag admission.",
                "No imports of downloaded code, installation, inference, controlled timing, source/OS closure or rights clearance."])
        save(output / "complete.json", result)
        return result
    except BaseException as error:
        save(output / "failed.json", dict(status="wheel_acquisition_failed", inputs=inputs, acquisitions=acquisitions,
            type=type(error).__name__, error=str(error), attempts=1, retry=False, publication_ready=False))
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("inventory", "inference-lock", "reader-lock", "orthohmm-wheel", "pip-wheel", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--timeout", type=int, default=60)
    parser.add_argument("--acknowledge-historical-runtime", action="store_true")
    print(json.dumps(acquire(parser.parse_args()), indent=2, sort_keys=True))
