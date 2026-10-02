"""Reconstruct the historical publication base offline in a fresh private prefix."""

import argparse
from importlib import metadata
import json
import os
from pathlib import Path
import platform
import shutil
import sys
import sysconfig

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmark_tools import stage_base_archives as archives
from benchmark_tools.audit_frozen_overlay_install import local_install_wheels
from benchmark_tools.audit_recovery_install import installed_payload
from benchmark_tools.run_integrated_publication_workflow import check, record, save, stage

RECEIPT_SHA = "7af55a78192eea22cf978ff331c1c5556965e17862b565793321b14857479f59"
PIP_NAME = "pip-26.2.1-py3-none-any.whl"
PIP_SHA = "71138adf1f4ca900cdb7d289c21b7494329f2332b6d85f0e1c42108c0384ed3e"
PIP_BYTES = 1816632


def supported_host():
    return platform.system() == "Linux" and platform.machine() == "x86_64"


def regular(path):
    path = Path(path).absolute()
    if path.is_symlink() or not path.is_file():
        raise ValueError("Require a regular supplied input: " + str(path))
    return path


def preflight(args):
    output = args.output.absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    if output.resolve() != output or type(args.timeout) is not int or args.timeout < 1:
        raise ValueError("Require a canonical fresh destination and positive timeout")
    if args.acknowledge_historical_runtime is not True:
        raise ValueError("Explicit historical-runtime acknowledgement is required")
    if not supported_host():
        raise ValueError("Historical base installation requires Linux x86-64")
    receipt, conda, wheel = [regular(p) for p in (args.receipt, args.conda, args.pip_wheel)]
    if record(receipt)["sha256"] != RECEIPT_SHA:
        raise ValueError("Changed frozen reconstruction receipt")
    if record(conda)["sha256"] != args.conda_sha256 or not os.access(conda, os.X_OK):
        raise ValueError("Supplied Conda entrypoint pin or executable permission differs")
    pip = record(wheel)
    if wheel.name != PIP_NAME or pip["sha256"] != PIP_SHA or pip["bytes"] != PIP_BYTES:
        raise ValueError("Changed pip bootstrap wheel")
    packages = archives.packages_from_receipt(json.loads(receipt.read_bytes()))
    if len(packages) != 19:
        raise ValueError("Require the complete retained 19-package base")
    watched = [record(p) for p in (receipt, conda, wheel, Path(__file__),
        Path(archives.__file__), Path(sys.modules[stage.__module__].__file__),
        Path(sys.modules[installed_payload.__module__].__file__))]
    for item in packages:
        path = regular(args.cache / item["package"]["fn"])
        if archives.identity(path) != item["expected"]:
            raise ValueError("Changed supplied base archive")
        watched.append(record(path))
    return output, packages, watched


def conda_inventory(prefix, packages):
    expected = {(r["package"]["name"], r["package"]["version"], r["package"]["build"])
                for r in packages}
    rows, observed = [], set()
    for path in sorted((prefix / "conda-meta").glob("*.json")):
        regular(path)
        item = json.loads(path.read_bytes())
        key = (item["name"], item["version"], item["build"])
        if key in observed:
            raise ValueError("Duplicate installed Conda package")
        observed.add(key)
        rows.append(dict(name=key[0], version=key[1], build=key[2], metadata=record(path)))
    if observed != expected or len(rows) != 19:
        raise ValueError("Installed Conda package inventory differs")
    return rows


def snapshot(output):
    site = Path(sysconfig.get_paths()["purelib"]).resolve()
    distributions = []
    for distribution in metadata.distributions(path=[str(site)]):
        if Path(distribution.locate_file("")).resolve() != site:
            raise ValueError("Distribution metadata escapes the installed site")
        distributions.append(dict(name=distribution.metadata["Name"], version=distribution.version))
    save(output, dict(python=record(sys.executable), prefix=sys.prefix, site=str(site),
        version=platform.python_version(), implementation=platform.python_implementation(),
        machine=platform.machine(), distributions=distributions))


def validate_snapshot(value, prefix):
    site = Path(value["site"])
    if (value["prefix"] != str(prefix) or value["version"] != "3.10.13"
            or value["implementation"] != "CPython" or value["machine"] != "x86_64"
            or value["python"]["path"] != str(prefix / "bin/python")
            or not site.is_relative_to(prefix) or site.resolve() != site
            or value["distributions"] != [dict(name="pip", version="26.2.1")]):
        raise ValueError("Installed Python/pip runtime differs from retained base")
    regular(prefix / "bin/python")
    check(value["python"])
    return site


def commands(conda, output, timeout):
    prefix, staging = output / "python-runtime", output / "staging"
    python, wheel = str(prefix / "bin/python"), output / "bootstrap_wheels" / PIP_NAME
    bootstrap = "import sys,runpy;sys.path.insert(0," + repr(str(wheel)) + ");runpy.run_module('pip',run_name='__main__')"
    return [
        ("conda_install", [str(conda), "create", "--offline", "--copy", "--no-default-packages",
            "--solver", "classic", "--yes", "--json", "--prefix", str(prefix),
            "--file", str(staging / "explicit.txt")], timeout),
        ("pip_bootstrap", [python, "-I", "-S", "-B", "-c", bootstrap, "--isolated",
            "--disable-pip-version-check", "install", "--prefix", str(prefix), "--no-index", "--no-deps",
            "--require-hashes", "--only-binary=:all:", "--no-cache-dir", "--find-links",
            str(wheel.parent), "--report", str(output / "bootstrap_install.json"),
            "-r", str(output / "bootstrap_pip.txt")], timeout),
        ("pip_check", [python, "-I", "-B", "-m", "pip", "check"], timeout),
        ("runtime_snapshot", [python, "-I", "-B", str(Path(__file__).absolute()), "snapshot",
            "--output", str(output / "runtime.json")], timeout),
    ]


def run(args):
    output, packages, watched = preflight(args)
    output.mkdir(parents=True, exist_ok=False)
    environment = dict(HOME=str(output / "home"), PATH="/usr/bin:/bin", LANG="C.UTF-8",
        CONDA_NO_PLUGINS="true", CONDA_PKGS_DIRS=str(output / "cache"),
        CONDA_ENVS_PATH=str(output / "envs"), CONDARC="/dev/null", PIP_CONFIG_FILE="/dev/null",
        PYTHONDONTWRITEBYTECODE="1", PYTHONNOUSERSITE="1", PYTHONHASHSEED="0",
        CONDA_EXTRACT_THREADS="1", CONDA_VERIFY_THREADS="1")
    (output / "home").mkdir()
    save(output / "started.json", dict(inputs=watched, environment=environment, attempts=1,
        historical_runtime=True, security_clearance=False, publication_ready=False))
    outcomes = []
    try:
        staged = archives.stage(args.receipt, args.cache, output / "staging")
        wheels = output / "bootstrap_wheels"
        wheels.mkdir()
        wheel = wheels / PIP_NAME
        shutil.copyfile(args.pip_wheel, wheel)
        if record(wheel)["sha256"] != PIP_SHA or record(wheel)["bytes"] != PIP_BYTES:
            raise ValueError("Bootstrap wheel changed while copying")
        with (output / "bootstrap_pip.txt").open("x") as stream:
            stream.write("pip==26.2.1 --hash=sha256:" + PIP_SHA + "\n")
        staged_inputs = [record(output / "staging" / item["archive"]["path"])
                         for item in staged["packages"]]
        staged_inputs += [record(output / "staging/explicit.txt"), record(wheel),
                          record(output / "bootstrap_pip.txt")]
        prefix = output / "python-runtime"
        for name, command, timeout in commands(args.conda.absolute(), output, args.timeout):
            for item in watched + staged_inputs:
                check(item)
            outcomes.append(stage(output, name, command, environment, timeout))
            if name == "conda_install":
                conda_inventory(prefix, packages)
        installed = conda_inventory(prefix, packages)
        runtime = json.loads((output / "runtime.json").read_bytes())
        site = validate_snapshot(runtime, prefix)
        report = json.loads((output / "bootstrap_install.json").read_bytes())
        installed_wheels = local_install_wheels(report, wheels)
        if (len(installed_wheels) != 1 or installed_wheels[0]["name"] != "pip"
                or installed_wheels[0]["version"] != "26.2.1"
                or installed_wheels[0]["wheel"]["sha256"] != PIP_SHA):
            raise ValueError("Bootstrap installation report differs")
        payload = installed_payload(wheel, site)
        for item in watched + staged_inputs:
            check(item)
        for item in staged["packages"]:
            if archives.identity(output / "staging" / item["archive"]["path"]) != {
                    k: item["archive"][k] for k in ("bytes", "sha256", "md5")}:
                raise ValueError("Staged base archive changed during installation")
        if record(wheel)["sha256"] != PIP_SHA or record(wheel)["bytes"] != PIP_BYTES:
            raise ValueError("Staged pip wheel changed during installation")
        result = dict(status="historical_publication_base_installed", inputs=watched, stages=outcomes,
            staging=record(output / "staging/staging.json"), installed_conda_packages=installed,
            runtime=runtime, bootstrap_report=record(output / "bootstrap_install.json"),
            pip_payload=payload, historical_runtime=True, native_inference_executed=False,
            controlled_timing=False, publication_ready=False, security_clearance=False,
            redistribution_cleared=False, limitations=[
                "Historical reproduction base, not a patched general-purpose installation recommendation.",
                "Caller supplies the trusted Conda bootstrap; its entrypoint pin is not its full dependency closure.",
                "19 Conda metadata records are checked, not every transformed base payload or OS library.",
                "Pip payload exclusions follow the existing audit: generated RECORD/bytecode and non-site data.",
                "Native tools, scientific wheel sets, inference and cross-host validation remain separate."])
        save(output / "complete.json", result)
        return result
    except Exception as error:
        save(output / "failed.json", dict(status="publication_base_installation_failed", attempts=1,
            retry=False, error_type=type(error).__name__, error=str(error), completed_stages=outcomes,
            publication_ready=False))
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    subcommands = parser.add_subparsers(dest="command", required=True)
    installer = subcommands.add_parser("install")
    for name in ("receipt", "cache", "conda", "pip-wheel", "output"):
        installer.add_argument("--" + name, type=Path, required=True)
    installer.add_argument("--conda-sha256", required=True)
    installer.add_argument("--timeout", type=int, default=600)
    installer.add_argument("--acknowledge-historical-runtime", action="store_true")
    worker = subcommands.add_parser("snapshot")
    worker.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "snapshot":
        snapshot(args.output)
    else:
        result = run(args)
        print(json.dumps(dict(status=result["status"], complete=record(args.output / "complete.json")), indent=2))
