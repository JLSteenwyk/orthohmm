"""Acquire a pinned Miniforge bootstrap and install it in a fresh private prefix."""

import argparse
import json
import os
from pathlib import Path
import platform
import sys
from urllib.parse import urlsplit
from urllib.request import build_opener, HTTPRedirectHandler, Request

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmark_tools.run_integrated_publication_workflow import record, save, stage

NAME = "Miniforge3-25.3.1-0-Linux-x86_64.sh"
URL = "https://github.com/conda-forge/miniforge/releases/download/25.3.1-0/" + NAME
EXPECTED = dict(bytes=93870801, sha256="376b160ed8130820db0ab0f3826ac1fc85923647f75c1b8231166e3d559ab768")
CHECKSUM = dict(bytes=104, sha256="57be9d8415cd75326aff2a518bd6eda8d3a87ef28747256d51a8b77a47cc0a14")


def provider(url):
    parsed = urlsplit(url)
    if (parsed.scheme != "https" or parsed.hostname not in {"github.com", "release-assets.githubusercontent.com"}
            or parsed.username or parsed.password or parsed.port not in {None, 443} or parsed.fragment
            or parsed.hostname == "github.com" and parsed.query or any(c.isspace() for c in url)):
        raise ValueError("Unsafe bootstrap provider/redirect")
    return parsed.hostname


class BootstrapRedirect(HTTPRedirectHandler):
    def redirect_request(self, request, response, code, message, headers, url):
        provider(url)
        return super().redirect_request(request, response, code, message, headers, url)


def fresh(args):
    output = args.output.absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    if (output.resolve() != output or type(args.timeout) is not int or args.timeout < 1
            or args.acknowledge_historical_runtime is not True):
        raise ValueError("Require fresh canonical output, positive timeout and historical acknowledgement")
    return output


def identity(path):
    ref = record(path)
    return {key: ref[key] for key in ("bytes", "sha256")}


def download(opener, url, target, expected, timeout):
    provider(url)
    partial = target.with_name(target.name + ".partial")
    with opener.open(Request(url, headers={"User-Agent": "OrthoHMM-publication-bootstrap",
        "Accept-Encoding": "identity"}), timeout=timeout) as response, partial.open("xb") as stream:
        host = provider(response.geturl())
        if response.headers.get("Content-Encoding", "identity") != "identity":
            raise ValueError("Encoded bootstrap response")
        length = response.headers.get("Content-Length")
        if length is not None and int(length) != expected["bytes"]:
            raise ValueError("Bootstrap response length differs")
        size = 0
        while True:
            block = response.read(min(1024 ** 2, expected["bytes"] - size + 1))
            if not block:
                break
            if len(block) > expected["bytes"] - size:
                raise ValueError("Bootstrap response exceeds bound")
            size += len(block)
            stream.write(block)
    if identity(partial) != expected:
        raise ValueError("Bootstrap artifact differs from fixed provider identity")
    if target.exists() or target.is_symlink():
        raise FileExistsError(target)
    partial.rename(target)
    target.chmod(0o644)
    # GitHub signed redirect queries are temporary credentials; never retain them.
    return dict(provider_url=url, final_host=host, file=record(target))


def acquire(args):
    output = fresh(args)
    source = record(__file__)
    output.mkdir(parents=True)
    save(output / "started.json", dict(source=source, installer=EXPECTED, checksum=CHECKSUM, attempts=1))
    downloads = []
    try:
        opener = build_opener(BootstrapRedirect())
        for suffix, expected in ((".sha256", CHECKSUM), ("", EXPECTED)):
            downloads.append(download(opener, URL + suffix, output / (NAME + suffix), expected, args.timeout))
        if record(__file__) != source:
            raise ValueError("Bootstrap acquisition source changed")
        result = dict(status="pinned_bootstrap_acquired", source=source, downloads=downloads, attempts=1,
            retry=False, installation_performed=False, publication_ready=False, security_clearance=False,
            redistribution_clearance=False,
            limitations=["Pinned historical official release, not a current-security recommendation or binary-signature proof.",
                "No package installation, scientific execution or OS/dependency/source-rights closure."])
        save(output / "complete.json", result)
        return result
    except BaseException as error:
        save(output / "failed.json", dict(status="bootstrap_acquisition_failed", type=type(error).__name__,
            error=str(error), downloads=downloads, attempts=1, retry=False, publication_ready=False))
        raise


def install(args):
    output = fresh(args)
    if platform.system() != "Linux" or platform.machine() != "x86_64":
        raise ValueError("Pinned bootstrap requires Linux x86-64")
    installer = args.installer.absolute()
    if (installer.is_symlink() or not installer.is_file() or installer.name != NAME
            or installer.stat().st_size != EXPECTED["bytes"] or identity(installer) != EXPECTED):
        raise ValueError("Require the exact regular supplied bootstrap installer")
    watched = [record(installer), record(__file__), record(sys.modules[stage.__module__].__file__)]
    output.mkdir(parents=True)
    for name in ("home", "cache", "envs"):
        (output / name).mkdir()
    prefix = output / "prefix"
    environment = dict(HOME=str(output / "home"), PATH=str(prefix / "bin") + ":/usr/bin:/bin", LANG="C.UTF-8",
        CONDARC="/dev/null", CONDA_NO_PLUGINS="true", CONDA_PKGS_DIRS=str(output / "cache"),
        CONDA_ENVS_PATH=str(output / "envs"), CONDA_EXTRACT_THREADS="1", CONDA_VERIFY_THREADS="1",
        PYTHONDONTWRITEBYTECODE="1", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    save(output / "started.json", dict(inputs=watched, prefix=str(prefix), environment=environment,
        attempts=1, shell_init_requested=False, historical_runtime=True))
    outcomes = []
    try:
        outcomes.append(stage(output, "install", ["/bin/bash", str(installer), "-b", "-p", str(prefix)], environment, args.timeout))
        conda, python = prefix / "bin/conda", prefix / "bin/python"
        if (not conda.is_file() or not os.access(conda, os.X_OK) or not conda.resolve().is_relative_to(prefix)
                or not python.is_file() or not python.resolve().is_relative_to(prefix)):
            raise ValueError("Installed bootstrap entrypoints missing/escaping")
        outcomes.append(stage(output, "conda_version", [str(conda), "--version"], environment, args.timeout))
        if (output / "conda_version.log").read_text().strip() != "conda 25.3.1":
            raise ValueError("Installed bootstrap Conda version differs")
        packages = []
        for path in sorted((prefix / "conda-meta").glob("*.json")):
            if path.is_symlink():
                raise ValueError("Symlinked bootstrap metadata")
            data = json.loads(path.read_bytes())
            packages.append(dict(name=data["name"], version=data["version"], build=data["build"], metadata=record(path)))
        if (len({p["name"] for p in packages}) != len(packages)
                or [p["version"] for p in packages if p["name"] == "conda"] != ["25.3.1"]):
            raise ValueError("Installed bootstrap package inventory differs")
        for item in watched:
            if record(item["path"]) != item:
                raise ValueError("Supplied bootstrap/source changed")
        result = dict(status="pinned_private_bootstrap_installed", inputs=watched, outcomes=outcomes,
            conda=record(conda), python=record(python), packages=packages, environment=environment,
            attempts=1, retry=False, shell_init_requested=False, scientific_environment_installed=False,
            native_inference_executed=False, publication_ready=False, security_clearance=False,
            redistribution_clearance=False,
            limitations=["Fresh fixed Miniforge bootstrap, not the frozen scientific base or its native admission.",
                "Package metadata/entrypoint identity, not complete installed payload, bootstrap dependency, OS or rights closure.",
                "Batch mode and private HOME/config/cache; no shell initialization or shared installation requested."])
        save(output / "complete.json", result)
        return result
    except BaseException as error:
        save(output / "failed.json", dict(status="bootstrap_install_failed", type=type(error).__name__,
            error=str(error), outcomes=outcomes, attempts=1, retry=False, publication_ready=False))
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("acquire", "install"):
        cmd = commands.add_parser(name)
        cmd.add_argument("--output", type=Path, required=True)
        cmd.add_argument("--timeout", type=int, default=600 if name == "install" else 60)
        cmd.add_argument("--acknowledge-historical-runtime", action="store_true")
        if name == "install":
            cmd.add_argument("--installer", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps((acquire if args.command == "acquire" else install)(args), indent=2, sort_keys=True))
