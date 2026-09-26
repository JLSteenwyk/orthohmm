"""Build pinned MAFFT core privately; retain source and execution provenance."""

import argparse
import json
import os
from pathlib import Path, PurePosixPath
import signal
import subprocess
import tarfile

from benchmark_tools.acquire_publication_fasttree import identity, verify

URL = "https://mafft.cbrc.jp/alignment/software/mafft-7.525-with-extensions-src.tgz"
SHA = "2876f4adc1a2de4ed206bc40896763bf208bf1a02bda52f8bfdd91cf52d73e4a"
ROOT = "mafft-7.525-with-extensions"


def unpack(archive, destination):
    verify(archive, SHA)
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(destination)
    with tarfile.open(archive, "r:gz") as stream:
        members = stream.getmembers()
        names = [m.name.rstrip("/") for m in members]
        if len(names) != len(set(names)) or sum(m.size for m in members) > 20_000_000:
            raise ValueError("Duplicate or oversized archive")
        for m in members:
            p = PurePosixPath(m.name)
            if (p.is_absolute() or ".." in p.parts or not p.parts or p.parts[0] != ROOT
                    or not (m.isfile() or m.isdir())):
                raise ValueError("Unsafe archive member")
        destination.mkdir(parents=True)
        records = []
        for m in members:
            target = destination / m.name
            if m.isdir():
                target.mkdir(parents=True, exist_ok=True)
            else:
                target.parent.mkdir(parents=True, exist_ok=True)
                with stream.extractfile(m) as source, target.open("xb") as sink:
                    sink.write(source.read())
                target.chmod(m.mode & 0o777)
                records.append(identity(target))
    return records


def execute(command, cwd, env, log, timeout=600):
    with log.open("x") as stream:
        child = subprocess.Popen(command, cwd=cwd, env=env, stdout=stream,
                                 stderr=subprocess.STDOUT, start_new_session=True)
        try:
            code = child.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            os.killpg(child.pid, signal.SIGKILL)
            child.wait()
            raise
    if code:
        raise RuntimeError(f"Command exited {code}: {command}")


def run(repo, archive, installed, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    archive_record = verify(archive, SHA)
    output.mkdir(parents=True)
    report = dict(status="building", archive=archive_record, acquisition_url=URL,
                  script=identity(__file__), commands=[], extensions_built=False,
                  historical_installation_modified=False, publication_ready=False)
    try:
        sources = unpack(archive, output / "source")
        report["source_files"] = sources
        comparisons = []
        for source in sources:
            rel = Path(source["path"]).relative_to(output / "source" / ROOT)
            previous = identity(installed / rel)
            comparisons.append(dict(relative=str(rel), retained=previous,
                                    archive_sha256=source["sha256"],
                                    byte_equal=previous["sha256"] == source["sha256"]))
        report["retained_source_comparisons"] = comparisons
        prefix = output / "install"
        env = dict(os.environ, PATH="/usr/bin:/bin", MAKEFLAGS="", MFLAGS="",
                   OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
        env.pop("MAFFT_BINARIES", None)
        report["environment_overrides"] = {key: env[key] for key in
            ("PATH", "MAKEFLAGS", "MFLAGS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}
        report["environment_unset"] = ["MAFFT_BINARIES"]
        report["compiler"] = dict(file=identity(Path("/usr/bin/gcc").resolve()), version=subprocess.check_output(
            ["/usr/bin/gcc", "--version"], text=True, env=env))
        core = output / "source" / ROOT / "core"
        command = ["/usr/bin/make", "-j2", "CC=/usr/bin/gcc", "CFLAGS=-O3",
                   f"PREFIX={prefix}", "install"]
        report["commands"].append(dict(argv=command, cwd=str(core), timeout_seconds=600))
        execute(command, core, env, output / "build.log")
        helpers = sorted((prefix / "libexec/mafft").iterdir())
        report["built_helpers"] = [identity(p) for p in helpers]
        retained_helpers = sorted((installed / "libexec/mafft").iterdir())
        report["retained_helpers"] = [identity(p) for p in retained_helpers]
        report["same_helper_names"] = [p.name for p in helpers] == [p.name for p in retained_helpers]
        report["launcher"] = identity(prefix / "bin/mafft")
        phylogeny = output / "phylogeny"
        command = ["/usr/bin/python3", "-S", "-m", "benchmark_tools.verify_frozen_phylogeny_install",
                   "--repo", str(repo), "--artifact", str(repo / "benchmarks/work/publication_frozen_overlay_20260926"),
                   "--mafft", str(prefix / "bin/mafft"), "--fasttree",
                   str(installed.parent / "FastTree_v220/FastTree"), "--output", str(phylogeny)]
        env["MAFFT_BINARIES"] = str(prefix / "libexec/mafft")
        report["phylogeny_environment"] = dict(report["environment_overrides"], MAFFT_BINARIES=env["MAFFT_BINARIES"])
        report["commands"].append(dict(argv=command, cwd=str(repo), timeout_seconds=240))
        execute(command, repo, env, output / "phylogeny_driver.log", timeout=240)
        report["phylogeny_report"] = identity(phylogeny / "report.json")
        result = json.loads((phylogeny / "report.json").read_text())
        report["fixture_counts"] = dict(genes=result["result"]["genes"], pairs=result["result"]["pairs"])
        # Check retained inputs again, not generated build intermediates.
        verify(archive, SHA)
        for comparison in comparisons:
            verify(Path(comparison["retained"]["path"]), comparison["retained"]["sha256"])
        for record in report["retained_helpers"]:
            verify(Path(record["path"]), record["sha256"])
        report["status"] = "core_built_and_installed_phylogeny_fixture_passed"
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
    for name in ("repo", "archive", "installed", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    run(args.repo.resolve(), args.archive.resolve(), args.installed.resolve(), args.output.absolute())
