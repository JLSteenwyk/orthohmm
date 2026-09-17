"""Snapshot the existing local QfO scorer and references without downloading or scoring."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_qfo_recovered_pairs import record

JAVA = Path("/home/bizon/anaconda3/pkgs/openjdk-17.0.18-ha668962_0/lib/jvm")
NEXTFLOW = Path("/home/bizon/bin/nextflow")
SINGULARITY = Path("/usr/local/bin/singularity")
IMAGES = ("qfobenchmark-python-2022.1.img", "qfobenchmark-darwin-2022.1.img", "qfobenchmark-fas_benchmark-2022.1.img")


def inventory(directory):
    files = sorted(p for p in directory.rglob("*") if p.is_file())
    if not files:
        raise ValueError("Empty required runtime/reference directory")
    return [record(p) for p in files]


def environment():
    return {**os.environ, "JAVA_HOME": str(JAVA), "JAVA_CMD": str(JAVA / "bin/java"), "NXF_OFFLINE": "true"}


def freeze(root, output):
    if output.exists():
        raise FileExistsError(output)
    pipeline = root / "qfo_benchmark/benchmark-webservice"
    config = root / "benchmark_tools/results/qfo_recovered_assessment.config"
    cache = root / "qfo_benchmark/scoring/container_cache"
    tracked = subprocess.check_output(["git", "-C", str(pipeline), "ls-files", "-z"]).decode().split("\0")
    sources = [record(pipeline / name) for name in tracked if name]
    subprocess.run(["git", "-C", str(pipeline), "diff", "--exit-code", "HEAD", "--"], check=True)
    env = environment()
    versions = {}
    for label, command in (("nextflow", [str(NEXTFLOW), "-version"]),
                           ("java", [str(JAVA / "bin/java"), "-version"]),
                           ("singularity", [str(SINGULARITY), "--version"])):
        result = subprocess.run(command, env=env, text=True, capture_output=True, check=True)
        versions[label] = {"command": command, "stdout": result.stdout, "stderr": result.stderr}
    if "version 22.10.8 build 5860" not in versions["nextflow"]["stdout"]:
        raise ValueError("Installed Nextflow differs from historical scorer version")
    config_command = [str(NEXTFLOW), "-c", str(config), "config", str(pipeline), "-profile", "singularity", "-flat"]
    effective = subprocess.run(config_command, env=env, text=True, capture_output=True, check=True)
    for name in IMAGES:
        if str(cache / name) not in effective.stdout:
            raise ValueError("Effective scorer config lacks pinned local image")
    report = {"status": "local_qfo_assessment_environment_frozen", "accuracy_evaluated": False,
              "source": record(__file__), "pipeline": str(pipeline),
              "pipeline_commit": subprocess.check_output(["git", "-C", str(pipeline), "rev-parse", "HEAD"], text=True).strip(),
              "pipeline_files": sources, "reference_files": inventory(pipeline / "reference_data/2020"),
              "java_files": inventory(JAVA), "images": [record(cache / n) for n in IMAGES],
              "executables": [record(NEXTFLOW), record(SINGULARITY)],
              "singularity_support": [record(Path("/usr/local/libexec/singularity/bin") / name)
                                      for name in ("starter", "starter-suid", "squashfuse_ll")],
              "singularity_config": record("/usr/local/etc/singularity/singularity.conf"),
              "execution_config": record(config), "effective_config": effective.stdout, "config_command": config_command,
              "versions": versions, "environment_overrides": {k: env[k] for k in ("JAVA_HOME", "JAVA_CMD", "NXF_OFFLINE")},
              "limitations": ["Local reference/runtime identity, not a hermetic operating-system snapshot.",
                              "Host libraries and kernel are not fully captured; container images and Java runtime are checksum-pinned.",
                              "Historical source/image identity is not retroactively established by a present-day snapshot."]}
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    freeze(args.root.resolve(), args.output.resolve())
