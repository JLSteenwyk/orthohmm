"""Read-only runtime identity snapshot for the corrected FastOMA launch."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.prepare_qfo_corrected_fastoma_assets import JAVA, IMAGE_DIGEST, image_identity, validate_nextflow
from benchmark_tools.snapshot_runtime_trees import inventory, verify

DOCKER = "/usr/bin/docker"
BINARIES = (DOCKER, "/usr/bin/dockerd", "/usr/bin/containerd", "/usr/bin/runc",
            "/usr/libexec/docker/docker-init", "/usr/bin/bash", "/home/bizon/bin/nextflow")
INFO_FIELDS = ("ServerVersion", "Driver", "CgroupDriver", "CgroupVersion", "DefaultRuntime",
               "KernelVersion", "OperatingSystem", "OSType", "Architecture", "SecurityOptions")


def validate_java(text):
    lines = text.splitlines()
    if (not lines or lines[0] != 'openjdk version "17.0.18-internal" 2026-01-20'
            or "OpenJDK Runtime Environment (build 17.0.18-internal+0-adhoc.conda.src)" not in lines):
        raise ValueError("Unexpected Java runtime version/build")


def run(argv, env=None):
    result = subprocess.run(argv, text=True, capture_output=True, check=True, timeout=60, env=env)
    return {"argv": argv, "stdout": result.stdout, "stderr": result.stderr, "exit_code": result.returncode}


def docker_state():
    version = json.loads(run([DOCKER, "version", "--format", "{{json .}}"])["stdout"])
    raw = json.loads(run([DOCKER, "info", "--format", "{{json .}}"])["stdout"])
    info = {key: raw[key] for key in INFO_FIELDS}
    if (info["CgroupVersion"] != "2" or info["CgroupDriver"] != "systemd"
            or info["DefaultRuntime"] != "runc" or info["OSType"] != "linux"
            or version["Client"]["Context"] != "default"):
        raise ValueError("Unexpected Docker context, cgroup or runtime configuration")
    return {"version": version, "info": info}


def inspect():
    # A remote daemon would make local binary fingerprints misleading.
    if any(os.environ.get(key) for key in ("DOCKER_HOST", "DOCKER_CONTEXT")):
        raise ValueError("Require default local Docker connection")
    trees = inventory([JAVA, Path.home() / ".nextflow/capsule/apps/nextflow-all_22.10.8"])
    external_dirs = [r["path"] for r in trees["records"] if r["kind"] == "symlink"
                     and r["path"] in trees["external_symlinks"] and Path(r["resolved"]).is_dir()]
    if external_dirs:
        raise ValueError("Runtime has untraversed external directories: " + repr(external_dirs))
    binaries = [record(path) for path in BINARIES]
    state = docker_state()
    image = image_identity(json.loads(run([DOCKER, "image", "inspect", IMAGE_DIGEST])["stdout"]))
    env = {**os.environ, "JAVA_HOME": str(JAVA), "JAVA_CMD": str(JAVA / "bin/java"),
           "NXF_VER": "22.10.8", "NXF_OFFLINE": "true", "NXF_OPTS": "-Xmx1g"}
    java = run([str(JAVA / "bin/java"), "-version"], env)
    validate_java(java["stderr"])
    nextflow = run(["/home/bizon/bin/nextflow", "-version"], env)
    validate_nextflow(nextflow["stdout"])
    verify(trees)
    for item in binaries:
        check(item)
    if docker_state() != state:
        raise ValueError("Docker identity changed during inspection")
    if image_identity(json.loads(run([DOCKER, "image", "inspect", IMAGE_DIGEST])["stdout"])) != image:
        raise ValueError("Container image identity changed during inspection")
    return {"status": "fastoma_local_runtime_identity_observed", "source": record(__file__),
            "runtime_trees": trees, "binaries": binaries, "docker": state, "image": image,
            "java_probe": java, "nextflow_probe": nextflow,
            "execution_authorized": False, "accuracy_evaluated": False, "publication_ready": False,
            "limitations": ["Identity of selected installed runtime files, not a hermetic host snapshot.",
                "Host dynamic libraries, kernel and active daemon executable mappings are not fully inventoried.",
                "External file symlinks include target hashes; directory symlinks outside the roots are rejected.",
                "No container or inference job started, no Docker configuration changed.",
                "Resource enforcement remains supported by the separate actual container probe.",
                "Recheck runtime identity before and after full inference; shared-host timing is not matched timing."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = inspect()
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
