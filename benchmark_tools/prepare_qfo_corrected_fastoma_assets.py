"""Pin FastOMA assets without substituting the missing corrected species tree."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.prepare_qfo_corrected_proteinortho import PRIMARY_SHA
from benchmark_tools.run_qfo_corrected_sonic import verify_corrected_inputs
from benchmark_tools.snapshot_runtime_trees import inventory

WORKFLOW_COMMIT = "bf6dcbaa8cf516ab6f6e074dba37eceb59a9b80e"
IMAGE_ID = "sha256:dd9c319b85e65a38f45f033d85663211eedccd0bfcc9c25df8048924bb4862dd"
IMAGE_DIGEST = "dessimozlab/fastoma@sha256:67c802a71c2d150ba8825a8c5dbc5a8a68e1dd42554e6695038043d666ac358b"
JAVA = Path("/home/bizon/anaconda3/pkgs/openjdk-17.0.18-ha668962_0/lib/jvm")


def image_identity(data):
    if len(data) != 1:
        raise ValueError("Require one FastOMA image")
    row = data[0]
    if (row["Id"] != IMAGE_ID or IMAGE_DIGEST not in row["RepoDigests"]
            or row["Architecture"] != "amd64" or row["Os"] != "linux"
            or row["Config"]["Labels"]["org.opencontainers.image.version"] != "0.3.5"):
        raise ValueError("Changed FastOMA image identity")
    return {key: row[key] for key in ("Id", "RepoDigests", "Architecture", "Os", "Config", "RootFS")}


def validate_nextflow(text):
    if "version 22.10.8 build 5860" not in text:
        raise ValueError("Unexpected Nextflow runtime")


def prepare(root, destination):
    if destination.exists():
        raise FileExistsError(destination)
    results = root / "benchmark_tools/results"
    primary_path = results / "qfo_corrected_primary_commands_20260918.json"
    primary = read_frozen(primary_path, PRIMARY_SHA)
    inputs = json.loads((Path(primary["input_directory"]) / "staging_manifest.json").read_text())["input_fastas"]
    verify_corrected_inputs(primary, inputs)
    workflow = root.parents[1] / "SOFTWARE/FastOMA"
    commit = subprocess.check_output(["git", "-C", str(workflow), "rev-parse", "HEAD"], text=True).strip()
    if commit != WORKFLOW_COMMIT:
        raise ValueError("Changed FastOMA workflow revision")
    subprocess.run(["git", "-C", str(workflow), "diff", "--exit-code", "HEAD", "--", "FastOMA.nf", "nextflow.config", "conf"], check=True)
    records = [record(workflow / name) for name in ("FastOMA.nf", "nextflow.config", "conf/base.config")]
    records += [record(root / name) for name in (
        "benchmark_tools/fastoma_qfo.config", "benchmark_tools/run_fastoma_qfo.slurm",
        "benchmark_tools/recover_fastoma_qfo_collection.slurm",
        "qfo_benchmark/results/fastoma/work/06/276078d39a60b51341cb350ebb63ad/.command.sh")]
    records += [record("/home/bizon/bin/nextflow"), record(JAVA / "bin/java"), record(primary_path)]
    db = record(root / "qfo_benchmark/data/LUCA.h5")
    image = image_identity(json.loads(subprocess.check_output(["docker", "image", "inspect", "dessimozlab/fastoma:0.3.5"], text=True)))
    package = subprocess.check_output(["docker", "run", "--rm", "--network", "none", "--read-only", IMAGE_ID,
        "python3", "-c", "import importlib.metadata; print(importlib.metadata.version('FastOMA'))"], text=True, timeout=60).strip()
    if package != "0.3.5":
        raise ValueError("Container package version differs")
    env = os.environ.copy()
    env.update(NXF_OFFLINE="true", NXF_VER="22.10.8", JAVA_HOME=str(JAVA), JAVA_CMD=str(JAVA / "bin/java"))
    probe = subprocess.run(["/home/bizon/bin/nextflow", "-version"], env=env, cwd=root,
                           text=True, capture_output=True, check=True, timeout=60)
    validate_nextflow(probe.stdout)
    capsule = inventory([Path.home() / ".nextflow/capsule/apps/nextflow-all_22.10.8"])
    for item in [*records, db]:
        check(item)
    report = {"status": "corrected_fastoma_assets_pinned_waiting_for_corrected_tree",
              "execution_authorized": False, "accuracy_admitted": False, "source": record(__file__),
              "workflow_commit": commit, "workflow_root": str(workflow), "sources": records,
              "omamer_database": db, "image": image, "container_package_version": package,
              "nextflow_probe": {"stdout": probe.stdout, "stderr": probe.stderr, "exit_code": probe.returncode},
              "nextflow_capsule_inventory": capsule, "input_fastas": inputs,
              "species_tree": None, "diagnostic": "Supplied species tree from corrected full OrthoFinder; not independent tree inference.",
              "required_tree_origin": primary["methods"]["orthofinder_full"]["output"],
              "remaining_gates": [
                  "Admit corrected full OrthoFinder output and bind its rooted species tree with exact 78-species scope.",
                  "Freeze fresh proteome/tree copies, exact Nextflow commands, immutable container digest selection and effective configuration.",
                  "Pin Java/Docker/Nextflow runtime and resource enforcement; isolate cache, logs and working directories.",
                  "Validate native completion, IDs and pair-generation semantics; never silently discard malformed pair rows.",
                  "Native output admission, conversion and corrected QfO scoring remain separate.",
              ],
              "limitations": ["No input tree substituted from the original release.",
                              "Current image identity does not alone prove historical container identity.",
                              "Asset identity is not a complete runtime or executable command freeze.",
                              "The retained workflow omits min-sequence-length in infer-roothogs, matching the historical native task."]}
    with destination.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.root.resolve(), args.output.resolve())
