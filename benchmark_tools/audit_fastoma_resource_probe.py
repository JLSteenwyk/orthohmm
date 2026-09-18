"""Verify a tiny Nextflow/Docker resource probe, not FastOMA inference."""

import argparse
import csv
import io
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.prepare_qfo_corrected_fastoma_assets import IMAGE_DIGEST, JAVA, image_identity
from benchmark_tools.run_simulation_methods import read_frozen

ASSETS_SHA = "bdb6878dbab9396c6cb5360d35625ebe44149cf422b1d18dbf9b0e693a2b7975"


def validate_probe(trace, output, wrapper, exit_code):
    rows = list(csv.DictReader(io.StringIO(trace), delimiter="\t"))
    if (len(rows) != 1 or rows[0]["name"] != "resource_probe"
            or rows[0]["status"] != "COMPLETED" or rows[0]["exit"] != "0" or exit_code.strip() != "0"):
        raise ValueError("Probe did not complete exactly once")
    limits = json.loads(output)
    fields = limits["cpu_max"].split()
    if len(fields) != 2 or any(not v.isdigit() for v in fields):
        raise ValueError("Missing finite cgroup CPU quota")
    quota, period = map(int, fields)
    if not quota == period > 0 or limits["memory_max"] != "268435456":
        raise ValueError("Unexpected observed CPU/memory limit")
    commands = [shlex.split(s.strip()) for s in wrapper.splitlines() if s.strip().startswith("docker run ")]
    if len(commands) != 1 or commands[0].count(IMAGE_DIGEST) != 1:
        raise ValueError("Missing unique immutable Docker invocation")
    argv = commands[0]
    expected = {"--cpus": "1.0", "--memory": "256m", "--network": "none"}
    for flag, value in expected.items():
        if argv.count(flag) != 1 or argv[argv.index(flag) + 1] != value:
            raise ValueError("Unexpected Docker resource/network option: " + flag)
    if "--privileged" in argv:
        raise ValueError("Unexpected privileged container")
    return {"trace": rows[0], "observed_limits": limits, "docker_argv": argv}


def validate_config(text):
    # Nextflow 22.10.8 exposes flat/properties output, not JSON configuration export.
    values = {}
    for line in text.splitlines():
        key, separator, value = line.partition(" = ")
        if separator:
            if key in values:
                raise ValueError("Duplicate effective configuration key")
            values[key] = value
    expected = {"manifest.version": "'0.3.5'", "process.container": repr(IMAGE_DIGEST),
        "process.executor": "'local'", "executor.cpus": "180", "executor.memory": "'700 GB'",
        "params.max_cpus": "180", "params.max_memory": "'700.GB'",
        "params.fasta_header_id_transformer": "'UniProt'", "params.force_pairwise_ortholog_generation": "true",
        "params.filter_method": "'col-row-threshold'", "params.filter_gap_ratio_row": "0.3",
        "params.filter_gap_ratio_col": "0.5", "params.nr_repr_per_hog": "5",
        "params.min_sequence_length": "40", "docker.enabled": "true", "docker.runOptions": "'--network none'",
        "process.'withName:omamer_run'.memory": "'24 GB'", "process.'withName:omamer_run'.maxForks": "28",
        "process.'withName:collect_subhogs'.cpus": "16", "process.'withName:collect_subhogs'.memory": "'280 GB'",
        "process.'withName:extract_pairwise_ortholog_relations'.cpus": "16",
        "process.'withName:extract_pairwise_ortholog_relations'.memory": "'280 GB'", "trace.enabled": "true"}
    for key, value in expected.items():
        if values.get(key) != value:
            raise ValueError("Effective FastOMA configuration differs: " + key)
    return expected


def audit(root, probe):
    assets_path = root / "benchmark_tools/results/qfo_corrected_fastoma_assets_20260918.json"
    assets = read_frozen(assets_path, ASSETS_SHA)
    config = root / "benchmark_tools/fastoma_corrected_execution.config"
    script = root / "benchmark_tools/fastoma_resource_probe.nf"
    checked = [record(assets_path), record(config), record(script), assets["omamer_database"]]
    checked += [r for r in assets["sources"] if Path(r["path"]).name in
                ("nextflow.config", "base.config", "fastoma_qfo.config", "nextflow", "java")]
    wrappers = list((probe / "work").glob("*/*/.command.run"))
    if len(wrappers) != 1:
        raise ValueError("Require one probe task directory")
    task = wrappers[0].parent
    paths = [probe / "trace.txt", probe / ".nextflow.log",
             *[task / name for name in (".command.run", ".command.sh", ".command.out", ".exitcode")]]
    checked.extend(record(p) for p in paths)
    for item in checked:
        check(item)
    result = validate_probe((probe / "trace.txt").read_text(), (task / ".command.out").read_text(),
                            (task / ".command.run").read_text(), (task / ".exitcode").read_text())
    if "Version: 22.10.8 build 5860" not in (probe / ".nextflow.log").read_text():
        raise ValueError("Probe Nextflow version differs")
    image = image_identity(json.loads(subprocess.check_output(["docker", "image", "inspect", IMAGE_DIGEST], text=True)))
    command = ["/home/bizon/bin/nextflow", "-C", str(config), "config", assets["workflow_root"], "-flat"]
    env = {**os.environ, "JAVA_HOME": str(JAVA), "JAVA_CMD": str(JAVA / "bin/java"),
           "NXF_OFFLINE": "true", "NXF_VER": "22.10.8", "NXF_OPTS": "-Xmx1g"}
    completed = subprocess.run(command, cwd=probe.parent, env=env, capture_output=True, text=True, check=True, timeout=60)
    selected = validate_config(completed.stdout)
    db_path = root / "qfo_benchmark/data/LUCA.h5"
    db_setting = "params.omamer_db = " + repr(str(db_path))
    if (db_setting not in completed.stdout.splitlines()
            or str(db_path.resolve()) != assets["omamer_database"]["path"]):
        raise ValueError("Configuration points to a different OMAmer database")
    for item in checked:
        check(item)
    return {"status": "fastoma_corrected_resource_probe_verified", "source": record(__file__),
        "helper": record(Path(__file__).with_name("prepare_qfo_corrected_fastoma_assets.py")),
        "checked_records": checked, "probe": result, "image": image,
        "configuration_command": command, "effective_configuration": completed.stdout,
        "configuration_stderr": completed.stderr, "validated_settings": selected,
        "inference_authorized": False, "accuracy_evaluated": False, "publication_ready": False,
        "limitations": ["A one-task cgroup-v2 probe on bizon, not biological inference or aggregate resource validation.",
            "Flat configuration validates selected settings, not evaluation of all dynamic task-resource closures.",
            "Production must reserve 20 GiB beyond the 700 GiB task pool for the controller and overhead.",
            "Corrected admitted tree, full runtime freeze, fresh input copies and native-output admission remain required."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("root", "probe", "output"):
        parser.add_argument("--" + flag, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.root.resolve(), args.probe.resolve())
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
