"""Recover retained historical commands without implying execution attestation."""

import argparse
import json
import math
from pathlib import Path

from benchmark_tools.assemble_benchmark_provenance import Reader
from benchmark_tools.audit_ob_orthofinder_provenance import log_command
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.summarize_matched_resources import verbose_time


REGISTER_SHA = "48be4b216f9f3776354246894c8c2a38f790fd4f075528e88f56b230a0a3c47c"
PARITY_SHA = "beba94cba53eef41c76c9d525eed48fb52287b98fc7e632c11c2894319e8a75a"
TEXT_FILES = (
    "source_commit.txt", "tool_version.txt", "runtime_note.txt",
    "inference_wall_time_s.txt", "wall_time_s.txt", "runner_wall_time_s.txt",
    "blast_wall_time_s.txt", "bpo_conversion_wall_time_s.txt",
    "bpo_index_wall_time_s.txt", "downstream_wall_time_s.txt",
)
TIME_SCOPES = {
    "time.log": "native command; FastOMA CPU/RSS are driver/waited descendants, not aggregate Docker tasks",
    "conversion_time.log": "checkpoint conversion only; not sequence-only inference",
    "bpo_converter.time.log": "BLAST-to-BPO conversion only",
    "downstream_time.log": "recovered mode4 downstream only; excludes BLAST/BPO",
}


def finite_nonnegative(value):
    if type(value) not in (int, float) or not math.isfinite(value) or value < 0:
        raise ValueError("Invalid resource number")
    return value


def metadata_events(text):
    # Recovery keys repeat intentionally; a dict would erase earlier attempts.
    events = []
    for line in text.splitlines():
        fields = line.split("\t")
        if len(fields) != 2 or not all(fields):
            raise ValueError("Invalid historical metadata line")
        events.append(dict(key=fields[0], value=fields[1]))
    return events


def metrics_details(metrics, selected):
    harness = metrics["harness"]
    if harness["exit_code"] != 0 or metrics["status"] != "complete":
        raise ValueError("Historical metrics did not complete")
    bound = [p for p in harness["output_manifest"]
             if p["sha256"] == selected["sha256"] and p["bytes"] == selected["bytes"]]
    if not bound:
        raise ValueError("Metrics output does not bind selected prediction")
    wall = finite_nonnegative(metrics["wall_s"])
    start = finite_nonnegative(metrics["started_at_epoch_s"])
    finish = finite_nonnegative(metrics["finished_at_epoch_s"])
    if not math.isclose(finish - start, wall, rel_tol=0, abs_tol=1e-6):
        raise ValueError("Metrics wall interval differs")
    for command in (metrics["command"], harness["command"]):
        if not isinstance(command, list) or not command or not all(isinstance(x, str) and x for x in command):
            raise ValueError("Invalid recorded metrics argv")
    return dict(
        entrypoint_argv=metrics["command"], harness_argv=harness["command"],
        measurement={"elapsed_seconds": wall,
            "user_seconds": finite_nonnegative(metrics["user_cpu_s"]),
            "system_seconds": finite_nonnegative(metrics["system_cpu_s"])},
        memory={"value": finite_nonnegative(metrics["peak_process_tree_rss_bytes"]),
                "unit": "bytes", "scope": metrics["rss_measurement"]},
        recorded_cwd=metrics["cwd"], harness_revision_at_completion=harness["git_commit"],
        harness_git_dirty=harness["git_dirty"], matched_output_manifest=bound,
        inherited_input_manifest=harness["input_manifest"],
        inherited_source_manifest=harness["source_manifest"],
        cpu_scope="Historical metrics collector's recorded CPU totals; no fresh process-tree/cgroup accounting validation",
        limitations=["Manifests were recorded by the historical harness; raw inputs and source bytes are not re-admitted here.",
                     "Launch source revision and harness-at-completion revision are distinct historical fields."])


def audit(repo, output):
    repo, output = Path(repo).resolve(), Path(output)
    if output.exists():
        raise FileExistsError(output)
    reader = Reader()
    base = repo / "benchmark_tools/results"
    register_pin = record(base / "all_benchmark_provenance_20261004_v2/register.json")
    parity_pin = record(base / "three_kingdoms_parity_20260907.json")
    if register_pin["sha256"] != REGISTER_SHA or parity_pin["sha256"] != PARITY_SHA:
        raise ValueError("Changed selected historical source")
    register, parity = reader.read(register_pin), reader.read(parity_pin)
    selected = {r["key"]: r for r in register["rows"] if r["dataset"] == "ThreeKingdoms"}
    if len(selected) != 8 or len(parity["methods"]) != 8:
        raise ValueError("Incomplete historical inventory")
    rows = []
    for historical in parity["methods"]:
        key = historical["key"]
        if key == "sonicparanoid_2_0_9":
            continue
        row = selected[key]
        prediction = row["output_records"][0]
        if historical["provenance"]["orthogroups_sha256"] != prediction["sha256"]:
            raise ValueError("Historical directory belongs to a different prediction")
        directory = repo / historical["result_dir"]
        observed = record(directory / "orthogroups.txt")
        if observed != prediction:
            raise ValueError("Retained selected prediction changed")
        reader.checked[observed["path"]] = observed
        values, evidence, timing = {}, {}, []
        for name in TEXT_FILES:
            path = directory / name
            if path.exists():
                pin = record(path)
                values[name], evidence[name] = reader.text(pin).strip(), pin
        events = []
        if (directory / "run_metadata.tsv").exists():
            pin = record(directory / "run_metadata.tsv")
            events, evidence["run_metadata.tsv"] = metadata_events(reader.text(pin)), pin
        metrics = None
        if (directory / "metrics.json").exists():
            pin = record(directory / "metrics.json")
            metrics, evidence["metrics.json"] = metrics_details(reader.read(pin), prediction), pin
        empty_logs = []
        for name, scope in TIME_SCOPES.items():
            path = directory / name
            if not path.exists():
                continue
            pin = record(path)
            text = reader.text(pin)
            evidence[name] = pin
            if not text.strip():
                empty_logs.append(name)
                continue
            timing.append(dict(scope=scope, argv=log_command(text, "Command being timed: "),
                               measurement=verbose_time(text), evidence=pin,
                               memory_scope="GNU-time maximum process RSS; not aggregate concurrent workflow memory"))
        rows.append(dict(key=key, selected_prediction=prediction,
            declared_version=row["declared_version"], retained_text=values,
            evidence=evidence, metadata_events=events, metrics=metrics,
            timing_records=timing, empty_timing_logs=empty_logs,
            previously_reported_resources=row["resources"],
            binding_scope="Retained result-directory association plus selected prediction identity; not immutable historical execution attestation"))
    if len(rows) != 7 or len({r["key"] for r in rows}) != 7:
        raise ValueError("Require seven non-Sonic historical rows")
    # Recheck bytes without parsing large prediction text or duplicating manifests.
    for pin in reader.checked.values():
        check(pin)
    result = dict(status="historical_three_kingdoms_metadata_supplement", rows=rows,
        checked_records=list(reader.checked.values()), source=record(__file__),
        helpers=[record(Path(__file__).with_name(n)) for n in
                 ("assemble_benchmark_provenance.py", "audit_ob_orthofinder_provenance.py",
                  "prepare_ob_candidate_neighborhood.py", "summarize_matched_resources.py", "gnu_time_companion.py")],
        scores_recomputed=False, native_runs_repeated=False, publication_ready=False,
        historical_consumption_proven=False, controlled_comparative_resources=False,
        limitations=["New direct metadata readback; no inference, conversion or scoring is re-executed.",
            "Contemporary matched Sonic evidence is already bound in the selected register; older Sonic logs are deliberately excluded.",
            "Fractional GNU-time durations and integer wrapper/stage intervals remain separate observations, not replacements or independent repeats.",
            "OrthoFinder MCL checkpoint timestamp, checkpoint conversion and full pipeline retain different scopes.",
            "Repeated recovery metadata keys preserve all ordered events, including later-invalidated attempts; exit0 alone does not validate their science.",
            "Executable names, version text and source revisions do not prove complete historical executable or dependency identity."])
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    result = audit(args.repo, args.output)
    print(json.dumps({"status": result["status"], "rows": len(result["rows"]),
                      "direct_records": len(result["checked_records"])}))
