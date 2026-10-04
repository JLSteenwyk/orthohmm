"""Consolidate all retained benchmark rows without re-admitting raw inference."""

import argparse
import csv
import json
import math
from pathlib import Path

from benchmark_tools.audit_ob_orthofinder_provenance import log_command
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.summarize_matched_resources import verbose_time


SOURCES = {
    "scores": ("current_benchmark_scores_20260926_v2/manifest.json", "8d56dae74721c4530ab2b86c5008b2d0ebac07878e574fd85b1ad84e4b908f07"),
    "ob": ("orthobench_provenance_register_20260927/register.json", "21b1a736521ddfd87e3def36387258b8f509a73ea7ab58af8678bb101aca490d"),
    "qfo": ("qfo_corrected_comparison_20260926_v7/manifest.json", "042aa221d01554114a1b0b413ca3e2ff56dad5f8c06362c69a1bf49f2de642fc"),
    "three": ("three_kingdoms_comparison_matched_20260918/comparison.json", "ba1ed664ad8ee89ba65b72e7d4d4341eed864957b0720ab1ff4fbe91a398945a"),
    "parity": ("three_kingdoms_parity_20260907.json", "beba94cba53eef41c76c9d525eed48fb52287b98fc7e632c11c2894319e8a75a"),
    "sonic": ("three_kingdoms_sonic_matched_assessment_21796.json", "aeddd422a5c733ac40bef0c4953a6bda12ab2a751baeeae9e9c28cac472e9653"),
    "three_inputs": ("three_kingdoms_method_inputs_20260918.json", "5dfd585193a6a576f2685388827bd6f8e3f5bad2ffcc715aaf665247cd2345dc"),
    "mcl_inputs": ("three_kingdoms_orthomcl_inputs_20260927.json", "3f169fd4f8036ba68238dfb7e386b7decbbafbe9f47356a1a53ac2926c2c3f84"),
    "of": ("qfo_orthofinder_provenance_consolidated_20260926.json", "2333aabd32e39fe9c845b05398b1fd8b0efd466a3b54d7dd8adc1f34940d27ea"),
}
KEYS = (
    "orthohmm_high_sensitivity", "orthohmm_phylogeny_satellite_v2",
    "orthofinder_3_1_5_full", "orthofinder_3_1_5_sequence_only",
    "sonicparanoid_2_0_9", "proteinortho_6_3_6", "fastoma_0_3_5", "orthomcl_1_4",
)


def inventory(rows, field="key"):
    result = {r[field]: r for r in rows}
    if len(rows) != 8 or set(result) != set(KEYS):
        raise ValueError("Require exactly eight unique retained methods")
    return result


def same_score(actual, expected):
    if (type(actual) not in (int, float) or type(expected) not in (int, float)
            or not math.isfinite(actual) or not math.isfinite(expected)
            or not 0 <= actual <= 1 or not math.isclose(actual, expected, rel_tol=0, abs_tol=1e-12)):
        raise ValueError("Score row differs from selected current evidence")


class Reader:
    def __init__(self):
        self.checked = {}

    def text(self, pin):
        check(pin)
        text = Path(pin["path"]).read_text()
        check(pin)
        old = self.checked.get(pin["path"])
        if old is not None and old != pin:
            raise ValueError("Conflicting evidence identity")
        self.checked[pin["path"]] = pin
        return text

    def read(self, pin):
        return json.loads(self.text(pin))


def execution_pin(native):
    if "execution" in native:
        return native["execution"]
    candidates = [r for r in native.get("checked_records", []) if Path(r["path"]).name == "execution.json"]
    unique = {r["path"]: r for r in candidates}
    if len(unique) != 1:
        raise ValueError("Missing or ambiguous native execution binding")
    return next(iter(unique.values()))


def timed_execution(reader, pin, scope, scheduler=None):
    execution = reader.read(pin)
    recovered = "native_timing" in execution
    command = execution["command"] if recovered else execution["native_argv"]
    exit_code = execution["native_exit_code"] if recovered else execution["exit_code"]
    if exit_code != 0:
        raise ValueError("Native execution failed")
    if scheduler is not None and (
        scheduler["JobIDRaw"] != execution["job_id"]
        or scheduler["State"] != "COMPLETED" or scheduler["ExitCode"] != "0:0"
    ):
        raise ValueError("Scheduler/native execution mismatch")
    timing = execution["native_timing"] if recovered else execution["timing"]
    text = reader.text(timing)
    if log_command(text, "Command being timed: ") != command:
        raise ValueError("Timing command differs from bound native execution")
    return dict(scope=scope, evidence=pin, timing_evidence=timing, command=command,
        job_id=execution["job_id"], measurement=verbose_time(text),
        memory_scope="GNU-time maximum process RSS; not aggregate concurrent process-tree/cgroup memory",
        input_records=execution.get("copied_inputs", execution.get("staged_inputs", [])),
        runtime_records=[execution[k] for k in ("runtime", "runtime_before", "runtime_after") if k in execution],
        full_inference=scope == "full native inference", scheduler=scheduler)


def bind_conversion(row, conversion):
    if (conversion["semantics"] != row["prediction_semantics"]
            or conversion["total_pairs"] != row["submitted_pairs"]
            or conversion["retained_pairs"] != row["retained_pairs"]
            or conversion["removed_mapping_pairs"] != row["removed_mapping_pairs"]
            or conversion["participant"] != row.get("participant", conversion["participant"])):
        raise ValueError("QfO conversion binding differs")
    return conversion["filtered_pairs"]


def qfo_details(row, reader, orthofinder):
    conversion = reader.read(row["conversion"])
    pair_file = bind_conversion(row, conversion)
    resources, commands = [], []
    inputs = conversion.get("input_fastas", [])
    native_pin = conversion.get("admission", conversion.get("native_admission", conversion.get("candidate_admission")))
    gaps = ["Complete transitive executable/input-consumption provenance is not certified here."]
    if row["key"].startswith("orthofinder"):
        bound = next(r for r in orthofinder["rows"] if r["key"] == row["key"])
        if bound["pair_file"] != pair_file or bound["native_admission"] != native_pin or bound["scores"] != row["scores"]:
            raise ValueError("OrthoFinder consolidated run/score binding differs")
        commands.append(dict(scope="shared full-run native command", argv=orthofinder["native_command"]))
        if row["key"].endswith("_full"):
            resources.append(dict(scope="full native inference", full_inference=True,
                measurement=orthofinder["full_native_resources"], scheduler=orthofinder["native_scheduler"],
                memory_scope="GNU-time maximum process RSS, not aggregate memory"))
        else:
            gaps.append("No separately timed sequence-only inference; same full-run checkpoint.")
        resources.append(dict(scope="pair conversion", full_inference=False,
            measurement={"elapsed_seconds": bound["conversion_wall_seconds"]}))
        inputs = [r for r in orthofinder["checked_records"] if r["path"].endswith(".fasta")]
    elif row["key"] in {"sonicparanoid_2_0_9", "proteinortho_6_3_6", "fastoma_0_3_5", "orthomcl_1_4"}:
        native = reader.read(native_pin)
        recovered = row["key"] == "orthomcl_1_4"
        scope = "recovered downstream mode4 only; excludes BLAST/BPO" if recovered else "full native inference"
        resource = timed_execution(reader, execution_pin(native), scope, native.get("scheduler"))
        resources.append(resource)
        commands.append(dict(scope=scope, argv=resource["command"]))
        inputs = native.get("input_fastas", [r for r in native.get("checked_records", []) if r["path"].endswith(".fasta")])
        if not inputs:
            inputs = resource["input_records"]
        if recovered:
            gaps.append("No uninterrupted full search-to-clusters runtime; recovered downstream time is not full inference.")
        if row["key"] == "fastoma_0_3_5":
            resource["memory_scope"] = "GNU-time driver maximum RSS; Docker tasks not aggregate workflow memory"
            resource["cpu_scope"] = "GNU-time driver/waited descendants; not aggregate Docker task CPU"
            gaps.append("Supplied OrthoFinder species tree; GNU-time CPU/RSS describe the driver and waited descendants, not all Docker workflow tasks.")
    else:
        gaps.append("Cached corrected factorial output; full search-to-this-output inference time is unestablished.")
        if row["key"] == "orthohmm_phylogeny_satellite_v2":
            native = reader.read(native_pin)
            if native["native_pairs"] != conversion["native_input"]:
                raise ValueError("Phylogenetic native-pair binding differs")
            gaps.append("Reconciliation admission is incremental; do not substitute its scheduler elapsed time for full inference.")
    if "started_epoch" in conversion:
        seconds = conversion["finished_epoch"] - conversion["started_epoch"]
        if not math.isfinite(seconds) or seconds < 0:
            raise ValueError("Invalid conversion interval")
        resources.append(dict(scope="pair conversion", full_inference=False, measurement={"elapsed_seconds": seconds}))
    return dict(output_records=[pair_file], native_admission=native_pin, input_records=inputs,
        commands=commands, resources=resources, gaps=gaps,
        input_note="Corrected-release inputs; retained pins/admissions, not fresh raw-input revalidation.")


def three_details(row, historical, inputs, sonic, reader):
    if row["key"] == "sonicparanoid_2_0_9":
        if row["groups"] != sonic["normalized"] or row["counts"] != sonic["counts"]:
            raise ValueError("Matched Sonic score/output binding differs")
        resource = timed_execution(reader, execution_pin(sonic), "full native inference", sonic["scheduler"])
        return dict(resources=[resource], commands=[dict(scope=resource["scope"], argv=resource["command"])],
            input_records=resource["input_records"], input_note="Contemporary matched12-proteome native admission; not historical Sonic timing.",
            gaps=["Maximum process RSS is not aggregate concurrent-process memory."])
    if historical["provenance"]["orthogroups_sha256"] != row["groups"]["sha256"]:
        raise ValueError("Historical timing would bind a different prediction")
    same_score(historical["score"]["f_score"], row["counts"]["f_score"])
    performance = historical["performance"]
    scope = performance["runtime_kind"]
    checkpoint = row["key"].endswith("sequence_only")
    resources = [dict(scope=scope, full_inference=not checkpoint and scope == "measured",
        measurement={"elapsed_seconds": performance["wall_s"], "max_process_rss_kib": performance["peak_rss_kib"]},
        memory_scope=performance["memory_measurement"], requested_cpus=performance["cpus_requested"])]
    return dict(resources=resources, commands=[], input_records=inputs["evidence"] + [r["file"] for r in inputs["copies"]],
        input_note=row["input_status"], gaps=["Historical native argv and complete runtime identity are not consolidated.",
            "Historical input consumption is not proven by current copies/hash manifests."],
        historical_run_metadata=historical["provenance"])


def assemble(repo, output):
    repo, output = Path(repo).resolve(), Path(output)
    if output.exists():
        raise FileExistsError(output)
    reader, documents = Reader(), {}
    for key, (name, sha) in SOURCES.items():
        pin = record(repo / "benchmark_tools/results" / name)
        if pin["sha256"] != sha:
            raise ValueError("Changed retained source: " + key)
        documents[key] = reader.read(pin)
    selected = inventory(documents["scores"]["rows"])
    ob = inventory(documents["ob"]["rows"])
    qfo = inventory(documents["qfo"]["methods"])
    three = inventory([r for r in documents["three"]["rows"] if r["use"] == "comparison"])
    parity = inventory(documents["parity"]["methods"])
    three_inputs = inventory(documents["three_inputs"]["methods"], "method")
    rows = []
    for dataset in ("OrthoBench", "QfO", "ThreeKingdoms"):
        for key in KEYS:
            score = selected[key]["scores"]
            version = parity[key]["version"]
            if dataset == "OrthoBench":
                source = ob[key]
                same_score(source["weighted_refog_f1"], score[dataset])
                if source["prediction"] != selected[key]["orthobench_retained_evidence"]["prediction_provenance"]:
                    raise ValueError("OrthoBench prediction binding differs")
                details = dict(output_records=[source["prediction"]], input_records=[],
                    input_note=source["input_note"], commands=source.get("command"),
                    retained_ob_metadata=source, gaps=source.get("limitations", []) + ["Historical transitive provenance remains partial."],
                    resources=[dict(scope=source["timing_basis"], full_inference=source["inference_wall_seconds"] is not None,
                        measurement={"elapsed_seconds": source["inference_wall_seconds"], "max_process_rss_kib": source["peak_process_rss_kib"]},
                        retained_measurements=source.get("resources"), replay_wall_seconds=source.get("replay_wall_seconds"),
                        conversion_wall_seconds=source["conversion_wall_seconds"])])
                values, semantics = {"weighted_refog_F1": score[dataset]}, source["output_semantics"]
            elif dataset == "QfO":
                source = qfo[key]
                if source["status"] != "admitted":
                    raise ValueError("QfO row not admitted")
                for endpoint, value in source["scores"].items():
                    same_score(value, score[endpoint])
                same_score(source["secondary_mean"], score["QfO_secondary_mean"])
                details = qfo_details(source, reader, documents["of"])
                details.update(score_admission=source["admission"], conversion=source["conversion"],
                    submitted_pairs=source["submitted_pairs"], retained_pairs=source["retained_pairs"],
                    secondary_mean=source["secondary_mean"], score_details=source["details"])
                values, semantics = source["scores"], source["prediction_semantics"]
            else:
                source = three[key]
                same_score(source["counts"]["f_score"], score[dataset])
                if source["groups"] != selected[key]["three_kingdoms"]["groups"]:
                    raise ValueError("Three Kingdoms prediction binding differs")
                details = three_details(source, parity[key], three_inputs[key], documents["sonic"], reader)
                details.update(output_records=[source["groups"]], counts=source["counts"], run=source["run"])
                if key == "orthomcl_1_4":
                    details["supplemental_input_audit"] = record(repo / "benchmark_tools/results" / SOURCES["mcl_inputs"][0])
                    details["input_note"] += "; later native copies/all.fa/species-map content audit available"
                values, semantics = {"BUSCO_reference_pair_F1": score[dataset]}, source["semantics"]
            rows.append(dict(dataset=dataset, key=key, label=selected[key]["label"], declared_version=version,
                version_scope="Retained method/version declaration; not universal runtime attestation", scores=values,
                output_semantics=semantics, inherited_output_input_identities=True,
                controlled_timing=False, **details))
    result = dict(status="all_benchmark_provenance_register_partial", rows=rows,
        checked_direct_records=list(reader.checked.values()), source=record(__file__),
        helpers=[record(Path(__file__).with_name(name)) for name in (
            "audit_ob_orthofinder_provenance.py", "prepare_ob_candidate_neighborhood.py",
            "summarize_matched_resources.py", "gnu_time_companion.py")],
        complete_transitive_provenance=False, publication_ready=False, scores_recomputed=False,
        controlled_comparative_resources=False,
        limitations=["24 selected rows, not a new raw-input/prediction/native/scoring admission.",
            "Output/input hashes copied from retained admissions are explicitly inherited, not rechecked here.",
            "Unknown is not zero; mixed resource scopes and allocations do not establish speed/memory rankings.",
            "QfO GO/EC/FAS are similarities; the six-score mean is project-defined/secondary.",
            "The completed replacement Threadripper timing panel is separate and is not substituted for historical score-run resources."])
    for pin in reader.checked.values():
        check(pin)
    output.mkdir(parents=True, exist_ok=False)
    (output / "register.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    fields = ["dataset", "method", "version", "wall_seconds", "scope", "peak_process_rss_kib", "memory_scope"]
    with (output / "resources.tsv").open("x") as stream:
        writer = csv.DictWriter(stream, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            for resource in row["resources"] or [dict(scope="Unestablished", measurement={})]:
                values = resource["measurement"]
                writer.writerow(dict(dataset=row["dataset"], method=row["label"], version=row["declared_version"],
                    wall_seconds=values.get("elapsed_seconds"), scope=resource["scope"],
                    peak_process_rss_kib=values.get("max_process_rss_kib"), memory_scope=resource.get("memory_scope", "See retained metadata")))
    lines = ["# All-Tool Benchmark Provenance", "", "Descriptive retained evidence; no comparative speed ranking.", "",
        "| Dataset | Method | Declared version | Resource scopes | Input evidence |", "|---|---|---|---|---|"]
    for row in rows:
        scopes = "; ".join(r["scope"] for r in row["resources"]) or "Full inference unestablished"
        lines.append(f"| {row['dataset']} | {row['label']} | {row['declared_version']} | {scopes} | {row['input_note']} |")
    lines.extend(["", "See register.json for exact scores, prediction/admission pins, commands, CPU/RSS scope and explicit gaps.", "",
                  *["- " + text for text in result["limitations"]]])
    (output / "register.md").write_text("\n".join(lines) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = assemble(args.repo, args.output)
    print(json.dumps({"status": result["status"], "rows": len(result["rows"]), "direct_records": len(result["checked_direct_records"])}))
