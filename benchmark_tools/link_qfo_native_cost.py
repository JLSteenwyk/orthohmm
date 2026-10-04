"""Associate retained full native QfO inference with a verified factorial partition."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmark_tools.assemble_benchmark_provenance import Reader
from benchmark_tools.audit_ob_orthofinder_provenance import log_command
from benchmark_tools.collect_factorial_resources import FIXED_INPUTS, finite, require
from benchmark_tools.link_factorial_scaling_resources import COMMON, STAGES, partition
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.summarize_matched_resources import verbose_time


CORE = "7f3a9e40dd7e79f842cc2c11fb8b548f9a802806"
CELL = "p1_c0_r0"
FIXED = {
    "native": ("benchmark_tools/results/qfo_corrected_high_sensitivity_admission_21720.json",
               "845874021057da54cd61c20763b00280c1a5569c5f451676e52f16043dacccff"),
    "prepared": FIXED_INPUTS["qfo_preparation"],
    "selected": ("benchmark_tools/results/qfo_orthohmm_stage_metadata_20261004/register.json",
                 "a1cf40e172614f7c9a5e7906cb4121dce17ac7f6c8ff516977a421de05ac867a"),
    "costs": ("benchmark_tools/results/factorial_retained_resources_20261004/resources.json",
              "4c7801197a9bc89927964f952d03325187341a0ce2fe2cc5f6c4aa7700e581a6"),
}


def unique_pin(records, path):
    found = {json.dumps(r, sort_keys=True): r for r in records if r["path"] == str(path)}
    require(len(found) == 1, "Missing/conflicting admitted record: " + str(path))
    return next(iter(found.values()))


def validate_execution(native, execution, plan, metrics, timing, timed_command):
    require(native["status"] == "corrected_high_sensitivity_native_evidence_admitted",
            "Native evidence not admitted")
    scheduler = native["scheduler"]
    require(all(scheduler[k] == v for k, v in {
        "JobIDRaw": "21707", "State": "COMPLETED", "ExitCode": "0:0",
        "NodeList": "bizon", "AllocCPUS": "32"}.items()), "Wrong native scheduler outcome")
    require(execution["job_id"] == scheduler["JobIDRaw"]
            and execution["method"] == "orthohmm_high_sensitivity"
            and execution["status"] == "process_succeeded_pending_native_admission"
            and execution["exit_code"] == 0, "Native execution differs from admission")
    config = plan["methods"]["orthohmm_high_sensitivity"]
    argv = config["native_argv"]
    require(config["search_reuse"] is False and execution["native_argv"] == argv == timed_command,
            "Wrong command or search reuse")
    require(argv[1:3] == ["-m", "orthohmm"]
            and metrics["command"] == [argv[0], str(Path(config["cwd"]) / "orthohmm/__main__.py"), *argv[3:]]
            and metrics["cwd"] == execution["cwd"] == config["cwd"], "Metrics command differs")
    require(argv[argv.index("--refinement_profile") + 1] == "default"
            and argv[argv.index("--stop") + 1] == "infer"
            and "--start" not in argv and "--phylogeny" not in argv,
            "Wrong fresh inference mode")
    require(metrics["status"] == "complete" and set(metrics["stages"]) == STAGES,
            "Missing or unexpected native stages")
    require(all(metrics["metadata"].get(k) == v for k, v in COMMON.items()),
            "Native scientific settings differ")
    require(metrics["counts"]["species"] == 78 and metrics["counts"]["genes"] == 984137
            and metrics["counts"]["orthogroups"] == 391908, "Wrong corrected input/output size")
    require(metrics["rss_measurement"] == "sampled_sum_of_linux_proc_tree_rss",
            "Unexpected historical memory convention")
    for key in ("wall_s", "user_cpu_s", "system_cpu_s", "peak_process_tree_rss_bytes"):
        finite(metrics[key], key, positive=True)
    require(type(metrics["peak_process_tree_rss_bytes"]) is int, "Noninteger memory bytes")
    for stage in metrics["stages"].values():
        for key in ("wall_s", "user_cpu_s", "system_cpu_s", "peak_process_tree_rss_bytes"):
            finite(stage[key], key, positive=key in {"wall_s", "peak_process_tree_rss_bytes"})
    require(timing["exit_status"] == 0, "Failed GNU-time command")


def validate_join(prepared, replay, replay_plan, native_pin, execution, selected, costs):
    require(prepared["status"] == "corrected_qfo_four_candidate_arms_prepared_unscored"
            and replay["status"] == "corrected_checked_replay_admitted", "Candidate/replay not admitted")
    require(prepared["plan"] == replay["plan"] and replay_plan["admission"] == native_pin
            and replay_plan["primary_plan"] == execution["manifest"], "Native/replay plan chain differs")
    require(prepared["input_fastas"] == replay_plan["input_fastas"] == selected["input_records"]
            and len(prepared["input_fastas"]) == 78, "Corrected FASTA manifests differ")
    arm = prepared["candidate_arms"]["p1_c0"]
    require(arm["candidate_expansion"] is False
            and all(arm["seed_partition"][k] == arm["candidate_partition"][k]
                    for k in ("bytes", "sha256")), "Candidate expansion present in C-off arm")
    require(selected["dataset"] == "QfO" and selected["key"] == "orthohmm_high_sensitivity"
            and selected["qfo_stage_provenance"]["cell"] == CELL
            and selected["qfo_stage_provenance"]["candidate_partition"] == arm["candidate_partition"]
            and selected["qfo_stage_provenance"]["replay_admission"] == prepared["admission"],
            "Selected score uses different arm or replay")
    seed = [s["output"] for s in replay["coverage"] if s["label"] == "strict_profiles_refined"]
    require(seed == [arm["seed_partition"]], "Wrong replay seed")
    require(costs["dataset"] == "Corrected QfO" and costs["cell"] == CELL
            and all(costs[k] is None for k in (
                "full_pipeline_wall_s", "full_pipeline_cpu_s", "full_pipeline_peak_memory_bytes")),
            "Original cached-execution full costs must remain unavailable")


def collect(root):
    root, reader = Path(root).resolve(), Reader()
    pins, docs = {}, {}
    for name, (path, digest) in FIXED.items():
        pins[name] = record(root / path)
        require(pins[name]["sha256"] == digest, "Fixed input changed: " + name)
        docs[name] = reader.read(pins[name])
    native, prepared = docs["native"], docs["prepared"]
    admitted = native["checked_records"]
    output = Path(native["content"]["native_groups"]["path"])
    execution = reader.read(unique_pin(admitted, output.parent.parent / "execution/orthohmm_high_sensitivity/status.json"))
    plan = reader.read(execution["manifest"])
    require(execution["manifest"] == unique_pin(admitted, execution["manifest"]["path"]),
            "Primary plan not admitted")
    metrics = reader.read(native["content"]["metrics"])
    require(native["content"]["metrics"] == unique_pin(admitted, plan["methods"]["orthohmm_high_sensitivity"]["metrics"]),
            "Metrics not admitted")
    require(execution["timing"] == unique_pin(admitted, execution["timing"]["path"]), "Timing not admitted")
    text = reader.text(execution["timing"])
    timing = verbose_time(text)
    validate_execution(native, execution, plan, metrics, timing, log_command(text, "Command being timed: "))
    replay = reader.read(prepared["admission"])
    replay_plan = reader.read(prepared["plan"])
    replay_native = reader.read(replay_plan["admission"])
    require(replay_native == native, "Different native admission content")
    require(native["checkpoint_manifest"] == replay_plan["checkpoint_manifest"]
            == prepared["numeric_checkpoint"]["manifest"], "Different upstream numeric checkpoint")
    reader.read(native["checkpoint_manifest"])
    selected = [r for r in docs["selected"]["rows"] if r["dataset"] == "QfO"
                and r["key"] == "orthohmm_high_sensitivity"]
    costs = [r for r in docs["costs"]["rows"] if r["dataset"] == "Corrected QfO" and r["cell"] == CELL]
    require(len(selected) == len(costs) == 1, "Missing/duplicate selected configuration")
    validate_join(prepared, replay, replay_plan, replay_plan["admission"], execution, selected[0], costs[0])
    baseline = reader.read(plan["baseline"])
    require(baseline["core_commit"] == prepared["runtime_before"]["core_commit"] == CORE
            and prepared["runtime_before"] == prepared["runtime_after"]
            and baseline["core_root"] == execution["cwd"], "Frozen core binding differs")
    core_records = []
    for ref in baseline["core_sources"]:
        pin = {"path": ref["absolute_path"], "bytes": ref["bytes"], "sha256": ref["sha256"]}
        check(pin)
        reader.checked[pin["path"]] = pin
        core_records.append(pin)
    require(len(core_records) == 40 and len({p["path"] for p in core_records}) == 40,
            "Wrong frozen source inventory")
    universe = set()
    for pin in prepared["input_fastas"]:
        require(pin == unique_pin(admitted, pin["path"]) and pin in plan["inputs"], "Input not admitted/planned")
        for line in reader.text(pin).splitlines():
            if line.startswith(">"):
                fields = line[1:].split()
                require(bool(fields) and fields[0] not in universe, "Empty/duplicate FASTA identifier")
                universe.add(fields[0])
    require(len(universe) == 984137, "Wrong corrected FASTA universe")
    candidate_pin = prepared["candidate_arms"]["p1_c0"]["candidate_partition"]
    native_pin = native["content"]["native_groups"]
    for pin in (native_pin, candidate_pin):
        check(pin)
        reader.checked[pin["path"]] = pin
    left = partition(Path(native_pin["path"]), "named_groups")
    right = partition(Path(candidate_pin["path"]), "space_separated_groups")
    require(left[0] == right[0] == universe and len(left[1]) == len(right[1]) == 391908,
            "Partitions differ from input/output inventory")
    require(left[1] == right[1], "Native/factorial partition differs; refuse exact-output association")
    for pin in reader.checked.values():
        check(pin)
    return {
        "schema": "qfo_native_configuration_cost_v1", "status": "retained_native_cost_partition_linked",
        "dataset": "Corrected QfO", "cell": CELL, "job_id": 21707, "repeats": 1,
        "frozen_core_commit": CORE, "native_python": metrics["python"],
        "sources": pins, "checked_records": sorted(reader.checked.values(), key=lambda r: r["path"]),
        "input_fastas": prepared["input_fastas"], "partition": {
            "native": native_pin, "factorial": candidate_pin, "genes": len(universe),
            "groups": len(right[1]), "partition_equal": True},
        "native_metrics": native["content"]["metrics"], "scheduler": native["scheduler"],
        "native_argv": execution["native_argv"], "metrics": {k: metrics[k] for k in (
            "wall_s", "user_cpu_s", "system_cpu_s", "peak_process_tree_rss_bytes", "rss_measurement")},
        "recorded_stages": metrics["stages"],
        "companion": {"evidence": execution["timing"], "measurement": timing},
        "original_cached_full_costs": {k: costs[0][k] for k in (
            "full_pipeline_wall_s", "full_pipeline_cpu_s", "full_pipeline_peak_memory_bytes")},
        "scoring_repeated": False, "native_inference_repeated": False,
        "controlled_comparative_resources": False, "publication_ready": False,
        "limitations": [
            "One retained historical shared-host native observation, not a new repeat or the original cached factorial execution's cost.",
            "Native metrics cover initial search through group materialization; exclude preparation, parent validation, conversion and scoring. Stage measurements are recorded, not summed into synthetic totals.",
            "Sampled summed process-tree RSS can miss short-lived peaks and double-count shared pages. GNU-time maximum process RSS is a different scope; neither is a lifetime cgroup memory peak.",
            "Retained GNU-time launch-to-exit and native internal intervals differ; keep their values/scopes separate. No background-workload series or isolated/corrected estimate is imputed.",
            "Exact FASTA bytes, frozen source bytes, command settings and whole partition are checked. Reuse historical admission; this is not a new historical process trace, search-completeness or pairwise-output audit.",
            "Do not pool with the later private-runtime Threadripper scaling panel or infer isolated tool rankings, causal component overhead or matched full costs for unmeasured configurations.",
        ],
    }


def render(report):
    m, c = report["metrics"], report["companion"]["measurement"]
    return "\n".join([
        "# Retained Native QfO Cost Association", "",
        "Historical shared host; unknown, potentially method-dependent contention. One observation.", "",
        "| Dataset / cell | Native wall, s | User CPU, s | System CPU, s | Sampled tree RSS, GiB |",
        "| --- | ---: | ---: | ---: | ---: |",
        "| Corrected QfO / %s | %.6f | %.6f | %.6f | %.6f |" %
        (CELL, m["wall_s"], m["user_cpu_s"], m["system_cpu_s"], m["peak_process_tree_rss_bytes"] / 2**30),
        "", "78 FASTAs, 984,137 genes and 391,908 groups checked; native/factorial partitions equal.",
        "Original cached factorial full costs remain unavailable; this is a separate native-run association.",
        "", "GNU-time companion: %.3f s launch-to-exit, %.2f user CPU s, %.2f system CPU s, %d KiB maximum process RSS."
        % (c["elapsed_seconds"], c["user_seconds"], c["system_seconds"], c["max_process_rss_kib"]),
        "These are distinct intervals/memory scopes, not additional independent repeats.",
        "", "## Limitations", "", *["- " + x for x in report["limitations"]], "",
    ])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    args = parser.parse_args()
    require(not args.output_directory.exists(), "Refusing existing output directory")
    report = collect(args.root)
    args.output_directory.mkdir(parents=True, exist_ok=False)
    (args.output_directory / "association.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    (args.output_directory / "association.md").write_text(render(report))


if __name__ == "__main__":
    main()
