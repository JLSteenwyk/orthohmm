"""Verify and rescore preserved OrthoBench profile/refinement controls."""

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.score_orthobench_partition import read_named_groups, score_partition


STAGES = ("multipass", "multipass_refined", "strict_profiles", "strict_profiles_refined")


def verify_file(path, expected):
    actual = file_provenance(path)
    if any(actual[k] != expected[k] for k in ("bytes", "sha256")):
        raise ValueError(f"Changed artifact: {path}")
    return actual


def read_partition(path, universe):
    groups = []
    seen = set()
    for line in path.read_text().splitlines():
        genes = line.split()
        if not genes:
            continue
        group = set(genes)
        if len(group) != len(genes) or group & seen or not group <= universe:
            raise ValueError(f"Invalid partition membership: {path}")
        groups.append(group)
        seen.update(group)
    if seen != universe:
        raise ValueError(f"Partition does not cover input universe: {path}")
    return groups


def audit(replay_path, fresh_metrics_path, fresh_output, fasta_directory, refogs):
    replay = json.loads(replay_path.read_text())
    fresh = json.loads(fresh_metrics_path.read_text())
    if fresh["status"] != "complete" or fresh["harness"]["exit_code"] != 0:
        raise ValueError("Fresh production run did not complete successfully")
    stage_index = {r["label"]: r for r in replay["stages"]}
    if set(stage_index) != set(STAGES) or len(replay["stages"]) != len(STAGES):
        raise ValueError("Unexpected historical ablation panel")
    parameters = replay["parameters"]
    expected = {"accuracy_profile": "high_sensitivity", "cpm_resolution": 0.1,
                "leiden_seed": 4, "matrix": "BLOSUM62", "profile_expansion": True}
    if any(parameters.get(k) != v for k, v in expected.items()):
        raise ValueError("Unexpected replay settings")
    input_records = []
    universe = set()
    input_names = {r["path"] for r in fresh["harness"]["input_manifest"]}
    observed_names = {p.name for p in fasta_directory.iterdir() if p.suffix.lower() in {".fa", ".faa", ".fasta", ".fsa"}}
    if observed_names != input_names:
        raise ValueError("FASTA directory differs from fresh-run manifest")
    for item in fresh["harness"]["input_manifest"]:
        path = fasta_directory / item["path"]
        input_records.append(verify_file(path, item))
        for record in SeqIO.parse(path, "fasta"):
            if record.id in universe:
                raise ValueError("Duplicate input gene")
            universe.add(record.id)
    cache = verify_file(Path(replay["input"]["path"]), replay["input"])
    fresh_relative = "orthohmm_working_res/orthohmm_edges_clustered.txt"
    expected_outputs = [r for r in fresh["harness"]["output_manifest"] if r["path"] == fresh_relative]
    if len(expected_outputs) != 1:
        raise ValueError("Fresh partition manifest missing or ambiguous")
    fresh_record = verify_file(fresh_output / fresh_relative, expected_outputs[0])
    partitions = {}
    provenance = {}
    for label in STAGES:
        item = stage_index[label]["output"]
        path = Path(item["path"])
        provenance[label] = verify_file(path, item)
        partitions[label] = read_partition(path, universe)
    if provenance["strict_profiles_refined"]["sha256"] != fresh_record["sha256"]:
        raise ValueError("Replay does not exactly reproduce corrected fresh partition")
    fresh_partition = read_partition(fresh_output / fresh_relative, universe)
    if {frozenset(g) for g in fresh_partition} != {frozenset(g) for g in partitions["strict_profiles_refined"]}:
        raise ValueError("Partition equality check failed")
    names = sorted(p.name for p in refogs.glob("RefOG*.txt"))
    if len(names) != 70:
        raise ValueError("Expected the complete 70-RefOG OrthoBench panel")
    references = read_named_groups(refogs, names)
    uncertain_paths = [refogs / "low_certainty_assignments" / name for name in names]
    uncertain = {p.name: set(p.read_text().split()) if p.exists() else set() for p in uncertain_paths}
    scores = {label: score_partition(partitions[label], references, uncertain) for label in STAGES}
    effects = {}
    for name, treated, control in (
        ("profile_without_refinement", "strict_profiles", "multipass"),
        ("profile_with_refinement", "strict_profiles_refined", "multipass_refined"),
        ("refinement_without_profile", "multipass_refined", "multipass"),
        ("refinement_with_profile", "strict_profiles_refined", "strict_profiles"),
    ):
        effects[name] = {"treated": treated, "control": control,
                         "difference_percentage_points": {k: scores[treated][k] - scores[control][k]
                                                          for k in ("f_score", "precision", "recall")}}
    return {"schema_version": 1, "generated_at": datetime.now(timezone.utc).isoformat(),
            "historical_evidence_only": True, "fresh_partition_byte_identical": True,
            "genes": len(universe), "scores": scores, "effects": effects,
            "cluster_counts": {k: len(v) for k, v in partitions.items()},
            "parameters": parameters, "historical_replay_git": replay["git"],
            "historical_fresh_harness_git": fresh["harness"]["git_commit"],
            "historical_replay_resources": {k: replay.get(k) for k in ("wall_s", "timings", "peak_process_rss_gib")},
            "fresh_wall_s": fresh["wall_s"], "runtime_comparison_controlled": False,
            "inputs": {"replay": file_provenance(replay_path), "fresh_metrics": file_provenance(fresh_metrics_path),
                       "cache": cache, "fastas": input_records, "fresh_partition": fresh_record,
                       "stage_partitions": provenance,
                       "references": [file_provenance(refogs / n) for n in names],
                       "uncertain": [file_provenance(p) for p in uncertain_paths if p.exists()]},
            "source": file_provenance(Path(__file__)),
            "scorer_source": file_provenance(Path(__file__).with_name("score_orthobench_partition.py")),
            "limitations": [
                "Historical development-exposed controls; descriptive effects, not prospective confirmation.",
                "Profile omission retains the HMM-based initial search; this is not an HMM-free comparison.",
                "Byte equality validates this historical endpoint, not all intermediate stages or future source versions.",
                "The cache hash is verified without deserializing or rechecking each normalized search hit.",
                "Current-source replay, full expansion/reconciliation factorial, and matched sequence-search control remain required.",
                "Historical source dirtiness is preserved in provenance; no clean-source equivalence is inferred.",
                "Replay timings are cumulative incremental costs, not independently timed ablation arms.",
            ]}


def render_report(result):
    lines = ["# Historical Profile And Refinement Ablation", "",
             "Development-exposed, descriptive controls; not the completed publication factorial.", "",
             f"All four partitions cover {result['genes']:,} input genes. The final replay partition",
             "is byte-identical to the corrected fresh high-sensitivity production output.", "",
             "| Stage | Groups | F1 (%) | Precision (%) | Recall (%) |",
             "| --- | ---: | ---: | ---: | ---: |"]
    for label in STAGES:
        scores = result["scores"][label]
        lines.append(f"| {label} | {result['cluster_counts'][label]} | "
                     + " | ".join(f"{scores[k]:.6f}" for k in ("f_score", "precision", "recall")) + " |")
    lines += ["", "| Descriptive effect | F1 difference (pp) | Precision difference (pp) | Recall difference (pp) |",
              "| --- | ---: | ---: | ---: |"]
    for name, effect in result["effects"].items():
        lines.append("| " + name + " | " + " | ".join(
            f"{effect['difference_percentage_points'][k]:+.6f}" for k in ("f_score", "precision", "recall")) + " |")
    lines += ["", "## Interpretation", "",
              "Profile expansion has a modest positive observed effect in this historical panel.",
              "The larger refinement effect must not be attributed to profile HMM expansion.",
              "All arms retain the HMM-based initial search. No significance or independent",
              "generalization claim is made from these descriptive differences.", "",
              "## Limitations", "", *["- " + x for x in result["limitations"]]]
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("replay", "fresh-metrics", "fresh-output", "fasta-directory", "refogs", "json"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--markdown", type=Path)
    args = parser.parse_args()
    if args.json.exists():
        raise ValueError("Refusing to overwrite ablation audit")
    if args.markdown is not None and args.markdown.exists():
        raise ValueError("Refusing to overwrite ablation report")
    result = audit(args.replay, args.fresh_metrics, args.fresh_output, args.fasta_directory, args.refogs)
    result["command"] = [sys.executable, *sys.argv]
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    if args.markdown is not None:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(render_report(result))
    print(json.dumps({"scores": {k: {m: v[m] for m in ("f_score", "precision", "recall")} for k, v in result["scores"].items()},
                      "effects": result["effects"]}, indent=2))


if __name__ == "__main__":
    main()
