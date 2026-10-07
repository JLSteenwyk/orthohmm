"""Test complete accepted-union reconstruction and localize changed VGNC pairs."""

import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools import replay_native_candidate_trace as accepted
from benchmark_tools.export_native_qfo_vgnc_blocks import unique_record
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.validate_native_factorial_outputs import require
from qfo_benchmark.og_to_pairwise import _strip_to_uniprot

ROOT = Path(__file__).resolve().parent.parent
PROTOCOL = ROOT / "benchmark_tools/results/NATIVE_QFO_CANDIDATE_GROUP_TRACE_PROTOCOL_20261006.md"
HEADER = ["protein_left", "protein_right", "baseline_state", "candidate_state",
    "baseline_group_left", "baseline_group_right", "candidate_group", "first_connected_iteration",
    "first_direct_event", "connection_path"]
DECOMPOSITION_SHA = "0c83ab917417ac4e97570b8a68b4bbb0bdfd72eeb301b930bf24446e1a623aea"
READBACK_SHA = "f75cc5775eb63f5a4aa9c72d844782f285d39b00481a380dca6450669c109c57"


def index(groups):
    return {gene: min(group) for group in groups for gene in group}


def localize(initial, final, trace, changes):
    universe, reconstructed = accepted.replay(initial, trace)
    require(reconstructed == final, "Accepted unions do not reconstruct the complete candidate partition")
    require(len(universe) == sum(map(len, final)), "Candidate universe or overlap differs")
    aliases = {}
    for gene in universe:
        accession = _strip_to_uniprot(gene)
        require(accession and accession not in aliases, "Noninjective accession normalization")
        aliases[accession] = gene
    baseline, candidate = index(initial), index(final)
    watched, seen = set(), set()
    for a, b, left, right in changes:
        require(a < b and (a, b) not in seen and a in aliases and b in aliases,
                "Duplicate, unordered or unmapped changed pair")
        require((left, right) in {("FN", "TP"), ("not_scored", "FP")}, "Unsupported changed state")
        seen.add((a, b))
        ga, gb = aliases[a], aliases[b]
        require(baseline[ga] != baseline[gb] and candidate[ga] == candidate[gb],
                "Changed pair not separated in baseline and co-grouped in candidate")
        watched.update((ga, gb))
    require(bool(changes), "No changed pairs")
    histories = {}
    iterations = sorted({r["iteration"] for r in trace})
    for iteration in iterations:
        groups = reconstructed if iteration == iterations[-1] else accepted.replay(initial,
            [r for r in trace if r["iteration"] <= iteration])[1]
        histories[iteration] = {g: min(group) for group in groups for g in group if g in watched}
    events = defaultdict(lambda: defaultdict(set))
    for event, row in enumerate(trace):
        for role in ("source", "target"):
            for gene in row[role + "_genes"]:
                if gene in watched:
                    events[gene][role].add(event)
    ledger, summary = [], Counter()
    for a, b, left, right in changes:
        ga, gb = aliases[a], aliases[b]
        first = next(i for i in iterations if histories[i][ga] == histories[i][gb])
        direct = (events[ga]["source"] & events[gb]["target"]) | (events[ga]["target"] & events[gb]["source"])
        event = min(direct, key=lambda e: (trace[e]["iteration"], e)) if direct else None
        require(event is None or trace[event]["iteration"] == first, "Direct event differs from first connection round")
        path = "direct_cross_endpoint" if event is not None else "transitive_union"
        ledger.append([a, b, left, right, baseline[ga], baseline[gb], candidate[ga], first,
                       event if event is not None else "", path])
        summary[right, first, path] += 1
    return ledger, [dict(candidate_state=s, iteration=i, connection_path=p, pairs=n)
        for (s, i, p), n in sorted(summary.items())], len(universe)


def collect(decomposition_ref, readback_ref):
    checked = []
    def load(ref):
        check(ref); checked.append(ref)
        return json.loads(Path(ref["path"]).read_text())
    decomposition, readback = load(decomposition_ref), load(readback_ref)
    require(decomposition.get("schema") == "native_qfo_candidate_vgnc_decomposition_v1"
        and decomposition.get("status") == "new_candidate_scored_rows_decomposed"
        and readback.get("schema") == "native_qfo_candidate_vgnc_readback_v1"
        and readback.get("status") == "candidate_decomposition_independently_verified"
        and readback.get("report") == decomposition_ref
        and all(decomposition.get(k) is False and readback.get(k) is False for k in
            ("uncertainty_admitted", "new_scoring_or_admission", "publication_ready"))
        and decomposition.get("failed_r1_timing_remains_ineligible") is True
        and readback.get("failed_r1_timing_remains_ineligible") is True, "Changed decomposition/readback scope")
    source_paths = [Path(__file__), PROTOCOL, Path(accepted.__file__),
        ROOT / "benchmark_tools/diagnose_candidate_trace_variation.py",
        ROOT / "benchmark_tools/link_factorial_scaling_resources.py",
        ROOT / "benchmark_tools/score_ygob_groups.py", ROOT / "qfo_benchmark/og_to_pairwise.py",
        ROOT / "benchmark_tools/prepare_ob_candidate_neighborhood.py",
        ROOT / "benchmark_tools/export_native_qfo_vgnc_blocks.py",
        ROOT / "benchmark_tools/probe_dgx_step_separation.py",
        ROOT / "benchmark_tools/validate_native_factorial_outputs.py",
        ROOT / "benchmark_tools/readback_native_qfo_candidate_groups.py",
        ROOT / "benchmark_tools/readback_native_qfo_vgnc_blocks.py"]
    for source in [decomposition["source"], readback["source"], *[record(p) for p in source_paths]]:
        check(source); checked.append(source)
    methods = []
    for index_value, cell, name in ((6, "p0_c0_r0", "baseline"), (8, "p0_c1_r0", "candidate")):
        method = decomposition[name]
        admission = load(method["admission"])
        require(admission.get("schema") == "full_native_factorial_qfo_admission_v1"
            and admission.get("status") == "full_native_factorial_qfo_assessment_admitted"
            and admission.get("accuracy_admitted") is True and admission.get("publication_ready") is False
            and (admission["native_index"], admission["cell"], admission["participant"])
                == (index_value, cell, method["participant"]), "Native admission identity differs")
        conversion_ref = admission["pairs_manifest"]
        require(conversion_ref in admission["checked_records"], "Conversion not originally inventoried")
        conversion = load(conversion_ref)
        require(conversion.get("schema") == "full_native_factorial_qfo_conversion_v1"
            and conversion.get("status") == "full_native_factorial_qfo_pairs_prepared_unscored"
            and conversion.get("conversion_kind") == "group"
            and (conversion["native_index"], conversion["cell"], conversion["native_job_id"])
                == (index_value, cell, admission["native_job_id"]), "Conversion semantics/identity differs")
        review = load(conversion["terminal_review"])
        require(review.get("schema") == "native_factorial_terminal_review_v1" and review.get("status") == "native_success"
            and review.get("native_outputs_validated") is True
            and (review["index"], review["cell"], review["job_id"])
                == (index_value, cell, admission["native_job_id"]), "Native terminal review differs")
        output_ref = review["reviews"]["outputs_or_failure"]
        outputs = load(output_ref)
        require(outputs.get("schema") == "native_factorial_output_review_v1"
            and outputs.get("status") == "native_outputs_validated"
            and outputs.get("native_outputs_validated") is True
            and (outputs["index"], outputs["cell"], outputs["job_id"])
                == (index_value, cell, admission["native_job_id"])
            and conversion["native_input"] in outputs["checked_files"], "Original output review binding differs")
        check(conversion["native_input"]); checked.append(conversion["native_input"])
        methods.append(dict(cell=cell, index=index_value, admission=method["admission"], conversion=conversion_ref,
            terminal_review=conversion["terminal_review"], output_review=output_ref,
            partition=conversion["native_input"], input_genes=outputs["input_genes"], outputs=outputs))
    original_files = methods[1]["outputs"]["checked_files"]
    refs = {name: unique_record(original_files, suffix) for name, suffix in (
        ("trace", "/phylogeny_candidate_merges.json"), ("metrics", "/metrics.json"),
        ("checkpoint", "/phylogeny_candidate_superfamilies.txt"))}
    for ref in refs.values():
        check(ref); checked.append(ref)
    trace, metrics = load(refs["trace"]), load(refs["metrics"])
    candidate_ref = methods[1]["partition"]
    require(all(candidate_ref[k] == refs["checkpoint"][k] for k in ("bytes", "sha256")),
            "Candidate checkpoint differs from converted partition")
    initial = accepted.partition(Path(methods[0]["partition"]["path"]), "space_separated_groups")
    final = accepted.partition(Path(candidate_ref["path"]), "space_separated_groups")
    profile = metrics["metadata"]["phylogeny_candidate_profile"]
    require(metrics.get("status") == "complete" and profile["profile"] == "satellite_v2"
        and initial[0] == final[0] and len(initial[0]) == methods[0]["input_genes"] == methods[1]["input_genes"]
        and len(initial[1]) == profile["seed_families"] and len(final[1]) == profile["candidate_families"]
        and len(trace) == profile["merges"], "Universe or declared candidate counts differ")
    pair_ref = decomposition["transition_table"]
    check(pair_ref); checked.append(pair_ref)
    changes, transition_counts, seen_pairs = [], Counter(), set()
    with Path(pair_ref["path"]).open(newline="") as stream:
        rows = csv.DictReader(stream, delimiter="\t")
        require(rows.fieldnames == ["protein_left", "protein_right", "block_left", "block_right", "p0_c0_r0", "p0_c1_r0"],
                "Unexpected scored transition header")
        for row in rows:
            require(None not in row and all(v is not None for v in row.values()), "Malformed scored transition row")
            pair = (row["protein_left"], row["protein_right"])
            require(pair[0] < pair[1] and pair not in seen_pairs, "Duplicate or unordered scored transition pair")
            seen_pairs.add(pair)
            left, right = row["p0_c0_r0"], row["p0_c1_r0"]
            transition_counts[left, right] += 1
            if left != right:
                changes.append((row["protein_left"], row["protein_right"], left, right))
    require(sum(transition_counts.values()) == decomposition["union_scored_pairs"]
        and [dict(baseline=a, candidate=b, pairs=n) for (a, b), n in sorted(transition_counts.items())]
            == decomposition["transition_counts"], "Changed scored-transition inventory")
    ledger, summary, genes = localize(initial[1], final[1], trace, changes)
    for ref in checked:
        check(ref)
    for method in methods:
        del method["outputs"]
    return dict(schema="native_qfo_candidate_group_trace_v1", status="complete_candidate_groups_and_changed_pairs_reconstructed",
        source=record(__file__), protocol=record(PROTOCOL), decomposition=decomposition_ref, decomposition_readback=readback_ref,
        methods=methods, native_trace=refs["trace"], metrics=refs["metrics"], candidate_checkpoint=refs["checkpoint"],
        scored_pair_table=pair_ref, genes=genes, baseline_groups=len(initial[1]), candidate_groups=len(final[1]),
        accepted_merges=len(trace), changed_pairs=len(ledger), localized_summary=summary, checked_records=checked,
        whole_baseline_replay_matches=True, accuracy_rescored=False, native_inference_reexecuted=False,
        uncertainty_admitted=False, failed_r1_timing_remains_ineligible=True, publication_ready=False,
        limitations=["Development-exposed software path, not true biological orthology or independent validation.",
            "Matching partitions and accepted unions do not establish identical upstream search histories or numeric seed ordering.",
            "Accepted decisions are replayed; eligibility, rejected alternatives and support-score computation are not.",
            "Direct versus transitive attachment is a software graph property, not presence or absence of homology evidence.",
            "Current-byte checks against original inventories are not uninterrupted integrity or repeated transitive admission.",
            "No VGNC CI, model-confidence calibration, default tuning, failed timing repair or isolated-performance claim."]), ledger


def execute(decomposition_ref, readback_ref, output):
    output = Path(output).absolute()
    require(output.resolve() == output and not output.exists() and not output.is_symlink(), "Require fresh direct output")
    output.mkdir(parents=True, exist_ok=False)
    try:
        result, ledger = collect(decomposition_ref, readback_ref)
        table = output / "changed_pair_groups.tsv"
        with table.open("x", newline="") as stream:
            writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
            writer.writerow(HEADER); writer.writerows(ledger)
        result["ledger"] = record(table)
        save(output / "report.json", result)
        return result
    except Exception as error:
        save(output / "failure.json", dict(schema="native_qfo_candidate_group_trace_failure_v1", source=record(__file__),
            decomposition=decomposition_ref, readback=readback_ref, error_type=type(error).__name__, error=str(error),
            automatic_retry=False, accuracy_rescored=False, native_inference_reexecuted=False, publication_ready=False))
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--decomposition", type=Path, default=ROOT / "benchmark_tools/results/native_qfo_candidate_vgnc_20261006_v1/report.json")
    parser.add_argument("--decomposition-sha256", default=DECOMPOSITION_SHA)
    parser.add_argument("--readback", type=Path, default=ROOT / "benchmark_tools/results/native_qfo_candidate_vgnc_readback_20261006_v1.json")
    parser.add_argument("--readback-sha256", default=READBACK_SHA)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    refs = [record(args.decomposition), record(args.readback)]
    require([r["sha256"] for r in refs] == [args.decomposition_sha256, args.readback_sha256], "Selected evidence changed")
    result = execute(*refs, args.output)
    print(json.dumps({k: result[k] for k in ("status", "genes", "baseline_groups", "candidate_groups", "accepted_merges", "changed_pairs", "localized_summary")}))
