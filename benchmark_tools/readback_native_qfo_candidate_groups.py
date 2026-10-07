"""Independent stdlib graph traversal of accepted native candidate unions."""

import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools import readback_native_qfo_vgnc_blocks as evidence

ROOT = Path(__file__).resolve().parent.parent
HEADER = ["protein_left", "protein_right", "baseline_state", "candidate_state",
    "baseline_group_left", "baseline_group_right", "candidate_group", "first_connected_iteration",
    "first_direct_event", "connection_path"]


def partitions(path):
    groups, genes = {}, set()
    with Path(path).open() as stream:
        for line in stream:
            tokens = line.split()
            evidence.need(tokens and len(set(tokens)) == len(tokens), "Empty or duplicate group")
            members = frozenset(tokens)
            evidence.need(not genes.intersection(members), "Overlapping partition")
            groups[min(members)] = members
            genes.update(members)
    evidence.need(groups, "Empty partition")
    return genes, groups


def traverse(groups, trace, watched):
    """Merge connected components by graph traversal, without union-find."""
    evidence.need(isinstance(trace, list) and trace, "Empty or malformed trace")
    rounds, semantic = defaultdict(list), set()
    for event, row in enumerate(trace):
        iteration = row["iteration"]
        evidence.need(type(iteration) is int and iteration in (0, 1), "Wrong iteration")
        for role in ("source", "target"):
            tokens = row[role + "_genes"]
            evidence.need(isinstance(tokens, list) and tokens
                and all(isinstance(g, str) and g and g.strip() == g for g in tokens)
                and len(set(tokens)) == len(tokens)
                and type(row[role + "_size"]) is int and row[role + "_size"] == len(tokens)
                and type(row[role + "_cluster"]) is int and row[role + "_cluster"] >= 0,
                "Malformed trace endpoint")
        a, b = frozenset(row["source_genes"]), frozenset(row["target_genes"])
        evidence.need(not a.intersection(b), "Overlapping trace endpoints")
        signature = (iteration, a, b)
        evidence.need(signature not in semantic, "Duplicate accepted event")
        semantic.add(signature)
        rounds[iteration].append((event, row, a, b))
    current, histories = dict(groups), {}
    for iteration, rows in sorted(rounds.items()):
        labels, graph = {}, defaultdict(set)
        for _, row, a, b in rows:
            for role, endpoint in (("source", a), ("target", b)):
                key = min(endpoint)
                evidence.need(current.get(key) == endpoint, "Endpoint is not a complete round-start group")
                label = row[role + "_cluster"]
                evidence.need(label not in labels or labels[label] == key, "Conflicting cluster label")
                labels[label] = key
            ka, kb = min(a), min(b)
            graph[ka].add(kb); graph[kb].add(ka)
        merged, visited = {}, set()
        for key in current:
            if key in visited:
                continue
            pending, members = [key], set()
            while pending:
                node = pending.pop()
                if node in visited:
                    continue
                visited.add(node)
                members.update(current[node])
                pending.extend(graph[node] - visited)
            merged[min(members)] = frozenset(members)
        evidence.need(len(current) - len(merged) == len(rows), "Cyclic or redundant accepted unions")
        current = merged
        histories[iteration] = {g: k for k, members in current.items() for g in members if g in watched}
    return current, histories


def localize(initial, final, trace, changes):
    initial_owners = {g: k for k, members in initial.items() for g in members}
    final_owners = {g: k for k, members in final.items() for g in members}
    evidence.need(set(initial_owners) == set(final_owners), "Different complete universes")
    aliases = {}
    for gene in initial_owners:
        tokens = gene.split("|")
        accession = tokens[1] if len(tokens) >= 2 and tokens[1] else gene
        evidence.need(accession not in aliases, "Noninjective accession normalization")
        aliases[accession] = gene
    watched, seen = set(), set()
    for a, b, left, right in changes:
        evidence.need(a < b and (a, b) not in seen and a in aliases and b in aliases,
                      "Duplicate, unordered or unmapped changed pair")
        evidence.need((left, right) in (("FN", "TP"), ("not_scored", "FP")), "Unsupported changed state")
        ga, gb = aliases[a], aliases[b]
        evidence.need(initial_owners[ga] != initial_owners[gb] and final_owners[ga] == final_owners[gb],
                      "Changed pair does not acquire shared membership")
        watched.update((ga, gb)); seen.add((a, b))
    evidence.need(changes, "No changed pairs")
    reconstructed, histories = traverse(initial, trace, watched)
    evidence.need(reconstructed == final, "Complete candidate partition differs from graph traversal")
    source_events, target_events = defaultdict(set), defaultdict(set)
    for event, row in enumerate(trace):
        for gene in watched.intersection(row["source_genes"]):
            source_events[gene].add(event)
        for gene in watched.intersection(row["target_genes"]):
            target_events[gene].add(event)
    result, counts = [], Counter()
    for a, b, left, right in changes:
        ga, gb = aliases[a], aliases[b]
        connected = [i for i, owners in sorted(histories.items()) if owners[ga] == owners[gb]]
        evidence.need(connected, "Pair never connected")
        first = connected[0]
        events = (source_events[ga] & target_events[gb]) | (source_events[gb] & target_events[ga])
        direct = sorted(events, key=lambda e: (trace[e]["iteration"], e))
        event = direct[0] if direct else ""
        evidence.need(not direct or trace[event]["iteration"] == first, "Direct-event iteration differs")
        path = "direct_cross_endpoint" if direct else "transitive_union"
        result.append([a, b, left, right, initial_owners[ga], initial_owners[gb], final_owners[ga], first, event, path])
        counts[right, first, path] += 1
    return result, [dict(candidate_state=s, iteration=i, connection_path=p, pairs=n)
        for (s, i, p), n in sorted(counts.items())]


def review(path, digest):
    ref = evidence.identity(path)
    evidence.need(ref["sha256"] == digest, "Changed primary report")
    report = evidence.read_json(ref)
    evidence.need(report.get("schema") == "native_qfo_candidate_group_trace_v1"
        and report.get("status") == "complete_candidate_groups_and_changed_pairs_reconstructed"
        and all(report.get(k) is False for k in ("accuracy_rescored", "native_inference_reexecuted",
            "uncertainty_admitted", "publication_ready"))
        and report.get("failed_r1_timing_remains_ineligible") is True
        and report.get("whole_baseline_replay_matches") is True, "Changed diagnostic scope")
    checked = report["checked_records"]
    for item in checked:
        evidence.verify(item)
    for key in ("source", "protocol", "decomposition", "decomposition_readback", "native_trace", "metrics",
                "candidate_checkpoint", "scored_pair_table"):
        evidence.need(report[key] in checked, "Unbound report artifact: " + key)
    evidence.need(report["source"] == evidence.identity(Path(__file__).with_name("trace_native_qfo_candidate_groups.py"))
        and report["protocol"] == evidence.identity(ROOT / "benchmark_tools/results/NATIVE_QFO_CANDIDATE_GROUP_TRACE_PROTOCOL_20261006.md")
        and evidence.identity(__file__) in checked and evidence.identity(evidence.__file__) in checked,
        "Unbound primary/protocol/reader sources")
    for name in ("replay_native_candidate_trace.py", "diagnose_candidate_trace_variation.py",
                 "link_factorial_scaling_resources.py", "score_ygob_groups.py",
                 "prepare_ob_candidate_neighborhood.py", "export_native_qfo_vgnc_blocks.py",
                 "probe_dgx_step_separation.py", "validate_native_factorial_outputs.py"):
        evidence.need(evidence.identity(Path(__file__).with_name(name)) in checked, "Unbound source: " + name)
    evidence.need(evidence.identity(ROOT / "qfo_benchmark/og_to_pairwise.py") in checked, "Unbound normalization source")
    decomposition = evidence.read_json(report["decomposition"])
    readback = evidence.read_json(report["decomposition_readback"])
    evidence.need(decomposition.get("schema") == "native_qfo_candidate_vgnc_decomposition_v1"
        and decomposition.get("status") == "new_candidate_scored_rows_decomposed"
        and readback.get("schema") == "native_qfo_candidate_vgnc_readback_v1"
        and readback.get("status") == "candidate_decomposition_independently_verified"
        and readback.get("report") == report["decomposition"]
        and all(r.get(k) is False for r in (decomposition, readback)
            for k in ("uncertainty_admitted", "new_scoring_or_admission", "publication_ready"))
        and all(r.get("failed_r1_timing_remains_ineligible") is True for r in (decomposition, readback))
        and report["scored_pair_table"] == decomposition["transition_table"]
        and decomposition["source"] in checked and readback["source"] in checked,
        "Prior decomposition/readback linkage differs")
    evidence.need(len(report["methods"]) == 2, "Wrong method cohort")
    outputs = []
    for method, expected, name in zip(report["methods"], ((6, "p0_c0_r0"), (8, "p0_c1_r0")), ("baseline", "candidate")):
        evidence.need((method["index"], method["cell"]) == expected
            and method["admission"] == decomposition[name]["admission"], "Method identity differs")
        for key in ("admission", "conversion", "terminal_review", "output_review", "partition"):
            evidence.need(method[key] in checked, "Method reference not checked: " + key)
        admission = evidence.read_json(method["admission"])
        evidence.need(admission.get("schema") == "full_native_factorial_qfo_admission_v1"
            and admission.get("status") == "full_native_factorial_qfo_assessment_admitted"
            and admission.get("accuracy_admitted") is True and admission.get("publication_ready") is False
            and (admission["native_index"], admission["cell"]) == expected
            and admission["participant"] == decomposition[name]["participant"]
            and admission["pairs_manifest"] == method["conversion"]
            and method["conversion"] in admission["checked_records"], "Admission/conversion binding differs")
        conversion = evidence.read_json(method["conversion"])
        terminal = evidence.read_json(method["terminal_review"])
        output = evidence.read_json(method["output_review"])
        job = admission["native_job_id"]
        evidence.need(conversion.get("schema") == "full_native_factorial_qfo_conversion_v1"
            and conversion.get("status") == "full_native_factorial_qfo_pairs_prepared_unscored"
            and conversion.get("conversion_kind") == "group"
            and (conversion["native_index"], conversion["cell"], conversion["native_job_id"]) == (*expected, job)
            and conversion["terminal_review"] == method["terminal_review"]
            and conversion["native_input"] == method["partition"], "Conversion semantics differ")
        evidence.need(terminal.get("schema") == "native_factorial_terminal_review_v1"
            and terminal.get("status") == "native_success" and terminal.get("native_outputs_validated") is True
            and (terminal["index"], terminal["cell"], terminal["job_id"]) == (*expected, job)
            and terminal["reviews"]["outputs_or_failure"] == method["output_review"], "Terminal review differs")
        evidence.need(output.get("schema") == "native_factorial_output_review_v1"
            and output.get("status") == "native_outputs_validated" and output.get("native_outputs_validated") is True
            and (output["index"], output["cell"], output["job_id"]) == (*expected, job)
            and method["partition"] in output["checked_files"]
            and method["input_genes"] == output["input_genes"], "Native output inventory differs")
        outputs.append(output)
    for key, suffix in (("native_trace", "/phylogeny_candidate_merges.json"), ("metrics", "/metrics.json"),
                        ("candidate_checkpoint", "/phylogeny_candidate_superfamilies.txt")):
        matches = {(r["path"], r["bytes"], r["sha256"]) for r in outputs[1]["checked_files"] if r["path"].endswith(suffix)}
        expected = report[key]
        evidence.need(matches == {(expected["path"], expected["bytes"], expected["sha256"])},
                      "Ambiguous or unbound original candidate artifact")
    evidence.need(all(report["candidate_checkpoint"][k] == report["methods"][1]["partition"][k]
        for k in ("bytes", "sha256")), "Candidate checkpoint content differs")
    initial_genes, initial = partitions(report["methods"][0]["partition"]["path"])
    final_genes, final = partitions(report["methods"][1]["partition"]["path"])
    trace, metrics = evidence.read_json(report["native_trace"]), evidence.read_json(report["metrics"])
    profile = metrics["metadata"]["phylogeny_candidate_profile"]
    evidence.need(metrics.get("status") == "complete" and profile["profile"] == "satellite_v2"
        and initial_genes == final_genes and len(initial_genes) == report["genes"]
            == outputs[0]["input_genes"] == outputs[1]["input_genes"]
        and len(initial) == report["baseline_groups"] == profile["seed_families"]
        and len(final) == report["candidate_groups"] == profile["candidate_families"]
        and len(trace) == report["accepted_merges"] == profile["merges"], "Declared complete counts differ")
    changes, counts, seen = [], Counter(), set()
    with Path(report["scored_pair_table"]["path"]).open(newline="") as stream:
        rows = csv.reader(stream, delimiter="\t")
        evidence.need(next(rows) == ["protein_left", "protein_right", "block_left", "block_right", "p0_c0_r0", "p0_c1_r0"],
                      "Scored transition header differs")
        for row in rows:
            evidence.need(len(row) == 6 and row[0] < row[1] and tuple(row[:2]) not in seen,
                          "Malformed or duplicate transition pair")
            seen.add(tuple(row[:2])); counts[row[4], row[5]] += 1
            if row[4] != row[5]:
                changes.append((row[0], row[1], row[4], row[5]))
    evidence.need(len(seen) == decomposition["union_scored_pairs"] and decomposition["transition_counts"]
        == [dict(baseline=a, candidate=b, pairs=n) for (a, b), n in sorted(counts.items())], "Transition inventory differs")
    ledger, summary = localize(initial, final, trace, changes)
    evidence.need(len(ledger) == report["changed_pairs"] and summary == report["localized_summary"], "Localized summary differs")
    evidence.table(report["ledger"], HEADER, ledger)
    for item in checked:
        evidence.verify(item)
    evidence.verify(ref)
    return dict(schema="native_qfo_candidate_group_trace_readback_v1", status="complete_graph_and_pair_localization_independently_verified",
        report=ref, source=evidence.identity(__file__), identity_kernel_source=evidence.identity(evidence.__file__),
        genes=len(initial_genes), baseline_groups=len(initial), candidate_groups=len(final), accepted_merges=len(trace),
        changed_pairs=len(ledger), localized_summary=summary, checked_input_records=len(checked),
        primary_or_accepted_replay_imported=False, full_numeric_support_recomputed=False,
        whole_candidate_partition_reconstructed=True, accuracy_rescored=False, native_inference_reexecuted=False,
        uncertainty_admitted=False, failed_r1_timing_remains_ineligible=True, publication_ready=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--report-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    evidence.need(not args.output.exists() and not args.output.is_symlink(), "Require fresh output")
    result = review(args.report, args.report_sha256)
    with args.output.open("x") as stream:
        json.dump(result, stream, sort_keys=True, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({k: result[k] for k in ("status", "genes", "changed_pairs", "localized_summary")}))
