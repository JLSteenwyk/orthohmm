"""Bounded direct-search-hit trace for every localized native SwissTrees exclusion."""

import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.export_native_factorial_progress import load, require
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

CELLS = ("p0_c0_r0", "p0_c0_r1")
FILES = ("gene_names.txt", "gene_to_species.npy", "hit_queries.npy", "hit_targets.npy", "hit_scores.npy")


def target_ids(path, targets, expected_genes):
    selected, previous, count = {}, None, 0
    with path.open() as stream:
        for line in stream:
            name = line.rstrip("\r\n")
            require(name and (previous is None or previous < name), "Invalid lexical gene-name universe")
            if name in targets:
                selected[name] = count
            previous, count = name, count + 1
    require(count == expected_genes and set(selected) == targets, "Missing endpoint or wrong gene universe")
    return selected


def scan_hits(queries, targets, scores, directed, genes, chunk_size=1000000):
    require(type(chunk_size) is int and chunk_size > 0 and type(genes) is int and genes > 0,
            "Invalid scan controls")
    require(queries.dtype == targets.dtype == np.dtype("int32") and scores.dtype == np.dtype("float64")
            and queries.ndim == targets.ndim == scores.ndim == 1
            and len(queries) == len(targets) == len(scores), "Changed checkpoint array shape/dtype")
    require(directed and all(type(a) is int and type(b) is int and 0 <= a < genes and 0 <= b < genes
            and a != b for a, b in directed), "Invalid directed-pair inventory")
    selected = np.zeros(genes, dtype=bool)
    for a, b in directed:
        selected[a] = selected[b] = True
    result = {pair: [] for pair in directed}
    for start in range(0, len(queries), chunk_size):
        q, t, s = (array[start:start + chunk_size] for array in (queries, targets, scores))
        require(((q >= 0) & (q < genes)).all() and ((t >= 0) & (t < genes)).all()
                and np.isfinite(s).all(), "Invalid checkpoint endpoint or score")
        mask = selected[q] & selected[t]
        for offset in np.flatnonzero(mask):
            pair = (int(q[offset]), int(t[offset]))
            if pair in result:
                result[pair].append(dict(row=start + int(offset), score=float(s[offset])))
    return result


def support(forward, reverse):
    return "both_directions" if forward and reverse else "one_direction" if forward or reverse else "no_direct_hit"


def binding(admission, index, evidence):
    stage = admission["conversion"]
    require(stage["cell"] == CELLS[index] and stage["native_index"] == index + 6, "Changed native identity")
    ref = stage["terminal_review"] if index == 0 else stage["scientific_recovery"]
    review, _ = load(ref["path"], ref["sha256"], evidence)
    output_ref = review["reviews"]["outputs_or_failure"] if index == 0 else review["outputs"]
    output, _ = load(output_ref["path"], output_ref["sha256"], evidence)
    require(output["native_outputs_validated"] is True and output["cell"] == CELLS[index]
            and stage["native_input"] in output["checked_files"], "Require admitted native output inventory")
    if index:
        require(admission["resources"] is None and admission["scientific_timings_admitted"] is False
                and admission["eligible_for_timing_comparison"] is False, "Recovered science relabels failed timing")
    native = Path(stage["native_input"]["path"])
    working = native.parent if index == 0 else native.parent.parent / "orthohmm_working_res"
    directory = working / "high_sensitivity_checkpoint"
    pins = {r["path"]: r for r in output["checked_files"]}
    refs = {}
    for name in ("manifest.json", *FILES):
        path = str(directory / name)
        require(path in pins, "Checkpoint file absent from original native validation: " + name)
        ref = pins[path]
        check(ref)
        evidence.append(ref)
        refs[name] = ref
    manifest = json.loads(Path(refs["manifest.json"]["path"]).read_text())
    require(manifest["schema_version"] == 1 and manifest["complete"] is True
            and type(manifest["genes"]) is int and manifest["genes"] == output["input_genes"]
            and type(manifest["hits"]) is int and manifest["hits"] >= 0 and set(manifest["files"]) == set(FILES),
            "Invalid checkpoint completion/universe/inventory")
    require(all(manifest["files"][name] == {k: refs[name][k] for k in ("bytes", "sha256")} for name in FILES),
            "Checkpoint manifest differs from independently admitted file pins")
    evidence.extend([review["source"], output["source"]])
    return dict(cell=CELLS[index], native_output_validation=output_ref, files=refs,
        genes=manifest["genes"], hits=manifest["hits"], previously_inventoried=True)


def run(localization_path, localization_sha, output, ledger, chunk_size=1000000):
    require(output != ledger and all(not p.exists() and not p.is_symlink() for p in (output, ledger)),
            "Require distinct fresh output paths")
    require(type(chunk_size) is int and 0 < chunk_size <= 10000000, "Invalid bounded chunk size")
    evidence = []
    localized, localized_ref = load(localization_path, localization_sha, evidence)
    require(localized["schema"] == "native_qfo_swiss_reconciliation_localization_v1"
            and localized["source"] == record(Path(__file__).with_name("trace_native_qfo_swiss_reconciliation.py"))
            and all(localized[k] is False for k in ("new_scoring_or_admission", "uncertainty_admitted",
                "scientific_timings_admitted", "independent_confirmation", "publication_ready")),
            "Wrong localized source/scope")
    transition, transition_ref = load(localized["transition"]["path"], localized["transition"]["sha256"], evidence)
    require(transition_ref == localized["transition"] and transition["schema"] == "native_qfo_swiss_pair_transitions_v1"
            and [r["cell"] for r in transition["cells"]] == list(CELLS), "Changed transition binding")
    for ref in (localized["source"], localized["pair_ledger"], transition["changed_relations_ledger"]):
        check(ref)
        evidence.append(ref)
    with open(localized["pair_ledger"]["path"], newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    with open(transition["changed_relations_ledger"]["path"], newline="") as stream:
        original = list(csv.DictReader(stream, delimiter="\t"))
    require([{k: row[k] for k in ("family", "protein_a", "protein_b", "before", "after")} for row in rows] == original
            and len(rows) == localized["changed_pairs_traced"] == transition["comparison"]["changed_relations"],
            "Incomplete or altered changed-pair inventory")
    require(rows and all(row["before"] in ("TP", "FP") and row["after"] in ("FN", "TN") for row in rows),
            "Require native excluded positive predictions")
    checkpoints = []
    for index, row in enumerate(transition["cells"]):
        admission, observed = load(row["admission"]["path"], row["admission"]["sha256"], evidence)
        require(observed == row["admission"] and admission["accuracy_admitted"] is True,
                "Changed scientific admission")
        checkpoints.append(binding(admission, index, evidence))
    require(all({k: checkpoints[0]["files"][name][k] for k in ("sha256", "bytes")} ==
                {k: checkpoints[1]["files"][name][k] for k in ("sha256", "bytes")}
                for name in ("gene_names.txt", "gene_to_species.npy")), "Different checkpoint gene/ownership bytes")
    targets = {row[key] for row in rows for key in ("gene_a", "gene_b")}
    scans, pair_ids = [], None
    for view in checkpoints:
        refs = view["files"]
        ids = target_ids(Path(refs["gene_names.txt"]["path"]), targets, view["genes"])
        pairs = [(ids[row["gene_a"]], ids[row["gene_b"]]) for row in rows]
        require(len(set(tuple(sorted(p)) for p in pairs)) == len(pairs), "Duplicate changed native pair")
        require(pair_ids is None or pairs == pair_ids, "Changed target gene indexing")
        pair_ids = pairs
        owners = np.load(refs["gene_to_species.npy"]["path"], mmap_mode="r", allow_pickle=False)
        require(owners.shape == (view["genes"],) and owners.dtype == np.dtype("int32")
                and all(owners[a] != owners[b] for a, b in pairs), "Wrong ownership or same-species pair")
        arrays = [np.load(refs[name]["path"], mmap_mode="r", allow_pickle=False)
                  for name in ("hit_queries.npy", "hit_targets.npy", "hit_scores.npy")]
        require(all(array.shape == (view["hits"],) and not array.flags.writeable for array in arrays),
                "Wrong hit length or writable checkpoint mapping")
        directed = {p for a, b in pairs for p in ((a, b), (b, a))}
        scans.append(scan_hits(*arrays, directed, view["genes"], chunk_size))
        del arrays, owners
    cases, summary = [], Counter()
    for row, (a, b) in zip(rows, pair_ids):
        views = {}
        for cell, scan in zip(CELLS, scans):
            forward, reverse = scan[a, b], scan[b, a]
            category = support(forward, reverse)
            summary[(row["before"], cell, category)] += 1
            views[cell] = dict(support=category, gene_a_to_b=forward, gene_b_to_a=reverse)
        identical = all(sorted(hit["score"] for hit in views[CELLS[0]][direction]) ==
                        sorted(hit["score"] for hit in views[CELLS[1]][direction])
                        for direction in ("gene_a_to_b", "gene_b_to_a"))
        cases.append(dict(family=row["family"], protein_a=row["protein_a"], protein_b=row["protein_b"],
            before=row["before"], after=row["after"], gene_a=row["gene_a"], gene_b=row["gene_b"],
            query_id_a=a, query_id_b=b, same_root_hog=row["same_root_hog"] == "True",
            direct_search=views, selected_directed_score_multisets_identical=identical))
    for ref in evidence:
        check(ref)
    with ledger.open("x", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(("family", "protein_a", "protein_b", "before", "after", "r0_support", "r1_support",
                         "r0_forward_records", "r0_reverse_records", "r1_forward_records", "r1_reverse_records",
                         "selected_score_multisets_identical"))
        for row in cases:
            left, right = (row["direct_search"][cell] for cell in CELLS)
            writer.writerow((row["family"], row["protein_a"], row["protein_b"], row["before"], row["after"],
                left["support"], right["support"], len(left["gene_a_to_b"]), len(left["gene_b_to_a"]),
                len(right["gene_a_to_b"]), len(right["gene_b_to_a"]), row["selected_directed_score_multisets_identical"]))
    result = dict(schema="native_qfo_swiss_direct_search_support_v1", localization=localized_ref,
        source=record(__file__), evidence=evidence, checkpoints=checkpoints, chunk_size=chunk_size,
        changed_pairs=len(cases), selected_genes=len(targets), cases=cases, support_ledger=record(ledger),
        summary=[dict(before=label, cell=cell, support=category, pairs=summary[label, cell, category])
            for label in ("TP", "FP") for cell in CELLS
            for category in ("no_direct_hit", "one_direction", "both_directions")],
        selected_pairs_with_different_score_multisets=sum(not r["selected_directed_score_multisets_identical"] for r in cases),
        new_scoring_or_admission=False, uncertainty_admitted=False, scientific_timings_admitted=False,
        independent_confirmation=False, publication_ready=False,
        limitations=[
            "Retrospective complete changed-pair trace, not all reference/submitted predictions or selected defaults.",
            "Checkpoint stores significant hits before RB-NH graph/group stages, not inferred orthology or raw prefilter candidates.",
            "No hit does not distinguish prefilter exclusion, insignificant scoring or indirect graph connectivity.",
            "Stored search scores are not calibrated orthology confidence; directed duplicates and exact row offsets are retained.",
            "Selected directional score multisets ignore row order, not scores; agreement is not whole-hit-set equality.",
            "Original checkpoint files were inventoried; direct digest/array checks are not repeated transitive scientific admission.",
            "Initial HMM search remains on; this is not a non-HMM sensitivity-matched comparator or new accuracy endpoint.",
            "Failed R1 inference timing remains ineligible; diagnostic shared-host contention effects unknown/tool-dependent."])
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("localization", "output", "ledger"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--localization-sha256", required=True)
    parser.add_argument("--chunk-size", type=int, default=1000000)
    args = parser.parse_args()
    result = run(args.localization, args.localization_sha256, args.output, args.ledger, args.chunk_size)
    print(json.dumps({k: result[k] for k in ("changed_pairs", "selected_genes", "summary",
        "selected_pairs_with_different_score_multisets")}))
