"""Independently enumerate selected significant hits by sorted integer-key lookup."""

import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.readback_native_qfo_swiss_transitions import record, require

CELLS = ("p0_c0_r0", "p0_c0_r1")
FILES = {"manifest.json", "gene_names.txt", "gene_to_species.npy", "hit_queries.npy", "hit_targets.npy", "hit_scores.npy"}


def scan_codes(q, t, s, pairs, genes, chunk_size=1000000):
    require(type(genes) is int and 1 < genes <= 3037000499 and type(chunk_size) is int and chunk_size > 0,
            "Invalid bounded integer-key controls")
    require(q.dtype == t.dtype == np.dtype("int32") and s.dtype == np.dtype("float64")
            and q.ndim == t.ndim == s.ndim == 1 and len(q) == len(t) == len(s), "Invalid hit array shape/dtype")
    require(pairs and all(type(a) is int and type(b) is int and 0 <= a < genes and 0 <= b < genes and a != b
                         for a, b in pairs), "Invalid directed query inventory")
    keys = np.asarray(sorted(a * genes + b for a, b in pairs), dtype=np.int64)
    result = {pair: [] for pair in pairs}
    for start in range(0, len(q), chunk_size):
        qb, tb, sb = (values[start:start + chunk_size] for values in (q, t, s))
        require(((qb >= 0) & (qb < genes)).all() and ((tb >= 0) & (tb < genes)).all()
                and np.isfinite(sb).all(), "Invalid hit endpoint or numeric score")
        codes = qb.astype(np.int64) * genes + tb
        positions = np.searchsorted(keys, codes)
        valid = positions < len(keys)
        matched = valid & (keys[np.minimum(positions, len(keys) - 1)] == codes)
        for offset in np.flatnonzero(matched):
            code = int(codes[offset])
            result[code // genes, code % genes].append(dict(row=start + int(offset), score=float(sb[offset])))
    return result


def verify(report_path, report_sha, chunk_size=1000000):
    report_ref = record(report_path)
    require(report_ref["sha256"] == report_sha, "Changed search-support report")
    report = json.loads(Path(report_path).read_text())
    require(report["schema"] == "native_qfo_swiss_direct_search_support_v1"
            and report["source"] == record(Path(__file__).with_name("trace_native_qfo_swiss_search_support.py"))
            and all(report[k] is False for k in ("new_scoring_or_admission", "uncertainty_admitted",
                "scientific_timings_admitted", "independent_confirmation", "publication_ready")),
            "Wrong direct-search diagnostic source/scope")
    checked = [report_ref, report["source"], report["localization"], report["support_ledger"]]
    for ref in checked:
        require(record(ref["path"]) == ref, "Changed direct search evidence")
    localized = json.loads(Path(report["localization"]["path"]).read_text())
    require(record(localized["pair_ledger"]["path"]) == localized["pair_ledger"], "Changed localized-pair inventory")
    checked.append(localized["pair_ledger"])
    with open(localized["pair_ledger"]["path"], newline="") as stream:
        original = list(csv.DictReader(stream, delimiter="\t"))
    cases = report["cases"]
    require(len(cases) == report["changed_pairs"] == len(original)
            and [{k: r[k] for k in ("family", "protein_a", "protein_b", "before", "after", "gene_a", "gene_b")}
                 for r in cases] == [{k: r[k] for k in ("family", "protein_a", "protein_b", "before", "after", "gene_a", "gene_b")}
                                   for r in original], "Changed or incomplete search-support cases")
    require([r["cell"] for r in report["checkpoints"]] == list(CELLS), "Changed checkpoint inventory")
    pairs = {(r["query_id_a"], r["query_id_b"]) for r in cases}
    require(len(pairs) == len(cases) and len({tuple(sorted(p)) for p in pairs}) == len(cases), "Duplicate pair case")
    directed = pairs | {(b, a) for a, b in pairs}
    targets = {r[k] for r in cases for k in ("gene_a", "gene_b")}
    require(len(targets) == report["selected_genes"], "Changed target gene count")
    scans = []
    for view in report["checkpoints"]:
        require(view["previously_inventoried"] is True and set(view["files"]) == FILES, "Uninventoried checkpoint")
        refs = [view["native_output_validation"], *view["files"].values()]
        for ref in refs:
            require(record(ref["path"]) == ref, "Changed native checkpoint or validation evidence")
        checked.extend(refs)
        validation = json.loads(Path(view["native_output_validation"]["path"]).read_text())
        require(validation["native_outputs_validated"] is True and validation["cell"] == view["cell"]
                and all(ref in validation["checked_files"] for ref in view["files"].values()),
                "Checkpoint not in original native validation")
        manifest = json.loads(Path(view["files"]["manifest.json"]["path"]).read_text())
        require(manifest["complete"] is True and manifest["genes"] == view["genes"]
                and manifest["hits"] == view["hits"] and all(manifest["files"][name] ==
                    {k: ref[k] for k in ("bytes", "sha256")} for name, ref in view["files"].items() if name != "manifest.json"),
                "Checkpoint manifest differs")
        names = {}
        count = 0
        with open(view["files"]["gene_names.txt"]["path"]) as stream:
            for index, line in enumerate(stream):
                name = line.rstrip("\r\n")
                if name in targets:
                    require(name not in names, "Duplicate target gene")
                    names[name] = index
                count += 1
        require(count == view["genes"] and set(names) == targets and all(
            names[r["gene_a"]] == r["query_id_a"] and names[r["gene_b"]] == r["query_id_b"] for r in cases),
            "Wrong target gene indexing")
        ownership = np.load(view["files"]["gene_to_species.npy"]["path"], mmap_mode="r", allow_pickle=False)
        require(ownership.shape == (view["genes"],) and ownership.dtype == np.dtype("int32")
                and all(ownership[a] != ownership[b] for a, b in pairs), "Invalid cross-species query")
        arrays = [np.load(view["files"][name]["path"], mmap_mode="r", allow_pickle=False)
                  for name in ("hit_queries.npy", "hit_targets.npy", "hit_scores.npy")]
        require(all(array.shape == (view["hits"],) and not array.flags.writeable for array in arrays), "Wrong hit length/mapping")
        scans.append(scan_codes(*arrays, directed, view["genes"], chunk_size))
        del arrays, ownership
    summary, different, positive_records = Counter(), 0, 0
    for row in cases:
        a, b = row["query_id_a"], row["query_id_b"]
        for cell, scan in zip(CELLS, scans):
            forward, reverse = scan[a, b], scan[b, a]
            category = ("both_directions" if len(forward) and len(reverse) else
                        "one_direction" if len(forward) + len(reverse) else "no_direct_hit")
            require(row["direct_search"][cell] == dict(support=category, gene_a_to_b=forward, gene_b_to_a=reverse),
                    "Independent complete hit enumeration differs")
            summary[row["before"], cell, category] += 1
            positive_records += len(forward) + len(reverse)
        equal = all(sorted(h["score"] for h in scans[0][p]) == sorted(h["score"] for h in scans[1][p])
                    for p in ((a, b), (b, a)))
        require(row["selected_directed_score_multisets_identical"] is equal, "Changed selected-score equality flag")
        different += not equal
    expected = [dict(before=label, cell=cell, support=category, pairs=summary[label, cell, category])
        for label in ("TP", "FP") for cell in CELLS
        for category in ("no_direct_hit", "one_direction", "both_directions")]
    require(expected == report["summary"] and different == report["selected_pairs_with_different_score_multisets"],
            "Changed direct-hit summary")
    with open(report["support_ledger"]["path"], newline="") as stream:
        ledger = list(csv.DictReader(stream, delimiter="\t"))
    require(len(ledger) == len(cases), "Incomplete support ledger")
    for row, actual in zip(cases, ledger):
        left, right = (row["direct_search"][cell] for cell in CELLS)
        wanted = {k: str(row[k]) for k in ("family", "protein_a", "protein_b", "before", "after")}
        wanted.update(r0_support=left["support"], r1_support=right["support"],
            r0_forward_records=str(len(left["gene_a_to_b"])), r0_reverse_records=str(len(left["gene_b_to_a"])),
            r1_forward_records=str(len(right["gene_a_to_b"])), r1_reverse_records=str(len(right["gene_b_to_a"])),
            selected_score_multisets_identical=str(row["selected_directed_score_multisets_identical"]))
        require(actual == wanted, "Support ledger differs")
    for ref in checked:
        require(record(ref["path"]) == ref, "Evidence changed during independent scan")
    return dict(schema="native_qfo_swiss_search_support_code_readback_v1", report=report_ref, source=record(__file__),
        checked_inputs=checked, chunk_size=chunk_size, changed_pairs_checked=len(cases), selected_genes_checked=len(targets),
        checkpoint_rows_checked=sum(view["hits"] for view in report["checkpoints"]),
        selected_directed_hit_records_checked=positive_records, summary=expected,
        different_selected_score_multisets=different, new_scoring_or_admission=False, uncertainty_admitted=False,
        scientific_timings_admitted=False, independent_confirmation=False, publication_ready=False,
        limitations=["Independent sorted integer-key enumeration, including complete absent directions; same retained data.",
                     "No inference, scoring, orthology confidence, biological replication or new transitive admission.",
                     "Selected-hit equivalence is not whole-search equivalence or proof of graph/orthology support."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--report-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(not args.output.exists() and not args.output.is_symlink(), "Output already exists")
    result = verify(args.report, args.report_sha256)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps({k: result[k] for k in ("changed_pairs_checked", "selected_genes_checked", "checkpoint_rows_checked",
        "selected_directed_hit_records_checked", "different_selected_score_multisets")}))
