"""Independent lookup oracle and synthetic full-readback refusals, not production admission."""

import csv
import json
from pathlib import Path

import numpy as np
import pytest

from benchmark_tools import readback_native_qfo_swiss_search_support as reader


def test_retained_actual_readback_receipt_without_repeating_large_scan():
    results = Path(reader.__file__).parent / "results"
    retained = json.loads((results / "native_qfo_swiss_search_support_code_readback_20261006.json").read_text())
    primary = json.loads(Path(retained["report"]["path"]).read_text())
    assert retained["source"] == reader.record(reader.__file__)
    assert retained["report"] == reader.record(retained["report"]["path"])
    assert retained["summary"] == primary["summary"]
    assert (retained["changed_pairs_checked"], retained["selected_genes_checked"], retained["checkpoint_rows_checked"],
            retained["selected_directed_hit_records_checked"], retained["different_selected_score_multisets"]) == (
                2023, 361, 181374654, 5884, 0)
    assert sum(row["before"] == "TP" for row in primary["cases"]) == 334
    assert sum(row["before"] == "FP" for row in primary["cases"]) == 1689
    assert all(row["selected_directed_score_multisets_identical"] is True for row in primary["cases"])
    for cell in reader.CELLS:
        hits = [hit for row in primary["cases"] for direction in ("gene_a_to_b", "gene_b_to_a")
                for hit in row["direct_search"][cell][direction]]
        assert len(hits) == len({h["row"] for h in hits}) == 2942
    for flag in ("new_scoring_or_admission", "uncertainty_admitted", "scientific_timings_admitted",
                 "independent_confirmation", "publication_ready"):
        assert retained[flag] is primary[flag] is False


@pytest.mark.parametrize("chunk", (1, 2, 3, 1000000))
def test_sorted_code_lookup_matches_exhaustive_oracle(chunk):
    q = np.array([0, 1, 2, 0], dtype=np.int32)
    t = np.array([1, 0, 1, 1], dtype=np.int32)
    s = np.array([10., 11., 5., 12.], dtype=np.float64)
    pairs = {(0, 1), (1, 0), (0, 2)}
    expected = {p: [] for p in pairs}
    for index in range(len(q)):
        pair = (int(q[index]), int(t[index]))
        if pair in expected:
            expected[pair].append(dict(row=index, score=float(s[index])))
    assert reader.scan_codes(q, t, s, pairs, 3, chunk) == expected


@pytest.mark.parametrize("fault", ("dtype", "length", "nan", "bounds", "query", "overflow", "chunk"))
def test_invalid_integer_scan_refuses(fault):
    q, t = np.array([0], dtype=np.int32), np.array([1], dtype=np.int32)
    s, pairs, genes, chunk = np.array([10.], dtype=np.float64), {(0, 1)}, 3, 1
    if fault == "dtype":
        q = q.astype(np.int64)
    elif fault == "length":
        t = t[:0]
    elif fault == "nan":
        s[0] = np.nan
    elif fault == "bounds":
        t[0] = 3
    elif fault == "query":
        pairs = {(0, 0)}
    elif fault == "overflow":
        genes = 3037000500
    elif fault == "chunk":
        chunk = True
    with pytest.raises(ValueError):
        reader.scan_codes(q, t, s, pairs, genes, chunk)


@pytest.mark.parametrize("fault", (None, "scope", "source", "case_inventory", "native_pin", "manifest_pin", "gene_id",
                                 "wrong_score", "missed_hit", "spurious_hit", "summary", "equality", "same_species", "ledger"))
def test_full_readback_on_synthetic_bound_arrays(tmp_path, fault):
    def dump(path, value):
        path.write_text(json.dumps(value))
        return reader.record(path)

    localization_ledger = tmp_path / "localized.tsv"
    localization_ledger.write_text("family\tprotein_a\tprotein_b\tbefore\tafter\tgene_a\tgene_b\n"
                                   "ref\ta\tb\tTP\tFN\ta\tb\n")
    localization = dump(tmp_path / "localization.json", dict(pair_ledger=reader.record(localization_ledger)))
    checkpoints = []
    for i, cell in enumerate(reader.CELLS):
        directory = tmp_path / str(i)
        directory.mkdir()
        (directory / "gene_names.txt").write_text("a\nb\nc\n")
        arrays = dict(gene_to_species=np.array([0, 0 if fault == "same_species" else 1, 2], dtype=np.int32),
            hit_queries=np.array([0, 1, 2, 0], dtype=np.int32), hit_targets=np.array([1, 0, 1, 1], dtype=np.int32),
            hit_scores=np.array([10., 11., 5., 12.], dtype=np.float64))
        for name, array in arrays.items():
            np.save(directory / (name + ".npy"), array, allow_pickle=False)
        refs = {name: reader.record(directory / name) for name in reader.FILES if name != "manifest.json"}
        manifest = dict(schema_version=1, complete=True, genes=3, hits=4,
                        files={name: {k: ref[k] for k in ("bytes", "sha256")} for name, ref in refs.items()})
        if fault == "manifest_pin":
            manifest["files"]["hit_queries.npy"]["sha256"] = "0" * 64
        refs["manifest.json"] = dump(directory / "manifest.json", manifest)
        validation_refs = list(refs.values())
        if fault == "native_pin":
            validation_refs.pop()
        validation = dump(directory / "validation.json", dict(native_outputs_validated=True, cell=cell,
                                                               checked_files=validation_refs))
        checkpoints.append(dict(cell=cell, files=refs, native_output_validation=validation,
                                previously_inventoried=True, genes=3, hits=4))
    forward = [dict(row=0, score=10.), dict(row=3, score=12.)]
    reverse = [dict(row=1, score=11.)]
    case = dict(family="ref", protein_a="a", protein_b="b", before="TP", after="FN", gene_a="a", gene_b="b",
        query_id_a=0, query_id_b=1, selected_directed_score_multisets_identical=True,
        direct_search={cell: dict(support="both_directions", gene_a_to_b=list(forward), gene_b_to_a=list(reverse))
                       for cell in reader.CELLS})
    summary = [dict(before=label, cell=cell, support=category,
                    pairs=int(label == "TP" and category == "both_directions"))
               for label in ("TP", "FP") for cell in reader.CELLS
               for category in ("no_direct_hit", "one_direction", "both_directions")]
    ledger = tmp_path / "support.tsv"
    ledger.write_text("family\tprotein_a\tprotein_b\tbefore\tafter\tr0_support\tr1_support\t"
        "r0_forward_records\tr0_reverse_records\tr1_forward_records\tr1_reverse_records\tselected_score_multisets_identical\n"
        "ref\ta\tb\tTP\tFN\tboth_directions\tboth_directions\t2\t1\t2\t1\tTrue\n")
    report = dict(schema="native_qfo_swiss_direct_search_support_v1",
        source=reader.record(Path(reader.__file__).with_name("trace_native_qfo_swiss_search_support.py")),
        localization=localization, support_ledger=reader.record(ledger), checkpoints=checkpoints,
        cases=[case], changed_pairs=1, selected_genes=2, summary=summary,
        selected_pairs_with_different_score_multisets=0, new_scoring_or_admission=False, uncertainty_admitted=False,
        scientific_timings_admitted=False, independent_confirmation=False, publication_ready=False)
    if fault == "scope":
        report["uncertainty_admitted"] = True
    elif fault == "source":
        report["source"] = reader.record(__file__)
    elif fault == "case_inventory":
        case["before"] = "FP"
    elif fault == "gene_id":
        case["query_id_a"] = 2
    elif fault == "wrong_score":
        case["direct_search"][reader.CELLS[0]]["gene_a_to_b"][0] = dict(row=0, score=99.)
    elif fault == "missed_hit":
        case["direct_search"][reader.CELLS[0]]["gene_a_to_b"].pop()
    elif fault == "spurious_hit":
        case["direct_search"][reader.CELLS[0]]["gene_a_to_b"].append(dict(row=2, score=5.))
    elif fault == "summary":
        summary[0]["pairs"] += 1
    elif fault == "equality":
        case["selected_directed_score_multisets_identical"] = False
    elif fault == "ledger":
        ledger.write_text(ledger.read_text().replace("\t2\t1\t2\t1\t", "\t1\t1\t2\t1\t"))
        report["support_ledger"] = reader.record(ledger)
    path = tmp_path / "report.json"
    report_ref = dump(path, report)
    if fault is None:
        result = reader.verify(path, report_ref["sha256"], chunk_size=2)
        assert result["checkpoint_rows_checked"] == 8 and result["selected_directed_hit_records_checked"] == 6
        assert result["different_selected_score_multisets"] == 0 and result["uncertainty_admitted"] is False
    else:
        with pytest.raises(ValueError):
            reader.verify(path, report_ref["sha256"], chunk_size=2)
