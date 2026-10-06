"""Bounded search lookup agrees with an exhaustive independent small-array oracle."""

import json
from pathlib import Path

import numpy as np
import pytest

from benchmark_tools import trace_native_qfo_swiss_search_support as trace


@pytest.fixture
def hits():
    return (np.array([0, 1, 3, 0, 0, 2, 4], dtype=np.int32),
            np.array([1, 0, 0, 1, 2, 0, 4], dtype=np.int32),
            np.array([10., 11., 12., 13., 14., 15., 16.], dtype=np.float64))


@pytest.mark.parametrize("chunk", (1, 2, 3, 4, 1000000))
def test_bounded_scan_matches_exhaustive_oracle(hits, chunk):
    pairs = {(0, 1), (1, 0), (0, 2), (2, 0), (1, 2)}
    expected = {p: [] for p in pairs}
    for i, (a, b, score) in enumerate(zip(*hits)):
        if (int(a), int(b)) in expected:
            expected[int(a), int(b)].append(dict(row=i, score=float(score)))
    assert trace.scan_hits(*hits, pairs, 5, chunk) == expected
    assert len(expected[0, 1]) == 2 and expected[1, 2] == []


@pytest.mark.parametrize("fault", ("qtype", "ttype", "stype", "length", "dimension", "negative", "outside", "nan", "inf"))
def test_bad_arrays_refuse(hits, fault):
    q, t, s = (a.copy() for a in hits)
    if fault == "qtype":
        q = q.astype(np.int64)
    elif fault == "ttype":
        t = t.astype(np.int64)
    elif fault == "stype":
        s = s.astype(np.float32)
    elif fault == "length":
        t = t[:-1]
    elif fault == "dimension":
        q = q.reshape(1, -1)
    elif fault == "negative":
        q[2] = -1
    elif fault == "outside":
        t[2] = 5
    elif fault == "nan":
        s[2] = np.nan
    elif fault == "inf":
        s[2] = np.inf
    with pytest.raises(ValueError):
        trace.scan_hits(q, t, s, {(0, 1)}, 5, 2)


@pytest.mark.parametrize("chunk", (0, -1, True, 1.5))
def test_bad_chunk_refuses(hits, chunk):
    with pytest.raises(ValueError, match="controls"):
        trace.scan_hits(*hits, {(0, 1)}, 5, chunk)


@pytest.mark.parametrize("pairs", (set(), {(0, 0)}, {(-1, 2)}, {(0, 5)}, {(True, 2)}))
def test_bad_directed_inventory_refuses(hits, pairs):
    with pytest.raises(ValueError, match="inventory"):
        trace.scan_hits(*hits, pairs, 5)


def test_empty_hit_checkpoint_retains_absence():
    result = trace.scan_hits(np.array([], dtype=np.int32), np.array([], dtype=np.int32),
                             np.array([], dtype=np.float64), {(0, 1)}, 2)
    assert result == {(0, 1): []}


@pytest.mark.parametrize("left,right,expected", (([], [], "no_direct_hit"), ([1], [], "one_direction"),
                                               ([], [1], "one_direction"), ([1], [2], "both_directions")))
def test_directional_support(left, right, expected):
    assert trace.support(left, right) == expected


def test_lexical_target_index(tmp_path):
    path = tmp_path / "genes.txt"
    path.write_text("a\nb\nc\n")
    assert trace.target_ids(path, {"a", "c"}, 3) == dict(a=0, c=2)
    with pytest.raises(ValueError, match="gene universe"):
        trace.target_ids(path, {"a", "missing"}, 3)
    with pytest.raises(ValueError, match="gene universe"):
        trace.target_ids(path, {"a"}, 4)
    path.write_text("b\na\n")
    with pytest.raises(ValueError, match="lexical"):
        trace.target_ids(path, {"a"}, 2)


@pytest.mark.parametrize("fault", (None, "missing_pin", "manifest_pin", "incomplete", "identity", "failed_timing"))
def test_checkpoint_binding_with_synthetic_original_review(tmp_path, fault):
    directory = tmp_path / "native" / "orthohmm_phylogeny"
    directory.mkdir(parents=True)
    native = directory / "orthohmm_pairwise_orthologs.tsv"
    native.write_text("native fixture\n")
    checkpoint = directory.parent / "orthohmm_working_res" / "high_sensitivity_checkpoint"
    checkpoint.mkdir(parents=True)
    (checkpoint / "gene_names.txt").write_text("a\nb\nc\n")
    arrays = dict(gene_to_species=np.array([0, 1, 2], dtype=np.int32),
                  hit_queries=np.array([0], dtype=np.int32), hit_targets=np.array([1], dtype=np.int32),
                  hit_scores=np.array([10.], dtype=np.float64))
    for name, array in arrays.items():
        np.save(checkpoint / (name + ".npy"), array, allow_pickle=False)
    pins = [trace.record(checkpoint / name) for name in trace.FILES]
    manifest = dict(schema_version=1, complete=fault != "incomplete", genes=3, hits=1,
                    files={Path(ref["path"]).name: {k: ref[k] for k in ("bytes", "sha256")} for ref in pins})
    if fault == "manifest_pin":
        manifest["files"]["hit_queries.npy"]["sha256"] = "0" * 64
    (checkpoint / "manifest.json").write_text(json.dumps(manifest))
    refs = [*pins, trace.record(checkpoint / "manifest.json"), trace.record(native)]
    if fault == "missing_pin":
        refs.pop(0)
    output = tmp_path / "outputs.json"
    output.write_text(json.dumps(dict(native_outputs_validated=True, cell="p0_c0_r1", input_genes=3,
                                    checked_files=refs, source=trace.record(__file__))))
    review = tmp_path / "review.json"
    review.write_text(json.dumps(dict(outputs=trace.record(output), source=trace.record(__file__))))
    admission = dict(resources=None, scientific_timings_admitted=False, eligible_for_timing_comparison=False,
        conversion=dict(cell="wrong" if fault == "identity" else "p0_c0_r1", native_index=7,
                        scientific_recovery=trace.record(review), native_input=trace.record(native)))
    if fault == "failed_timing":
        admission["scientific_timings_admitted"] = True
    if fault is None:
        result = trace.binding(admission, 1, [])
        assert result["previously_inventoried"] is True and result["hits"] == 1
    else:
        with pytest.raises(ValueError):
            trace.binding(admission, 1, [])
