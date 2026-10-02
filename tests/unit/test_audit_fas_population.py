import ast
import hashlib
import json
from pathlib import Path
import sqlite3

import numpy as np
import pytest

from benchmark_tools.audit_fas_population import QUERY, database_pairs, encode, lookup_index, match, summarize_population

SOURCE_DIRECTORY = Path(__file__).resolve().parents[1] / "samples/native_qfo_fas"
SOURCE_PINS = {
    "fas_benchmark.py": (9930, "1045c57f4d0f4787bec3d1f0690799df63c68dcc67ccfd5a925c338e3d33661d"),
    "LICENSE": (17510, "1658efc25eda746f9a30c3cbe8728134eb390034092d848e996080d6f52a1f6f"),
}


def checked_native_source(directory):
    for name, expected in SOURCE_PINS.items():
        path = directory / name
        if path.is_symlink():
            raise ValueError("Native source reference must not be a symlink")
        data = path.read_bytes()
        if (len(data), hashlib.sha256(data).hexdigest()) != expected:
            raise ValueError("Changed native source reference")
    return directory / "fas_benchmark.py"


@pytest.fixture
def retained_native_fas_source():
    return checked_native_source(SOURCE_DIRECTORY)


@pytest.mark.parametrize("name", list(SOURCE_PINS))
@pytest.mark.parametrize("change", ["bytes", "missing", "symlink"])
def test_changed_native_reference_is_rejected(tmp_path, name, change):
    for filename in SOURCE_PINS:
        (tmp_path / filename).write_bytes((SOURCE_DIRECTORY / filename).read_bytes())
    path = tmp_path / name
    if change == "bytes":
        path.write_bytes(path.read_bytes() + b"\n")
    else:
        path.unlink()
        if change == "symlink":
            path.symlink_to(SOURCE_DIRECTORY / name)
    with pytest.raises((ValueError, FileNotFoundError)):
        checked_native_source(tmp_path)


def test_native_source_provenance_is_explicit(retained_native_fas_source):
    manifest = json.loads((SOURCE_DIRECTORY / "provenance.json").read_text())
    assert manifest["upstream_repository"] == "https://github.com/qfo/benchmark-webservice"
    assert manifest["upstream_commit"] == "c0854a96c1a0fd7f2a891d971af0863002fabc90"
    assert {r["path"]: (r["bytes"], r["sha256"]) for r in manifest["files"]} == SOURCE_PINS
    assert manifest["modifications_to_upstream_files"] is False
    assert manifest["biological_data_included"] is False
    assert manifest["darwin_binary_included"] is False
    assert manifest["full_native_benchmark_execution"] is False


@pytest.mark.parametrize("a,b", [(0, 0), (0, 1), (2**32 - 1, 0), (4, 8)])
def test_canonical_encoding(a, b):
    assert encode(a, b) == encode(b, a) == (min(a, b) << 32) | max(a, b)


@pytest.mark.parametrize("a,b", [(-1, 0), (0, -1), (2**32, 0), (0, 2**32)])
def test_encoding_range(a, b):
    with pytest.raises(ValueError):
        encode(a, b)


def test_index_matches_actual_native_loader_overwrites_and_invalid_values(tmp_path, retained_native_fas_source):
    entries = [("A_B", [.2, .4]), ("B_A", [.6, .8]),
               ("A_C", [.5, .5]), ("C_A", ["NA", "NA"]),
               ("B_C", ["NA", "NA"])]
    ids = {}
    codes, scores, summary = lookup_index(entries, ids)
    source = retained_native_fas_source
    function = next(n for n in ast.parse(source.read_text()).body
                    if isinstance(n, ast.FunctionDef) and n.name == "load_precomputed_fas_scores")
    namespace = dict(numpy=np, load_json_file=lambda _: dict(entries))
    import logging
    namespace["logger"] = logging.getLogger("test")
    namespace["Path"] = Path
    exec(compile(ast.Module(body=[function], type_ignores=[]), str(source), "exec"), namespace)
    native = namespace["load_precomputed_fas_scores"](tmp_path / "unused.json")
    assert summary == dict(entries_scanned=5, invalid_value_entries=2,
                           valid_canonical_pairs=2, canonical_overwrites=1)
    for pair, value in native.items():
        present, observed = match(codes, scores, np.array([encode(*(ids[a] for a in pair))], dtype=np.uint64))
        assert present.tolist() == [True]
        assert observed.tolist() == pytest.approx([value])


@pytest.mark.parametrize("directions", [[float("nan"), 0], [0, float("inf")],
                                       [-.01, .5], [1.01, .5], [], [.5], [.2, .3, .4], "0.5"])
def test_unsupported_lookup_values_fail(directions):
    with pytest.raises(ValueError):
        lookup_index([("A_B", directions)], {})


@pytest.mark.parametrize("second", [[.5, .5], ["NA", "NA"]])
def test_exact_duplicate_keys_not_silently_streamed(second):
    with pytest.raises(ValueError, match="Duplicate JSON"):
        lookup_index([("A_B", [.3, .4]), ("A_B", second)], {})


def test_empty_lookup_and_new_identifiers():
    ids = {}
    codes, scores, summary = lookup_index([], ids)
    assert summary["valid_canonical_pairs"] == 0
    result = summarize_population([("A", "B")], ids, {"A", "B"}, codes, scores, 1)
    assert result["missing"] == 1 and result["precomputed_mean"] is None
    assert result["hypothetical_full_mean_bounds"] == [0., 1.]


@pytest.mark.parametrize("batch", [1, 2, 100])
def test_population_rules_match_native_precedence_and_bound(batch):
    ids = {}
    codes, scores, _ = lookup_index([("A_B", [.2, .4]), ("A_C", [.8, 1])], ids)
    # A/C remains precomputed even though C lacks annotations. Alias pairs skip.
    pairs = [("A", "B"), ("A", "C"), ("A", "D"), ("B", "X"), ("A", "G_ALIAS")]
    result = summarize_population(pairs, ids, {"A", "B", "D"}, codes, scores, batch)
    assert {k: result[k] for k in ("distinct_query_pairs", "skipped_alias_pairs", "precomputed", "missing", "unannotated", "eligible_pairs")} == dict(
        distinct_query_pairs=5, skipped_alias_pairs=1, precomputed=2, missing=1, unannotated=1, eligible_pairs=3)
    assert result["precomputed_score_sum"] == pytest.approx(1.2)
    assert result["hypothetical_full_mean_bounds"] == pytest.approx([.4, 2.2 / 3])


def test_no_eligible_population_has_no_invented_mean():
    codes, scores, _ = lookup_index([], {})
    with pytest.raises(ValueError, match="No eligible"):
        summarize_population([("A", "B")], {}, set(), codes, scores)


def test_query_literal_matches_retained_native_source(retained_native_fas_source):
    source = retained_native_fas_source
    node = next(n for n in ast.walk(ast.parse(source.read_text()))
                if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == "query" for t in n.targets)
                and isinstance(n.value, ast.Constant) and str(n.value.value).startswith("SELECT DISTINCT"))
    assert node.value.value == QUERY


def test_database_distinct_aliases_and_read_only(tmp_path):
    path = tmp_path / "predictions.db"
    with sqlite3.connect(path) as connection:
        connection.execute("CREATE TABLE proteomes(prot_nr INTEGER, uniprot_id TEXT)")
        connection.execute("CREATE TABLE orthologs(prot_nr1 INTEGER, prot_nr2 INTEGER)")
        connection.executemany("INSERT INTO proteomes VALUES (?, ?)", [(1, "A"), (2, "B"), (3, "A"), (4, "C_ALIAS")])
        connection.executemany("INSERT INTO orthologs VALUES (?, ?)", [(1, 2), (2, 1), (3, 2), (1, 2), (1, 4)])
    before = path.read_bytes()
    assert sorted(database_pairs(path)) == [("A", "B"), ("A", "C_ALIAS")]
    assert path.read_bytes() == before
    Path(str(path) + "-wal").touch()
    with pytest.raises(ValueError, match="sidecar"):
        list(database_pairs(path))
