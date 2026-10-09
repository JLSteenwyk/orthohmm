"""Check frozen selection independently of a preferred stage explanation."""

import copy
import hashlib
import json
from pathlib import Path

import pytest

from benchmark_tools import select_controlled_fragment_trace as current


def fixture():
    owners = dict(a="A", b="B", c="C", d="D", e="E")
    flags = dict(a=False, b=False, c=True, d=True, e=False)
    truth = {("a", "b"), ("a", "c"), ("c", "d"), ("a", "e")}
    baseline = {("a", "b"), ("a", "c"), ("b", "e"), ("a", "d")}
    fragment = {("a", "b"), ("c", "d"), ("b", "c"), ("b", "e")}
    return [dict(seed=s, owners=owners.copy(), flags=flags.copy(), truth=truth.copy(),
        predictions={"baseline": {m:baseline.copy() for m in current.METHODS},
                     "fragment": {m:fragment.copy() for m in current.METHODS}}) for s in current.SEEDS]


def test_exact_categories_and_overlaps():
    data = fixture()[0]
    output = current.categories(data["truth"], data["predictions"]["baseline"][current.METHODS[0]],
                                data["predictions"]["fragment"][current.METHODS[0]])
    assert output == dict(fragment_fn={("a", "c"), ("a", "e")}, fragment_fp={("b", "c"), ("b", "e")},
        new_fn={("a", "c")}, new_fp={("b", "c")}, recovered_fn={("c", "d")},
        removed_fp={("a", "d")}, retained_tp={("a", "b")})


def test_all_84_bins_hash_minimum_and_no_mutation():
    data = fixture()
    before = copy.deepcopy(data)
    bins, cases = current.select(data)
    assert data == before and len(bins) == 84
    assert len(cases) == 36 and sum(r["status"] == "empty_bin" for r in bins) == 48
    for row in bins:
        possible = []
        for d in data:
            category = current.categories(d["truth"], d["predictions"]["baseline"][row["method"]],
                                          d["predictions"]["fragment"][row["method"]])[row["category"]]
            for a, b in category:
                if sum(d["flags"][g] for g in (a,b)) == row["fragment_endpoints"]:
                    digest = hashlib.sha256(f"{row['method']}:{row['category']}:{d['seed']}:{a}:{b}".encode()).hexdigest()
                    possible.append((digest, d["seed"], a, b))
        assert row["eligible_records"] == len(possible)
        if possible:
            assert (row["selection_sha256"], row["seed"], row["gene_a"], row["gene_b"]) == min(possible)
        else:
            assert all(row[k] is None for k in ("seed", "gene_a", "gene_b", "selection_sha256"))
    assert current.select(list(reversed(data))) == (bins, cases)


def test_deterministic_collision_tie_break(monkeypatch):
    monkeypatch.setattr(current, "selection_key", lambda method, category, seed, pair: ("0"*64,seed,*pair))
    _, cases = current.select(fixture())
    assert all(row["seed"] == current.SEEDS[0] for row in cases)


@pytest.mark.parametrize("change", ["seed", "missing_method", "untyped_flag", "unknown_gene", "noncanonical", "same_species"])
def test_malformed_universe_refused(change):
    data = fixture()
    if change == "seed": data[0]["seed"] = data[1]["seed"]
    elif change == "missing_method": data[0]["predictions"]["baseline"].pop(current.METHODS[0])
    elif change == "untyped_flag": data[0]["flags"]["a"] = 0
    elif change == "unknown_gene": data[0]["truth"].add(("a","missing"))
    elif change == "noncanonical": data[0]["truth"].add(("c","a"))
    elif change == "same_species": data[0]["owners"]["b"] = "A"
    with pytest.raises(ValueError): current.select(data)


def test_hash_checked_original_reference(tmp_path):
    path = tmp_path / "data.json"
    path.write_text(json.dumps(dict(gene="a")))
    ref = current.record(path)
    evidence = {}
    assert current.checked(dict(absolute_path=ref["path"],bytes=ref["bytes"],sha256=ref["sha256"]), evidence) == path
    path.write_text("changed")
    with pytest.raises(ValueError): current.checked(ref, evidence)


def test_existing_selection_guard_precedes_loading(tmp_path):
    with pytest.raises(FileExistsError): current.run(tmp_path, tmp_path)
