import copy

import pytest

from benchmark_tools.prepare_blast_recovery_partition import partition


def example():
    blocks = [dict(query="a", input_ordinal_0based=0, start=0, end=10, rows=1,
                   sha256="a"*64, final_observed_query=False, reuse_authorized=False),
              dict(query="c", input_ordinal_0based=2, start=10, end=20, rows=1,
                   sha256="b"*64, final_observed_query=True, reuse_authorized=False)]
    boundary = dict(block_count=2, last_query="c", last_query_start=10, prefix_end=20)
    return ["a", "b", "c", "d"], blocks, boundary


def test_exhaustive_disjoint_partition_excludes_entire_final_query():
    rows = list(partition(*example()))
    assert [r["query"] for r in rows] == ["a", "b", "c", "d"]
    assert [r["disposition"] for r in rows] == ["candidate_prefix_not_admitted", "replay_absent", "replay_final_incomplete", "replay_absent"]
    assert rows[0]["retained_end"] == 10
    assert all(r["retained_start"] is r["retained_end"] is None for r in rows[1:])


@pytest.mark.parametrize("problem", ["duplicate_ids", "wrong_id", "wrong_order", "gap", "overlap", "empty_block", "hash", "authorized", "no_final", "early_final", "extra", "boundary", "count"])
def test_bad_partition_fails(problem):
    ids, blocks, boundary = example()
    if problem == "duplicate_ids":
        ids[-1] = "a"
    elif problem == "wrong_id":
        blocks[0]["query"] = "wrong"
    elif problem == "wrong_order":
        blocks[1]["input_ordinal_0based"] = 0
    elif problem in {"gap", "overlap"}:
        blocks[1]["start"] += 1 if problem == "gap" else -1
    elif problem == "empty_block":
        blocks[1]["rows"] = 0
    elif problem == "hash":
        blocks[0]["sha256"] = "bad"
    elif problem == "authorized":
        blocks[0]["reuse_authorized"] = True
    elif problem == "no_final":
        blocks[1]["final_observed_query"] = False
    elif problem == "early_final":
        blocks[0]["final_observed_query"] = True
    elif problem == "extra":
        blocks.append(dict(copy.deepcopy(blocks[1]), input_ordinal_0based=8))
    elif problem == "boundary":
        boundary["last_query_start"] = 11
    else:
        boundary["block_count"] = 3
    with pytest.raises(ValueError):
        list(partition(ids, blocks, boundary))
