import hashlib
import io

import pytest

from benchmark_tools.merge_blast_recovery_blocks import ordered_blocks, copy_block, merge


def row(query):
    return (query+"\tsubject\t100\t10\t0\t0\t1\t10\t1\t10\t1e-8\t50\n").encode()


def block(query, start=0, path="old"):
    data = row(query)
    return dict(query=query, start=start, end=start+len(data), rows=1,
                sha256=hashlib.sha256(data).hexdigest(), path=path, final_observed_query=False)


def test_interleaved_missing_and_boundary_query():
    prefix = [block("a"), block("c", len(row("a")))]
    replay = [block("b", path="new"), block("d", len(row("b")), "new")]
    planned = list(ordered_blocks(["a", "b", "c", "nohit", "d"], prefix, replay,
                                  ["b", "nohit", "d"], prefix[-1]["end"]))
    assert [b["query"] for b in planned] == ["a", "b", "c", "d"]
    assert [b["source"] for b in planned] == ["prefix", "replay", "prefix", "replay"]


@pytest.mark.parametrize("problem", ["duplicate_input", "missing_prefix", "bad_end", "final_block",
    "reordered_replay", "extra_replay", "overlap", "duplicate_replay"])
def test_partition_errors(problem):
    queries, prefix, replay, ids, end = ["a", "b"], [block("a")], [block("b", path="new")], ["b"], len(row("a"))
    if problem == "duplicate_input": queries.append("a")
    elif problem == "missing_prefix": prefix = []
    elif problem == "bad_end": end -= 1
    elif problem == "final_block": prefix[0]["final_observed_query"] = True
    elif problem == "reordered_replay": ids = ["b", "a"]
    elif problem == "extra_replay": replay.append(block("z", path="new"))
    elif problem == "overlap": ids = ["a", "b"]
    else: ids.append("b")
    with pytest.raises(ValueError):
        list(ordered_blocks(queries, prefix, replay, ids, end))


def test_exact_interleaved_merge(tmp_path):
    old, new = tmp_path/"old", tmp_path/"new"
    old.write_bytes(row("a")+row("c")+b"excluded broken tail")
    new.write_bytes(row("b")+row("d"))
    prefix = [block("a"), block("c", len(row("a")))]
    replay = [block("b", path="new"), block("d", len(row("b")), "new")]
    report = merge(ordered_blocks(["a", "b", "c", "d"], prefix, replay, ["b", "d"], prefix[-1]["end"]),
                   {"old": old, "new": new}, tmp_path/"merged")
    data = b"".join(row(q) for q in "abcd")
    assert (tmp_path/"merged/all.blast.candidate").read_bytes() == data
    assert report["sha256"] == hashlib.sha256(data).hexdigest()
    assert report["rows"] == 4 and report["search_admitted"] is False
    assert old.read_bytes().endswith(b"excluded broken tail")
    with pytest.raises(FileExistsError):
        merge([], {}, tmp_path/"merged")


@pytest.mark.parametrize("problem", ["hash", "rows", "query", "truncated", "nul", "columns"])
def test_block_corruption_rejected(problem):
    data, entry = row("a"), block("a")
    if problem == "hash": entry["sha256"] = "0"*64
    elif problem == "rows": entry["rows"] = 2
    elif problem == "query": entry["query"] = "b"
    elif problem == "truncated": data = data[:-1]
    elif problem == "nul": data = data.replace(b"subject", b"sub\0ect")
    else: data = data.replace(b"subject\t", b"")
    with pytest.raises(ValueError):
        copy_block(io.BytesIO(data), io.BytesIO(), entry, hashlib.sha256())


def test_late_inventory_failure_keeps_partial(tmp_path):
    source = tmp_path/"source"
    source.write_bytes(row("a"))
    def plan():
        yield block("a")
        raise ValueError("Late exhausted-inventory check failed")
    with pytest.raises(ValueError):
        merge(plan(), {"old": source}, tmp_path/"merged")
    assert (tmp_path/"merged/all.blast.partial").read_bytes() == row("a")
    assert not (tmp_path/"merged/all.blast.candidate").exists()


def test_source_change_during_merge_rejected(tmp_path):
    source = tmp_path/"source"
    source.write_bytes(row("a"))
    def plan():
        yield block("a")
        with source.open("ab") as stream:
            stream.write(b"changed")
    with pytest.raises(ValueError, match="Source changed"):
        merge(plan(), {"old": source}, tmp_path/"merged")
    assert not (tmp_path/"merged/all.blast.candidate").exists()


def test_repeated_query_rejected(tmp_path):
    source = tmp_path/"source"
    source.write_bytes(row("a"))
    with pytest.raises(ValueError, match="Repeated query"):
        merge([block("a"), block("a")], {"old": source}, tmp_path/"merged")
    assert not (tmp_path/"merged/all.blast.candidate").exists()


def test_fsync_failure_cannot_publish_candidate(tmp_path, monkeypatch):
    source = tmp_path/"source"
    source.write_bytes(row("a"))
    def fail(_):
        raise OSError("test fsync failure")
    monkeypatch.setattr("benchmark_tools.merge_blast_recovery_blocks.os.fsync", fail)
    with pytest.raises(OSError):
        merge([block("a")], {"old": source}, tmp_path/"merged")
    assert (tmp_path/"merged/all.blast.partial").exists()
    assert not (tmp_path/"merged/all.blast.candidate").exists()
