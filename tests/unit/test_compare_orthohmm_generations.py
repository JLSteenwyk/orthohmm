from copy import deepcopy
import math
from pathlib import Path

import pytest

from benchmark_tools import compare_orthohmm_generations as timing


def panel():
    return dict(runs=[dict(r, status="validated", elapsed_s=10+r["index"], search_s=5+r["index"])
                      for r in timing.schedule()])


def test_three_by_three_schedule_balances_each_order_position():
    rows = timing.schedule()
    assert len(rows) == 9
    assert [r["index"] for r in rows] == list(range(9))
    for j in range(3):
        assert {rows[r*3+j]["method"] for r in range(3)} == set(timing.METHODS)
    with pytest.raises(ValueError):
        timing.schedule(2)


def test_summary_keeps_every_observation_and_uses_median_not_fastest():
    data = panel()
    original = deepcopy(data)
    rows = timing.summarize(data)
    assert rows[0]["elapsed_s"] == dict(median=15, minimum=10, maximum=17, values=[10.,15.,17.])
    assert rows[0]["original_to_modern_elapsed_ratio"] == 1
    assert data == original


@pytest.mark.parametrize("kind", ["missing", "reordered", "failed", "timeout", "zero", "negative", "nan", "infinite"])
def test_failed_incomplete_or_invalid_measurements_cannot_produce_speed_ratios(kind):
    data = panel()
    if kind == "missing": data["runs"].pop()
    elif kind == "reordered": data["runs"].reverse()
    elif kind in {"failed", "timeout"}: data["runs"][0]["status"] = kind
    else:
        data["runs"][0]["elapsed_s"] = {"zero":0, "negative":-1, "nan":math.nan, "infinite":math.inf}[kind]
    with pytest.raises(ValueError):
        timing.summarize(data)


def test_partition_requires_exact_protein_ownership_not_only_success_exit_code(tmp_path):
    path = tmp_path / "groups.txt"
    path.write_text("OG0: a b\nOG1: c\n")
    assert timing.partition(path, ["a", "b", "c"])["proteins"] == 3
    for text in ("", "OG0: a a b c\n", "OG0: a b\n", "OG0: a b c d\n", "OG0:\nOG1: a b c\n"):
        path.write_text(text)
        with pytest.raises(ValueError): timing.partition(path, ["a", "b", "c"])


def test_original_search_requires_all_directional_completed_tables(tmp_path):
    for left in ("a.faa", "b.faa"):
        for right in ("a.faa", "b.faa"):
            (tmp_path / f"{left}_2_{right}.phmmerout.txt").write_text("hit\n# [ok]\n")
    timing.search_complete(tmp_path, ["a.faa", "b.faa"])
    (tmp_path / "a.faa_2_b.faa.phmmerout.txt").write_text("hit\n")
    with pytest.raises(ValueError): timing.search_complete(tmp_path, ["a.faa", "b.faa"])


def test_legacy_commands_keep_original_nested_parallelism_inside_fixed_affinity(tmp_path):
    plan = dict(python=dict(path="/python"), executables={k:dict(path="/"+k) for k in ("taskset", "phmmer", "mcl", "mafft", "FastTree")}, cpu_set=list(range(64,80)))
    _, cwd, command = timing.command(tmp_path, dict(index=0, method="original"), plan)
    assert command[:3] == ["/taskset", "-c", ",".join(map(str,range(64,80)))]
    assert command[command.index("-c",3)+1] == "4"
    assert cwd == tmp_path / "sources/original"
    assert "--start" not in command
    for method in ("sensitive", "phylogeny"):
        _, cwd, command = timing.command(tmp_path, dict(index=1, method=method), plan)
        assert command[command.index("-c",3)+1] == "16"
        assert command[command.index("--threads_per_worker")+1] == "4"
        assert "--start" not in command
        assert ("--phylogeny" in command) == (method == "phylogeny")


def test_fasta_identifier_compatibility_is_checked_for_both_methods(tmp_path):
    path = tmp_path / "a.faa"
    path.write_text(">one\nACD\n>two\nEFG\n")
    assert timing.sequences(tmp_path)[0][1] == ["one", "two"]
    path.write_text(">" + "a"*51 + "\nACD\n")
    with pytest.raises(ValueError): timing.sequences(tmp_path)
