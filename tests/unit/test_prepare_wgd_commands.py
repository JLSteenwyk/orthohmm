from pathlib import Path

from benchmark_tools.prepare_wgd_commands import commands


def panel():
    return commands(*map(Path, ("/core", "/input", "/out", "/python", "/of", "/sonic")))


def test_four_methods_with_preserved_search_parameters():
    rows = panel()
    assert len(rows) == 4
    high, satellite = rows[:2]
    assert high["argv"][5:] == ["--cpu", "32", "--threads-per-worker", "8",
                               "--accuracy-profile", "high_sensitivity"]
    assert satellite["argv"][5:11] == high["argv"][5:]
    assert satellite["argv"][-2:] == ["--species-tree-rooting", "min_variance"]
    assert satellite["output_semantics"] == "root_hogs"


def test_comparators_use_distinct_fresh_inputs_and_native_defaults():
    of, sonic = panel()[2:]
    assert of["argv"] == ["/of", "-f", "/out/orthofinder/input", "-t", "32", "-a", "8", "-S", "diamond"]
    assert "-og" not in of["argv"]
    assert sonic["argv"] == ["/sonic", "-i", "/out/sonicparanoid/input", "-o", "/out/sonicparanoid/output", "-t", "32"]
    assert of["copy_inputs_to"] != sonic["copy_inputs_to"]
