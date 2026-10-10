from copy import deepcopy
import json
from pathlib import Path

import pytest

from benchmark_tools import integrate_conditional_fas_manuscript as current


def inputs():
    root=Path(__file__).resolve().parents[2]/"benchmark_tools/results"
    return [json.loads((root/current.PINS[key][0]).read_text()) for key in ("report","execution")]


def test_scoped_insertions_restore_entire_parent_and_keep_separate_targets():
    report,execution=inputs()
    parent="Old evidence\n"+current.ANCHORS[0]+"prior results\n"+current.ANCHORS[1]+"prior limits\n"
    revised,inserts=current.manuscript(parent,report,execution)
    restored=revised
    for section in inserts: restored=restored.replace(section,"",1)
    assert restored==parent and len(inserts)==2
    assert "September eight-method" in revised and "not the later four-cell" in revised
    assert "not biological error bars" in revised and "Proteinortho exceeds both" in revised
    assert "newer\nnative factorial" in revised
    assert len([line for line in inserts[0].splitlines() if line.startswith("| ")])==10


@pytest.mark.parametrize("change",["methods","coverage","historical","partial","contrast","terminal"])
def test_changed_data_scope_or_observed_terminal_refused(change):
    report,execution=deepcopy(inputs())
    if change=="methods": report["methods"].reverse()
    elif change=="coverage": report["component_error"]*=2
    elif change=="historical": report["unconditional_historical_interval_admission"]=True
    elif change=="partial": report["methods"][0]["interval"]=None;report["all_eight_ranges_computed"]=False
    elif change=="contrast": report["contrasts"][0]["conditional_expected_difference_bounds"][0]+=.01
    else: execution["document_checks"][-1]["result"]["exit_code"]=1
    with pytest.raises(ValueError): current.sections(report,execution)


@pytest.mark.parametrize("parent",["missing",current.ANCHORS[0]*2+current.ANCHORS[1]])
def test_missing_or_ambiguous_anchors_refused(parent):
    with pytest.raises(ValueError): current.manuscript(parent,*inputs())


def test_all_three_integration_inputs_are_pinned():
    root=Path(__file__).resolve().parents[2]/"benchmark_tools/results"
    for _,(name,sha) in current.PINS.items(): assert current.record(root/name)["sha256"]==sha


def test_occupied_destinations_refused_without_generation(tmp_path):
    output=tmp_path/"benchmark_tools/results"
    output.mkdir(parents=True)
    existing=output/"exists.md";existing.touch()
    with pytest.raises(FileExistsError): current.run(tmp_path,existing,output/"new.json")
