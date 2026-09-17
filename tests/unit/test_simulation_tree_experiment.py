from copy import deepcopy

import pytest

from benchmark_tools.run_simulation_tree_experiment import (
    PANEL_SHA, METHODS, VARIANTS, require_mode_equivalence, select_tree, run,
)
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def mode_fixture():
    rows, records = [], []
    for index in range(70):
        unavailable = {METHODS[1]: {"status": "failed", "reason": "original failure"}} if index < 8 else {}
        rows.append({"label": str(index), "unavailable": unavailable})
        for method in METHODS:
            row = {"label": str(index), "method": method,
                   "status": "unavailable" if method in unavailable else "equivalent"}
            if method in unavailable:
                row["baseline_failure"] = unavailable[method]
            records.append(row)
    admission = {"status": "mode_panel_verified_unscored", "accuracy_evaluated": False,
                 "panel": {"sha256": PANEL_SHA}, "records": records}
    return admission, {"rows": rows}


@pytest.mark.parametrize("problem", [None, "running", "non_equivalence", "new_failure", "missing", "duplicate", "unavailable_reason"])
def test_all_available_modes_must_pass_before_main_experiment(problem):
    admission, panel = mode_fixture()
    if problem == "running":
        admission["status"] = "validating"
    elif problem == "non_equivalence":
        admission["records"][-1]["status"] = "not_equivalent"
    elif problem == "new_failure":
        admission["records"][-1]["status"] = "failed"
    elif problem == "missing":
        admission["records"].pop()
    elif problem == "duplicate":
        admission["records"][-1] = deepcopy(admission["records"][0])
    elif problem == "unavailable_reason":
        admission["records"][1]["baseline_failure"] = {"reason": "different"}
    if problem:
        with pytest.raises(ValueError):
            require_mode_equivalence(admission, panel)
    else:
        require_mode_equivalence(admission, panel)
        assert len(admission["records"]) == 140


def tree_fixture(tmp_path):
    original, plain = tmp_path / "original", tmp_path / "plain"
    text = "((a:1,b:1):1,(c:1,d:1):1):0;\n"
    original.write_text("[&R] " + text)
    plain.write_text(text)
    portable, prepared = {"trees": []}, {"datasets": []}
    for index in range(70):
        dataset = {"condition": str(index), "seed": 1, "taxa": ["a", "b", "c", "d"], "variants": []}
        for variant in VARIANTS:
            dataset["variants"].append({"label": variant, "tree": record(original)})
            portable["trees"].append({"condition": str(index), "seed": 1, "variant": variant,
                "taxa": dataset["taxa"], "source_tree": record(original), "tree": record(plain)})
        prepared["datasets"].append(dataset)
    return portable, prepared


@pytest.mark.parametrize("index", [0, 1, 209])
def test_fixed_tree_index_checks_source_and_portable_identity(tmp_path, index):
    portable, prepared = tree_fixture(tmp_path)
    assert select_tree(portable, prepared, index) == portable["trees"][index]


@pytest.mark.parametrize("problem", ["missing", "duplicate", "taxa", "content", "topology", "length"])
def test_changed_tree_manifest_or_tree_is_rejected(tmp_path, problem):
    portable, prepared = tree_fixture(tmp_path)
    if problem == "missing":
        portable["trees"].pop()
    elif problem == "duplicate":
        portable["trees"][-1] = deepcopy(portable["trees"][0])
    elif problem == "taxa":
        portable["trees"][0]["taxa"] = ["a", "b", "c", "e"]
    else:
        path = tmp_path / "plain"
        path.write_text("((a:1,c:1):1,(b:1,d:1):1):0;" if problem == "topology"
                        else "((a:2,b:1):1,(c:1,d:1):1):0;")
        if problem != "content":
            portable["trees"][0]["tree"] = record(path)
    with pytest.raises(ValueError):
        select_tree(portable, prepared, 0)


def test_existing_run_not_restarted(tmp_path):
    with pytest.raises(FileExistsError):
        run(tmp_path, 0, "unused", tmp_path)
