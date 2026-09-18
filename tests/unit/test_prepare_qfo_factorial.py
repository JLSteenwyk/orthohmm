from copy import deepcopy

import pytest

from benchmark_tools.prepare_qfo_factorial import select_seeds, verify_inputs
from benchmark_tools.prepare_orthobench_factorial import plan_cells


def stages():
    return [{"label": name, "output": {"path": name, "sha256": str(i)}} for i, name in enumerate(
        ("multipass", "multipass_refined", "strict_profiles", "strict_profiles_refined"))]


def test_refined_seed_selection_preserves_identity():
    rows = stages()
    assert select_seeds(rows) == [(False, rows[1]["output"]), (True, rows[3]["output"])]
    assert select_seeds(rows)[0][1] is rows[1]["output"]


@pytest.mark.parametrize("kind", ["missing", "reordered", "renamed"])
def test_wrong_stage_inventory_rejected(kind):
    rows = stages()
    if kind == "missing":
        rows.pop()
    elif kind == "reordered":
        rows.reverse()
    else:
        rows[1]["label"] = "multipass"
    with pytest.raises(ValueError):
        select_seeds(rows)


def test_fasta_inventory_requires_exact_content_and_unique_files():
    rows = [{"path": f"/{i}.fasta", "sha256": str(i), "bytes": i + 1} for i in range(78)]
    verify_inputs(rows, rows[::-1])
    changed = deepcopy(rows)
    changed[0]["sha256"] = "changed"
    with pytest.raises(ValueError):
        verify_inputs(changed, rows)
    with pytest.raises(ValueError):
        verify_inputs(rows[:-1], rows)
    with pytest.raises(ValueError):
        verify_inputs([rows[0]] * 78, [rows[0]] * 78)


def test_qfo_plan_keeps_eight_cells_and_constraints_only_on_expanded_arms(tmp_path):
    rows = plan_cells(tmp_path / "out", tmp_path / "fasta", tmp_path / "launcher.py", 32)
    assert len(rows) == 8
    assert len({r["label"] for r in rows}) == 8
    for row in rows:
        if row["reconciliation"]:
            assert ("--membership-constraints" in row["argv"]) == row["candidate_expansion"]
            assert row["argv"][row["argv"].index("--species-tree-mode") + 1] == "infer"
            assert row["argv"][row["argv"].index("--cpu") + 1] == "32"
        else:
            assert "argv" not in row
            assert row["prediction"] == row["candidate_partition"]
