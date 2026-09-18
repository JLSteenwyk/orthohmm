import json
from pathlib import Path

import pytest

from benchmark_tools.run_dgx_hierarchy_native_smoke import ROOT, SPEC_SHA, relocate, launch, read_pinned


def test_only_output_prefix_changes_in_complete_frozen_spec():
    root = Path(__file__).resolve().parents[2]
    spec = read_pinned(root / "benchmark_tools/results/dgx_launcher_smoke_spec_20260917.json", SPEC_SHA)
    old, new = str(ROOT / "launcher_smoke_v1"), str(ROOT / "hierarchy_native_smoke_v1")
    for row in spec["runs"]:
        observed = relocate(row)
        assert json.loads(json.dumps(observed).replace(new + "/", old + "/")) == row
        assert observed["cwd"] == row["cwd"]
        assert observed["environment_role"] == row["environment_role"]
    assert relocate(old + "_other/file") == old + "_other/file"


@pytest.mark.parametrize("index", [-1, 3, 9])
def test_invalid_method_index_rejected_before_execution(index):
    with pytest.raises(ValueError, match="three frozen"):
        launch(Path("missing"), index, Path("missing"), "bad")


def test_transferred_recipe_matches_local_sources():
    import hashlib

    root = Path(__file__).resolve().parents[2] / "benchmark_tools"
    recipe = json.loads((root / "results/dgx_hierarchy_native_recipe_v1_20260918.json").read_text())
    files = [r for r in recipe["records"] if r["kind"] == "file"]
    assert len(files) == 24
    assert "__init__.py" in {Path(r["path"]).name for r in files}
    for row in files:
        name = Path(row["path"]).name
        path = root / name
        if not path.exists():
            path = root / "results" / name
        assert len(path.read_bytes()) == row["bytes"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == row["sha256"]
