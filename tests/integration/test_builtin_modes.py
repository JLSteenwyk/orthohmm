"""Exercise application modes without research drivers or benchmark data."""

import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from orthohmm.search.msa_center_star import _load_pair_align
from tests.support import fixture


ROOT = Path(__file__).resolve().parents[2]
pytestmark = pytest.mark.integration


def run_fixture(tmp_path, options):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    proteins = fixture()
    for species, sequences in proteins.items():
        (inputs / (species + ".faa")).write_text(
            "".join(f">{gene}\n{sequence}\n" for gene, sequence in sequences.items())
        )
    before = {path: path.read_bytes() for path in inputs.iterdir()}
    output = tmp_path / "output"
    output.mkdir()
    metrics = tmp_path / "metrics.json"
    command = [sys.executable, "-m", "orthohmm", str(inputs), "-o", str(output),
               "-c", "1", "--threads_per_worker", "1", "--metrics_json", str(metrics)]
    for option, value in options.items():
        command.extend(["--" + option, str(value)])
    result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, timeout=180)
    assert result.returncode == 0, result.stdout + result.stderr
    assigned = []
    for line in (output / "orthohmm_orthogroups.txt").read_text().splitlines():
        group, members = line.split(":", 1)
        assert group and members.split()
        assigned.extend(members.split())
    assert sorted(assigned) == sorted(gene for genes in proteins.values() for gene in genes)
    assert before == {path: path.read_bytes() for path in inputs.iterdir()}
    report = json.loads(metrics.read_text())
    assert report["status"] == "complete"
    for name, value in options.items():
        if name in report["metadata"]:
            assert report["metadata"][name] == value
    return output, report


@pytest.fixture
def native_alignment():
    try:
        _load_pair_align()
    except RuntimeError:
        pytest.skip("Native pair_align kernel required; run make kernels")


@pytest.mark.parametrize("options", [
    {},
    {"profile_sequence_weighting": "henikoff"},
    {"profile_projection_policy": "residue_preserving"},
    {"profile_projection_policy": "residue_preserving", "profile_sequence_weighting": "henikoff_cluster"},
    {"profile_alignment_policy": "strict_span_guard"},
    {"search_candidate_policy": "protected_union", "search_result_order": "canonical",
     "search_max_candidates_per_query": 8},
])
def test_high_sensitivity_modes(tmp_path, native_alignment, options):
    _, report = run_fixture(tmp_path, dict(accuracy_profile="high_sensitivity", **options))
    assert report["metadata"]["accuracy_profile"] == "high_sensitivity"


@pytest.mark.parametrize("tree_mode", ["supplied", "infer"])
def test_phylogeny_cli(tmp_path, native_alignment, tree_mode):
    pytest.importorskip("dendropy")
    aligner = shutil.which("mafft")
    tree_builder = shutil.which("FastTree") or shutil.which("fasttree")
    if not aligner or not tree_builder:
        pytest.skip("Phylogeny integration requires MAFFT and FastTree")
    options = dict(accuracy_profile="high_sensitivity", phylogeny="reconcile",
                   species_tree_mode=tree_mode, phylogeny_candidates="satellite_v2",
                   phylogeny_root_rule="species_overlap", phylogeny_pair_rule="positive_paralogy",
                   aligner=aligner, tree_builder=tree_builder)
    if tree_mode == "supplied":
        tree = tmp_path / "species.nwk"
        tree.write_text("[&R] ((S1:1,S2:1):1,(S3:1,S4:1):1);\n")
        options["species_tree"] = str(tree)
    else:
        options["species_tree_rooting"] = "min_variance"
    output, report = run_fixture(tmp_path, options)
    assert report["metadata"]["phylogeny"] == "reconcile"
    summary = json.loads((output / "orthohmm_phylogeny/reconciliation_summary.json").read_text())
    assert summary["candidate_families"] > 0
    assert summary["reconciled_families"] > 0
