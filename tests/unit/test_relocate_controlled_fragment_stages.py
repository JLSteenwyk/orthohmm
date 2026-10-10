"""Relocation fixtures reuse exact frozen readers, without scientific reruns."""

from collections import Counter
import csv
import json
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import pytest

from benchmark_tools import readback_controlled_fragment_stage_table as reader
from benchmark_tools import relocate_controlled_fragment_stages as portable
from benchmark_tools.fragment_trace_artifact_access import ArtifactAccess


def put(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data if isinstance(data, bytes) else data.encode("ascii"))
    return portable.pin(path)


def dump(path, data):
    return put(path, portable.json_bytes(data))


@pytest.fixture
def panel(tmp_path):
    root = tmp_path / "original"
    source_root = Path(portable.__file__).resolve().parents[1]
    for name in (*portable.KERNELS, "fragment_trace_artifact_access", "relocate_controlled_fragment_stages", "__init__"):
        put(root / "benchmark_tools" / (name + ".py"), (source_root / "benchmark_tools" / (name + ".py")).read_bytes())
    input_dir = root / "input"
    inputs = [put(input_dir / (species + ".fasta"), ">" + gene + "\nACDEFGHIK\n")
              for species, gene in (("x", "a"), ("y", "b"))]
    execution_path = root / "execution.json"
    runs, outputs = {}, {}
    for method in reader.previous.METHODS[:2]:
        output = root / method
        outputs[method] = output
        working = output / "orthohmm_working_res"
        checkpoint = working / "high_sensitivity_checkpoint"
        refs = [put(output / "orthohmm_orthogroups.txt", "OG0: a b\n"),
                put(working / "orthohmm_edges.txt", "a\tb\t1\n")]
        cp_refs = {"gene_names.txt": put(checkpoint / "gene_names.txt", "a\nb\n")}
        for name, array in (("gene_to_species", [0, 1]), ("hit_queries", [0, 1]),
                            ("hit_targets", [1, 0]), ("hit_scores", [42., 42.])):
            np.save(checkpoint / (name + ".npy"), np.array(array))
            cp_refs[name + ".npy"] = portable.pin(checkpoint / (name + ".npy"))
        refs.extend(cp_refs.values())
        refs.append(dump(checkpoint / "manifest.json", dict(genes=2, hits=2, files=cp_refs)))
        if method == reader.previous.METHODS[1]:
            phy = output / "orthohmm_phylogeny"
            refs.extend([
                put(phy / "orthohmm_pairwise_orthologs.tsv", "gene_a\tgene_b\tspecies_a\tspecies_b\na\tb\tx\ty\n"),
                put(working / "phylogeny_candidate_superfamilies.txt", "a b\n"),
                put(phy / "orthohmm_root_hogs.tsv", "root_hog\tgenes\tsource_family\nR0\ta,b\tFamily0000000\n"),
                dump(phy / "provenance_manifest.json", dict(pair_orthology_rule="positive_paralogy",
                     root_duplication_rule="species_overlap", membership_reconciliation=None)),
                dump(working / "phylogeny_candidate_merges.json", []),
                put(phy / "orthohmm_reconciliation_nodes.tsv", "node_id\tsource_family\n"),
            ])
        runs[method] = dict(status="process_succeeded", exit_code=0, argv=["python", "run.py", str(input_dir)],
                            outputs=[dict(r, absolute_path=r["path"]) for r in refs])
    of = root / "orthofinder_full"
    of_input = of / "input"
    copied_inputs = [put(of_input / Path(r["path"]).name, Path(r["path"]).read_bytes()) for r in inputs]
    native = of_input / "OrthoFinder" / "Results_test"
    native_refs = [put(native / "WorkingDirectory" / "SequenceIDs.txt", "0_0: a\n1_0: b\n"),
                   put(native / "WorkingDirectory" / "clusters_OrthoFinder_I1.5.txt_id_pairs.txt",
                       "(mclmatrix\nbegin\n0 0_0 1_0 $\n)\n"),
                   put(native / "Orthologues" / "Orthologues_x" / "x__v__y.tsv", "Orthogroup\tx\ty\nOG0\ta\tb\n")]
    runs["orthofinder_full"] = dict(status="process_succeeded", exit_code=0,
                                    argv=["orthofinder", "-f", str(of_input)],
                                    outputs=[dict(r, absolute_path=r["path"]) for r in native_refs + copied_inputs])
    outputs["orthofinder_full"] = outputs[portable.CHECKPOINT] = of
    execution_ref = dump(execution_path, dict(methods=runs))
    binding = {}
    for method in reader.previous.METHODS:
        config = dict(argv=["python", "run.py", str(input_dir)]) if method in reader.previous.METHODS[:2] else (
            dict(copy_inputs_from=str(input_dir)) if method == "orthofinder_full" else dict(parent_method="orthofinder_full"))
        binding[method] = dict(execution=execution_ref, output=str(outputs[method]), configured=config)
    cases = [dict(case_id="Case" + str(i), method=method, seed=1, fragment_endpoints=0, category="TP", truth=True,
                  gene_a="a", gene_b="b", comparator_predictions={arm: dict.fromkeys(reader.previous.METHODS, True)
                                                                for arm in reader.previous.ARMS})
             for i, method in enumerate(reader.previous.METHODS)]
    selection_ref = dump(root / "selection.json", dict(cases=cases, bindings=[dict(seed=1, arms=dict.fromkeys(reader.previous.ARMS, binding))],
                                                      checked_inputs=inputs + [execution_ref]))
    table, counts, staged = [], Counter(), []
    for case in cases:
        stages = {}
        for arm in reader.previous.ARMS:
            row = reader.previous.observation(case["method"], reader.previous.read_context(binding[case["method"]], case["method"], {("a", "b")}), ("a", "b"))
            stages[arm] = row
            counts[case["method"], case["category"], row["observed_location"], arm] += 1
            table.append(dict(case_id=case["case_id"], method=case["method"], seed=1, fragment_endpoints=0,
                              category="TP", truth=True, arm=arm, **row))
        staged.append(dict(case, stages=stages))
    copied_paths = {r["path"] for r in copied_inputs}
    refs = inputs + [execution_ref] + [r for run in runs.values() for r in run["outputs"] if r["path"] not in copied_paths]
    report_ref = dump(root / "report.json", dict(schema="controlled_fragment_stage_trace_v1", status="retained_stages_verified",
                      new_inference_or_scoring=False, new_bootstrap_draws=0, publication_ready=False, uncertainty_admitted=False,
                      scientific_timings_admitted=False, independent_confirmation=False, selection=selection_ref,
                      checked_inputs=refs, cases=staged, stage_rows=len(table),
                      summary=[dict(method=m, category=c, observed_location=l, arm=a, representatives=n)
                               for (m, c, l, a), n in sorted(counts.items())]))
    with (root / "stages.tsv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, reader.FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(reader.serialized(r) for r in table)
    expected = reader.verify(Path(report_ref["path"]), report_ref["sha256"])
    # Source labels belong to the synthetic original tree; source bytes are unchanged.
    expected["source"] = portable.pin(root / "benchmark_tools" / (portable.KERNELS[2] + ".py"))
    expected["previous_reader"] = portable.pin(root / "benchmark_tools" / (portable.KERNELS[1] + ".py"))
    readback_ref = dump(root / "readback.json", expected)
    return root, Path(report_ref["path"]), Path(readback_ref["path"]), readback_ref["sha256"]


def prepare_panel(panel, output):
    return portable.prepare(*panel, output)


def test_plan_adds_only_inventoried_checkpoint_fastas(panel):
    planned = portable.plan(*panel)
    fastas = [r for r in planned["artifacts"] if "checkpoint_copied_fasta" in r["roles"]]
    assert len(fastas) == planned["added_checkpoint_fastas"] == 2
    assert {Path(r["logical_path"]).name for r in fastas} == {"x.fasta", "y.fasta"}


def test_relocated_subprocess_replays_all_methods_after_original_removed(panel, tmp_path):
    prepared = prepare_panel(panel, tmp_path / "component")
    moved = tmp_path / "elsewhere" / "component"
    moved.parent.mkdir()
    shutil.move(tmp_path / "component", moved)
    shutil.rmtree(panel[0])
    output = tmp_path / "readback.json"
    subprocess.run([sys.executable, "-I", "-B", str(moved / "runner/benchmark_tools/relocate_controlled_fragment_stages.py"),
                    "replay", "--component", str(moved), "--manifest-sha256", prepared["manifest"]["sha256"],
                    "--output", str(output)], cwd=moved.parent, check=True, capture_output=True, text=True)
    result = json.loads(output.read_text())
    assert result["status"] == "relocated_native_readback_matches_retained_result"
    assert result["native_readback"]["stage_rows"] == 8
    assert result["native_readback"]["native_contexts"] == 8
    assert not result["original_path_fallback"] and not result["historical_inference_reproduced"]


def test_private_loading_restores_existing_package_and_readers(panel, tmp_path):
    prepared = prepare_panel(panel, tmp_path / "component")
    manifest = json.loads(Path(prepared["manifest"]["path"]).read_text())
    access = ArtifactAccess(tmp_path / "component", manifest["artifacts"])
    before = {key: sys.modules.get(key) for key in ["benchmark_tools"] + ["benchmark_tools." + n for n in portable.KERNELS]}
    copied = portable.load_readers(access, manifest["sources"])
    assert copied is not reader and copied.previous is not reader.previous
    assert reader.Path is Path and reader.previous.Path is Path
    assert before == {key: sys.modules.get(key) for key in before}


def test_private_loading_also_restores_modules_on_import_failure(panel, tmp_path):
    prepared = prepare_panel(panel, tmp_path / "component")
    manifest = json.loads(Path(prepared["manifest"]["path"]).read_text())
    row = next(r for r in manifest["artifacts"] if r["logical_path"] == manifest["sources"]["reader"])
    path = tmp_path / "component" / row["target"]
    path.write_text("raise RuntimeError('synthetic import failure')\n")
    row.update(portable.identity(path))
    before = {key: sys.modules.get(key) for key in ["benchmark_tools"] + ["benchmark_tools." + n for n in portable.KERNELS]}
    with pytest.raises(RuntimeError, match="synthetic import failure"):
        portable.load_readers(ArtifactAccess(tmp_path / "component", manifest["artifacts"]), manifest["sources"])
    assert before == {key: sys.modules.get(key) for key in before}


def test_native_result_must_equal_expected_not_just_pass_file_hashes(panel, tmp_path):
    prepared = prepare_panel(panel, tmp_path / "component")
    manifest_path = Path(prepared["manifest"]["path"])
    manifest = json.loads(manifest_path.read_text())
    row = next(r for r in manifest["artifacts"] if "expected_native_readback" in r["roles"])
    path = tmp_path / "component" / row["target"]
    expected = json.loads(path.read_text())
    expected["stage_rows"] += 1
    path.write_bytes(portable.json_bytes(expected))
    old_size = row["bytes"]
    row.update(portable.identity(path))
    manifest["artifact_payload_bytes"] += row["bytes"] - old_size
    manifest_path.write_bytes(portable.json_bytes(manifest))
    with pytest.raises(ValueError, match="differs from retained result"):
        portable.replay(tmp_path / "component", portable.identity(manifest_path)["sha256"])


@pytest.mark.parametrize("target", ["manifest", "native", "kernel", "runtime"])
def test_changed_manifest_payload_or_source_is_refused(panel, tmp_path, target):
    prepared = prepare_panel(panel, tmp_path / "component")
    manifest_path = Path(prepared["manifest"]["path"])
    manifest = json.loads(manifest_path.read_text())
    if target == "manifest":
        manifest_path.write_text(manifest_path.read_text() + "\n")
    else:
        role = dict(native="stage_input", kernel="unchanged_kernel", runtime="relocation_runtime")[target]
        row = next(r for r in manifest["artifacts"] if role in r["roles"])
        (tmp_path / "component" / row["target"]).write_text("changed")
    with pytest.raises(ValueError, match="Changed|differs"):
        portable.replay(tmp_path / "component", prepared["manifest"]["sha256"])


def test_changed_expected_readback_is_refused_before_output(panel, tmp_path):
    panel[2].write_text(panel[2].read_text() + "\n")
    with pytest.raises(ValueError, match="Changed expected"):
        prepare_panel(panel, tmp_path / "component")
    assert not (tmp_path / "component").exists()


def test_uninventoried_checkpoint_fasta_is_refused(panel, tmp_path):
    put(panel[0] / "orthofinder_full/input/z.fasta", ">z\nAAAA\n")
    with pytest.raises(ValueError, match="execution inventory"):
        prepare_panel(panel, tmp_path / "component")
    assert not (tmp_path / "component").exists()


def test_existing_outputs_are_preserved(panel, tmp_path):
    output = tmp_path / "component"
    output.mkdir()
    (output / "canary").write_text("preserve")
    with pytest.raises(ValueError, match="Existing"):
        prepare_panel(panel, output)
    assert (output / "canary").read_text() == "preserve"


def test_payload_bound_is_checked_before_copy(panel, tmp_path, monkeypatch):
    monkeypatch.setattr(portable, "MAX_BYTES", 1)
    with pytest.raises(ValueError, match="bounded"):
        prepare_panel(panel, tmp_path / "component")
    assert not (tmp_path / "component").exists()
