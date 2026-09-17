import json
from pathlib import Path

import pytest

from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.validate_scaling_outputs import validate
from benchmark_tools.gnu_time_companion import command as time_command


def fixture(tmp_path, method="orthohmm_satellite_v2"):
    source = tmp_path / "source"
    source.mkdir()
    for species, gene in (("a", "g1"), ("b", "g2")):
        (source / (species + ".fa")).write_text(">" + gene + "\nAAA\n")
    inputs = [record(p) for p in sorted(source.iterdir())]
    output = tmp_path / "output"
    output.mkdir()
    config = {"output": str(output), "metrics": str(tmp_path / "metrics.json")}
    cwd = str(tmp_path / "core")
    if method == "orthofinder_full":
        copies = output / "input"
        copies.mkdir()
        for p in source.iterdir():
            (copies / p.name).write_bytes(p.read_bytes())
        config["copy_inputs_to"] = str(copies)
        argv = ["orthofinder", "-f", str(copies), "-t", "32", "-a", "32"]
        (output / "Log.txt").write_text("2026 : Started OrthoFinder version 3.1.5\nCommand Line: " + " ".join(argv) + "\n2026 : OrthoFinder run completed\n")
        (output / "OrthoFinder_graph.txt").write_text("0 1:0.9 $\n")
        (output / "SequenceIDs.txt").write_text("0_0: g1\n1_0: g2\n")
        (output / "clusters_OrthoFinder_I1.2.txt_id_pairs.txt").write_text("begin\n0 0_0 1_0 $\n)\n")
        for a, b, ga, gb in (("a", "b", "g1", "g2"), ("b", "a", "g2", "g1")):
            folder = output / "Orthologues" / ("Orthologues_" + a)
            folder.mkdir(parents=True)
            (folder / (a + "__v__" + b + ".tsv")).write_text(f"Orthogroup\t{a}\t{b}\nOG1\t{ga}\t{gb}\n")
    else:
        argv = ["python", "-m", "orthohmm", str(source), "-o", str(output)]
        metrics = {"status": "complete", "cwd": cwd,
                   "command": ["python", str(Path(cwd) / "orthohmm/__main__.py"), *argv[3:]],
                   "counts": {"genes": 2, "species": 2, "orthogroups": 1, "phylogeny_root_hogs": 1, "phylogeny_ortholog_pairs": 1}}
        Path(config["metrics"]).write_text(json.dumps(metrics))
        (output / "orthohmm_orthogroups.txt").write_text("OG1: g1 g2\n")
        phylo = output / "orthohmm_phylogeny"
        phylo.mkdir()
        (phylo / "orthohmm_root_hogs.tsv").write_text("root_hog\tsource_family\tgenes\nHOG1\tOG1\tg1,g2\n")
        (phylo / "orthohmm_pairwise_orthologs.tsv").write_text("gene_a\tspecies_a\tgene_b\tspecies_b\ng1\ta\tg2\tb\n")
    run = {"configuration": config, "native_argv": argv, "cwd": cwd, "native_method": method,
           "dataset": {"inputs": inputs, "proteomes": 2, "proteins": 2}}
    measurement = {"command": argv, "cwd": cwd, "exit_code": 0, "timed_out": False}
    return run, measurement


@pytest.mark.parametrize("method", ["orthohmm_high_sensitivity", "orthohmm_satellite_v2", "orthofinder_full"])
def test_native_outputs_without_harness_provenance(tmp_path, method):
    result = validate(*fixture(tmp_path, method))
    assert result["status"] == "native_scaling_outputs_checked"
    assert result["input_genes"] == 2 and result["checked_files"]
    assert not result["accuracy_evaluated"] and not result["resource_measurements_admitted"]


@pytest.mark.parametrize("method", ["orthohmm_high_sensitivity", "orthohmm_satellite_v2", "orthofinder_full"])
def test_explicit_gnu_time_keeps_native_command_checks(tmp_path, method):
    run, measurement = fixture(tmp_path, method)
    output = tmp_path / "time.tsv"
    run["gnu_time"] = {"executable": "/usr/bin/time", "output": str(output)}
    measurement["command"] = time_command(run["native_argv"], output)
    output.write_text("elapsed_seconds\t1\nuser_seconds\t0\nsystem_seconds\t0\nmax_process_rss_kib\t1\nexit_status\t0\n")
    result = validate(run, measurement)
    assert result["gnu_time_companion"]["wrapper_in_collector_wall_time"]
    measurement["command"][-1] = "changed_native_argument"
    with pytest.raises(ValueError, match="Measured command"):
        validate(run, measurement)


def test_gnu_time_exit_disagreement(tmp_path):
    run, measurement = fixture(tmp_path)
    output = tmp_path / "time.tsv"
    run["gnu_time"] = {"executable": "/usr/bin/time", "output": str(output)}
    measurement["command"] = time_command(run["native_argv"], output)
    output.write_text("elapsed_seconds\t1\nuser_seconds\t0\nsystem_seconds\t0\nmax_process_rss_kib\t1\nexit_status\t7\n")
    with pytest.raises(ValueError, match="statuses disagree"):
        validate(run, measurement)


@pytest.mark.parametrize("problem", ["command", "exit", "timeout", "missing_gene", "duplicate_gene", "pair_species", "duplicate_pair", "metrics_command", "metrics_count"])
def test_orthohmm_failure_and_corruption_rejected(tmp_path, problem):
    run, measurement = fixture(tmp_path)
    output = Path(run["configuration"]["output"])
    if problem == "command":
        measurement["command"] = ["other"]
    elif problem == "exit":
        measurement["exit_code"] = 1
    elif problem == "timeout":
        measurement["timed_out"] = True
    elif problem in {"missing_gene", "duplicate_gene"}:
        (output / "orthohmm_orthogroups.txt").write_text("OG1: g1\n" if problem == "missing_gene" else "OG1: g1 g2 g2\n")
    elif problem in {"pair_species", "duplicate_pair"}:
        p = output / "orthohmm_phylogeny/orthohmm_pairwise_orthologs.tsv"
        p.write_text(p.read_text().replace("g1\ta", "g1\tb") if problem == "pair_species" else p.read_text() + "g1\ta\tg2\tb\n")
    else:
        p = Path(run["configuration"]["metrics"])
        d = json.loads(p.read_text())
        if problem == "metrics_command":
            d["command"] = ["other"]
        else:
            d["counts"]["phylogeny_ortholog_pairs"] = 2
        p.write_text(json.dumps(d))
    with pytest.raises(ValueError):
        validate(run, measurement)


@pytest.mark.parametrize("problem", ["version", "incomplete", "copy", "nonfinite", "checkpoint", "orientation"])
def test_orthofinder_corruption_rejected_with_fa_inputs(tmp_path, problem):
    run, measurement = fixture(tmp_path, "orthofinder_full")
    output = Path(run["configuration"]["output"])
    if problem in {"version", "incomplete"}:
        p = output / "Log.txt"
        p.write_text(p.read_text().replace("3.1.5", "2.5.5") if problem == "version" else p.read_text().replace("2026 : OrthoFinder run completed\n", ""))
    elif problem == "copy":
        (output / "input/a.fa").write_text(">g1\nBBB\n")
    elif problem == "nonfinite":
        (output / "OrthoFinder_graph.txt").write_text("0 1:nan $\n")
    elif problem == "checkpoint":
        (output / "clusters_OrthoFinder_I1.2.txt_id_pairs.txt").write_text("begin\n0 0_0 $\n)\n")
    else:
        (output / "Orthologues/Orthologues_b/b__v__a.tsv").write_text("Orthogroup\tb\ta\n")
    with pytest.raises(ValueError):
        validate(run, measurement)
