"""Independent numeric, graph-edge and partition readback before orthology scoring."""

import argparse
import csv
import json
import math
from pathlib import Path
import subprocess

import numpy as np
from Bio import SeqIO

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_matched_graph_cell import MANIFEST_SHA


def partition(path, names):
    groups = [line.rstrip("\n").split("\t") for line in path.read_text().splitlines()]
    flat = [g for group in groups for g in group]
    if not groups or any(not g for g in flat) or len(flat) != len(set(flat)) or set(flat) != set(names):
        raise ValueError("Partition does not cover each input gene exactly once")
    return groups


def raw_numeric(cell, dataset):
    genes = {}
    for item in dataset["inputs"]:
        check(item)
        for seq in SeqIO.parse(item["path"], "fasta"):
            if seq.id in genes:
                raise ValueError("Duplicate input gene")
            genes[seq.id] = (Path(item["path"]).name, len(seq))
    names = sorted(genes)
    ids = {g: i for i, g in enumerate(names)}
    species = {s: i for i, s in enumerate(sorted({v[0] for v in genes.values()}))}
    hits = {}
    for item in cell["sources"]:
        check(item)
        with Path(item["path"]).open() as stream:
            rows = csv.reader(stream, delimiter="\t")
            if cell["arm"] == "hmm":
                if next(rows) != ["query_species", "target_species", "query_id", "target_id", "score", "evalue"]:
                    raise ValueError("HMM header mismatch")
            for row in rows:
                if cell["arm"] == "hmm":
                    qs, ts, q, t, score, e = row
                    if (qs, ts) != (genes[q][0], genes[t][0]) or not 0 <= float(e) < 1e-4:
                        raise ValueError("Invalid native HMM hit")
                    value = float(score)
                else:
                    q, t, qlen, tlen, raw, bits, e = row
                    if (int(qlen), int(tlen)) != (genes[q][1], genes[t][1]):
                        raise ValueError("DIAMOND length mismatch")
                    if float(e) > 1e-40:
                        continue
                    value = float(raw) / math.sqrt(int(qlen) * int(tlen))
                key = ids[q], ids[t]
                if key in hits or not math.isfinite(value) or value <= 0:
                    raise ValueError("Invalid duplicate/nonpositive graph hit")
                hits[key] = value
    keys = sorted(hits)
    return dict(gene_names=names, gene_to_species=[species[genes[g][0]] for g in names],
                hit_queries=[q for q, t in keys], hit_targets=[t for q, t in keys],
                hit_scores=[hits[k] for k in keys])


def checkpoint(directory, numeric):
    manifest = json.loads((directory / "manifest.json").read_text())
    expected = {"gene_names.txt", *(k + ".npy" for k in numeric if k != "gene_names")}
    if set(manifest["files"]) != expected or {p.name for p in directory.iterdir()} != expected | {"manifest.json"}:
        raise ValueError("Checkpoint inventory mismatch")
    for name, expected_record in manifest["files"].items():
        observed = record(directory / name)
        if any(observed[k] != expected_record[k] for k in ("bytes", "sha256")):
            raise ValueError("Checkpoint checksum mismatch")
    if (directory / "gene_names.txt").read_text().splitlines() != numeric["gene_names"]:
        raise ValueError("Checkpoint gene order differs")
    for key in numeric.keys() - {"gene_names"}:
        array = np.load(directory / (key + ".npy"), allow_pickle=False)
        expected_array = np.asarray(numeric[key], dtype=np.float64 if key == "hit_scores" else np.int32)
        if array.dtype != expected_array.dtype or array.shape != expected_array.shape or not np.array_equal(array, expected_array):
            raise ValueError("Checkpoint numeric array differs")
    if not manifest["complete"] or manifest["genes"] != len(numeric["gene_names"]) or manifest["hits"] != len(numeric["hit_scores"]):
        raise ValueError("Checkpoint counts differ")


def reference_edges(numeric, initial):
    species = numeric["gene_to_species"]
    hits = list(zip(numeric["hit_queries"], numeric["hit_targets"], numeric["hit_scores"]))
    best = {}
    for q, t, s in hits:
        key = q, species[t]
        if q != t and (key not in best or s > best[key][1]):
            best[key] = t, s
    thresholds = [math.inf] * len(species)
    for (q, _), (t, s) in best.items():
        reverse = best.get((t, species[q]))
        if reverse and reverse[0] == q:
            value = (s + reverse[1]) / 2
            thresholds[q] = min(thresholds[q], value)
            thresholds[t] = min(thresholds[t], value)
    rbnh = {}
    def add(edges, q, t, value):
        key = tuple(sorted((q, t)))
        edges[key] = max(edges.get(key, -math.inf), value)
    for q, t, s in hits:
        if q != t and s >= min(thresholds[q], thresholds[t]):
            add(rbnh, q, t, s)
    owners = {g: i for i, group in enumerate(initial) if len(group) > 1 for g in group}
    names = numeric["gene_names"]
    chosen = {}
    eligible = [(q, t, s) for q, t, s in hits if names[q] not in owners and names[t] in owners]
    for q, t, s in eligible:
        if q not in chosen or s > chosen[q][1]:
            chosen[q] = owners[names[t]], s
    multipass = dict(rbnh)
    for q, t, s in eligible:
        if owners[names[t]] == chosen[q][0]:
            add(multipass, q, t, s)
    return rbnh, multipass


def verify_edges(path, expected):
    with np.load(path, allow_pickle=False) as data:
        if set(data.files) != {"sources", "targets", "weights"}:
            raise ValueError("Unexpected edge arrays")
        q, t, w = [data[k] for k in ("sources", "targets", "weights")]
    if q.ndim != 1 or q.shape != t.shape or q.shape != w.shape:
        raise ValueError("Graph array shape mismatch")
    if q.dtype != np.int32 or t.dtype != np.int32 or w.dtype != np.float64:
        raise ValueError("Graph dtypes differ")
    observed = {(int(a), int(b)): float(s) for a, b, s in zip(q, t, w)}
    if len(observed) != len(q) or observed != expected:
        raise ValueError("Graph differs from independent edge derivation")


def execution_bindings(path, execution, native, python, installed):
    expected_env = dict(PATH="/usr/bin:/bin", HOME=str(path), LC_ALL="C",
                        OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                        MKL_NUM_THREADS="1", PYTHONHASHSEED="0")
    if execution["environment"] != expected_env:
        raise ValueError("Recorded execution environment differs")
    expected_records = dict(private_numeric=record(path / "numeric.json"),
                            native_receipt=record(path / "graph/receipt.json"))
    if any(execution[key] != item for key, item in expected_records.items()):
        raise ValueError("Execution artifact binding differs")
    if execution["logs"] != [record(path / name) for name in ("graph.log", "graph.time.txt")]:
        raise ValueError("Execution log inventory differs")
    if (native["numeric"] != expected_records["private_numeric"]
            or native["status"] != "native_graph_completed_pending_independent_readback"
            or native["executable"] != str(python) or native["prefix"] != str(python.parent.parent)):
        raise ValueError("Native input/runtime binding differs")
    expected_outputs = [record(path / "graph" / name) for name in
                        ("rbnh_edges.npz", "multipass_edges.npz", "initial.tsv", "multipass.tsv", "final.tsv")]
    expected_checkpoint = record(path / "graph/orthohmm_working_res/high_sensitivity_checkpoint/manifest.json")
    if native["outputs"] != expected_outputs or native["checkpoint_manifest"] != expected_checkpoint:
        raise ValueError("Native output inventory differs")
    # Membership alone permits empty or duplicate lists; require the full module set.
    modules = native["modules"]
    expected_names = {"accuracy.py", "externals.py", "refinement.py"}
    prefix = python.parent.parent.resolve()
    if (len(modules) != len(expected_names)
            or {Path(r["path"]).name for r in modules} != expected_names
            or any(r not in installed["checked_records"] or
                   not Path(r["path"]).resolve().is_relative_to(prefix) for r in modules)):
        raise ValueError("Native module inventory differs")


def audit(submission_path, output):
    if output.exists():
        raise FileExistsError(output)
    submission = json.loads(submission_path.read_text())
    for item in [submission["manifest"], submission["installation"], *submission["sources"]]:
        check(item)
    if submission["manifest"]["sha256"] != MANIFEST_SHA:
        raise ValueError("Unexpected numeric manifest")
    manifest = json.loads(Path(submission["manifest"]["path"]).read_text())
    installed = json.loads(Path(submission["installation"]["path"]).read_text())
    for item in manifest["checked_records"]:
        check(item)
    search = json.loads(Path(manifest["search_result"]["path"]).read_text())
    check(search["plan"])
    plan = json.loads(Path(search["plan"]["path"]).read_text())
    accounting = subprocess.check_output(["sacct", "-j", submission["job_id"], "--format=JobID,State,ExitCode", "-n", "-P"], text=True)
    states = {r[0]: r[1:3] for line in accounting.splitlines() if len(r := line.split("|")) >= 3}
    if any(states.get(f'{submission["job_id"]}_{i}') != ["COMPLETED", "0:0"] for i in range(70)):
        raise ValueError("Require all 70 successful terminal jobs")
    python = next(Path(r["path"]).parent / "venv_clean/bin/python" for r in installed["checked_records"] if Path(r["path"]).name == "staging.json")
    worker = Path(submission["executor"]) / "benchmark_tools/matched_graph_worker.py"
    rows = []
    for index, cell in enumerate(manifest["cells"]):
        path = Path(submission["output"]) / f"cell_{index}"
        execution = json.loads((path / "execution.json").read_text())
        if (execution["status"] != "native_completed_pending_independent_readback" or execution["cell"] != cell
                or execution["attempt"] != 1 or execution["index"] != index
                or execution["slurm_array_job_id"] != submission["job_id"] or execution["slurm_array_task_id"] != str(index)
                or execution["manifest"] != submission["manifest"] or execution["installation"] != submission["installation"]):
            raise ValueError("Execution provenance mismatch")
        for item in [*execution["checked_records"], execution["private_numeric"], execution["native_receipt"], *execution["logs"]]:
            check(item)
        if any(r not in execution["checked_records"] for r in installed["checked_records"]):
            raise ValueError("Incomplete installed-runtime checks")
        expected_command = [str(python), "-I", str(worker), "--numeric", str(path / "numeric.json"), "--output", str(path / "graph")]
        if len(execution["stages"]) != 1 or execution["stages"][0]["command"] != expected_command or execution["stages"][0]["returncode"] != 0:
            raise ValueError("Changed inference command")
        native = json.loads((path / "graph/receipt.json").read_text())
        execution_bindings(path, execution, native, python, installed)
        if (native["source"] != record(worker) or native["settings"] != manifest["graph_settings"]
                or not native["isolated"] or native["truth_loaded"]
                or any(r not in installed["checked_records"] for r in native["modules"])):
            raise ValueError("Native production identity/settings mismatch")
        for item in [*native["outputs"], native["checkpoint_manifest"]]:
            check(item)
        numeric = json.loads(Path(cell["numeric"]["path"]).read_text())
        if json.loads((path / "numeric.json").read_text()) != numeric:
            raise ValueError("Private numeric input differs")
        dataset = plan["datasets"][cell["dataset_index"]]
        if dataset["split"] != "reporting" or any(dataset[k] != cell[k] for k in ("condition", "seed")):
            raise ValueError("Wrong reporting dataset")
        if raw_numeric(cell, dataset) != numeric:
            raise ValueError("Raw-to-numeric score mapping differs")
        checkpoint(path / "graph/orthohmm_working_res/high_sensitivity_checkpoint", numeric)
        groups = {name: partition(path / "graph" / (name + ".tsv"), numeric["gene_names"]) for name in ("initial", "multipass", "final")}
        for name, expected in zip(("rbnh", "multipass"), reference_edges(numeric, groups["initial"])):
            verify_edges(path / "graph" / (name + "_edges.npz"), expected)
        if native["genes"] != len(numeric["gene_names"]) or native["groups"] != len(groups["final"]) or native["hits"] != len(numeric["hit_scores"]):
            raise ValueError("Native counts differ")
        rows.append(dict(index=index, condition=cell["condition"], seed=cell["seed"], arm=cell["arm"],
                         execution=record(path / "execution.json"), native_receipt=record(path / "graph/receipt.json"),
                         final_partition=record(path / "graph/final.tsv"), truth=dataset["truth"],
                         inputs=dataset["inputs"], genes=len(numeric["gene_names"]), groups=len(groups["final"])))
    result = dict(status="all_native_graphs_verified_pending_orthology_scoring", source=record(__file__),
                  submission=record(submission_path), accounting=accounting, cells=rows,
                  checks=["Raw-to-normalized numeric mapping", "Exact saved checkpoint arrays", "Independent RBNH and singleton edge derivation", "Complete unique-gene partitions", "Pinned settings/source and terminal success", "Exact recorded environment and artifact/runtime bindings"],
                  limitations=["Does not independently rerun Leiden or derive the refinement decisions.", "Environment is checked against the runner receipt, not independently observed inside the historical process.", "No orthology scores calculated; simulations remain development-exposed."],
                  accuracy_evaluated=False, publication_ready=False)
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--submission", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    audit(args.submission.resolve(), args.output.absolute())
