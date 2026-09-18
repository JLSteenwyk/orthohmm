"""Compare serial and pair-parallel native inference on the bundled small sample."""

import argparse
import io
import json
from pathlib import Path
import shutil
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_qfo_corrected_orthomcl_native import write_sources, TOOL, RUNTIME_SHA, HELPERS_SHA
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_qfo_corrected_blast import environment
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.snapshot_runtime_trees import verify
from benchmark_tools.normalize_three_kingdoms_orthogroups import iter_orthomcl
from benchmark_tools.orthomcl_matrix_to_pairwise import load_species, load_index, write_pairs


def partition(path, universe):
    groups, seen = [], set()
    for _, genes in iter_orthomcl(path):
        members = set(genes)
        if len(members) != len(genes) or not members <= universe or members & seen:
            raise ValueError("Unknown/duplicate native sample group membership")
        seen.update(members)
        groups.append(tuple(sorted(members)))
    if not groups:
        raise ValueError("Empty native sample groups")
    return sorted(groups)


def weighted_graph(matrix, index, gg):
    write_pairs(matrix, index, gg, io.StringIO())
    genes, _ = load_index(index, load_species(gg))
    edges, inside = [], False
    for line in matrix.read_text().splitlines():
        text = line.strip()
        if text == "begin":
            inside = True
        elif text == ")":
            inside = False
        elif inside and text:
            fields = text[:-1].split()
            source = genes[int(fields[0])]
            for token in fields[1:]:
                target, score = token.split(":")
                edges.append((source, genes[int(target)], score))
    return sorted(edges)


def probe(root, output):
    if output.exists():
        raise FileExistsError(output)
    results = root / "benchmark_tools/results"
    paths = [results / "qfo_corrected_orthomcl_perl_runtime_20260918.json",
             results / "qfo_corrected_orthomcl_system_helpers_20260918.json"]
    manifests = [read_frozen(path, sha) for path, sha in zip(paths, (RUNTIME_SHA, HELPERS_SHA))]
    for manifest in manifests:
        verify(manifest)
    sample = TOOL / "sample_data"
    checked = [*[record(path) for path in paths],
               *[record(sample / name) for name in ("AtHsSc.bpo", "AtHsSc.gg", "all_orthomcl.out")],
               *[record(Path(__file__).with_name(name)) for name in (
                   "probe_orthomcl_native_inference.py", "prepare_qfo_corrected_orthomcl_native.py",
                   "configure_orthomcl_1_4.py", "parallelize_orthomcl_pairs.py", "run_orthomcl_perl_script.pl",
                   "normalize_three_kingdoms_orthogroups.py", "orthomcl_matrix_to_pairwise.py")]]
    universe = set(load_species(sample / "AtHsSc.gg"))
    reference = partition(sample / "all_orthomcl.out", universe)
    output.mkdir(parents=True, exist_ok=False)
    report = {"status": "running", "checked_records": checked, "arms": {},
              "accuracy_admitted": False, "publication_ready": False}
    try:
        partitions, matrix_bytes, index_bytes, graphs = [], [], [], []
        for name, parallel, workers in (("native_serial", False, 1), ("patched_one", True, 1), ("patched_two", True, 2)):
            directory = output / name
            data = directory / "data"
            data.mkdir(parents=True)
            for old, new in (("AtHsSc.bpo", "all.bpo"), ("AtHsSc.gg", "all.gg")):
                shutil.copyfile(sample / old, data / new)
            sources = write_sources(directory / "tool", data, 2, parallel)
            argv = [str(TOOL / "venv_orthomcl/bin/perl"), "-I" + str(directory / "tool"),
                    str(Path(__file__).with_name("run_orthomcl_perl_script.pl")), str(directory / "tool/orthomcl.pl"),
                    "--mode", "4", "--bpo_file", str(data / "all.bpo"), "--gg_file", str(data / "all.gg")]
            env = {**environment(), "ORTHOMCL_PAIR_WORKERS": str(workers)}
            with (directory / "native.log").open("xb") as log:
                done = subprocess.run(argv, cwd=directory, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=120)
            arm = {"source_preparation": sources, "command": argv, "cwd": str(directory),
                   "environment": env, "exit_code": done.returncode}
            report["arms"][name] = arm
            if done.returncode:
                raise ValueError("Native fixture failed: " + name)
            outputs = list((directory / "tool").glob("*/all_orthomcl.out"))
            if len(outputs) != 1:
                raise ValueError("Require one native fixture result")
            native = outputs[0].parent
            groups = partition(outputs[0], universe)
            partitions.append(groups)
            matrix_bytes.append((native / "tmp/all_ortho.mtx").read_bytes())
            index_bytes.append((native / "tmp/all_ortho.idx").read_bytes())
            graph = weighted_graph(native / "tmp/all_ortho.mtx", native / "tmp/all_ortho.idx", data / "all.gg")
            graphs.append(graph)
            arm.update(groups=len(groups), grouped_proteins=sum(map(len, groups)),
                       bundled_partition_equal=groups == reference,
                       groups_only_in_current=[list(g) for g in groups if g not in reference],
                       groups_only_in_bundled=[list(g) for g in reference if g not in groups],
                       weighted_directed_edges=len(graph),
                       outputs=[record(p) for p in sorted(directory.rglob("*")) if p.is_file()])
        report["partition_equal"] = all(value == partitions[0] for value in partitions)
        report["matrix_bytes_equal"] = all(value == matrix_bytes[0] for value in matrix_bytes)
        report["matrix_index_bytes_equal"] = all(value == index_bytes[0] for value in index_bytes)
        report["weighted_graph_equal"] = all(value == graphs[0] for value in graphs)
        for manifest in manifests:
            verify(manifest)
        for item in checked:
            check(item)
        for arm in report["arms"].values():
            for item in [*arm["source_preparation"]["originals"], *arm["outputs"]]:
                check(item)
        report["status"] = ("native_inference_fixture_parity_verified" if all(report[key] for key in
            ("partition_equal", "weighted_graph_equal", "matrix_index_bytes_equal"))
            else "native_inference_fixture_difference_requires_review")
        report["limitations"] = [
            "Small bundled sample only, not corrected QfO or independent biological validation.",
            "Serial/one-worker/two-worker checks do not prove behavior for all datasets or 64-worker schedules.",
            "Bundled expected groups are tool regression evidence, not orthology ground truth.",
            "Raw matrix order differences are retained; graph comparison checks every directed endpoint and printed weight after native matrix validation.",
            "Native intermediate and cache outputs are retained in fresh per-arm directories; no implicit resume."]
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        with (output / "report.json").open("x") as stream:
            json.dump(report, stream, indent=2, sort_keys=True)
            stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = probe(args.root.resolve(), args.output.resolve())
    print(json.dumps({key: report.get(key) for key in ("status", "partition_equal", "matrix_bytes_equal", "matrix_index_bytes_equal")}))
