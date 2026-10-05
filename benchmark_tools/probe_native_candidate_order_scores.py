"""Prespecified fixed-seed hit-order/score factorial for scored native candidates."""

import argparse
import json
import os
from pathlib import Path
import pickle
import platform
import sys
import time

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmark_tools.diagnose_candidate_trace_variation import load_engine
from benchmark_tools.link_factorial_scaling_resources import compare, partition
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_ob_candidate_order_scores import LABELS, factorial_arrays


def save(path, value):
    with Path(path).open("x") as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")


def runtime():
    return {"python": platform.python_version(), "numpy": np.__version__,
            "executable": sys.executable, "python_binary": record(sys.executable),
            "numpy_init": record(np.__file__),
            "numpy_extension": record(np._core._multiarray_umath.__file__)}


def prepare(preparation_ref, score_ref, output):
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(output)
    check(preparation_ref)
    check(score_ref)
    prep = json.loads(Path(preparation_ref["path"]).read_text())
    score = json.loads(Path(score_ref["path"]).read_text())
    if (score.get("schema") != "native_factorial_orthobench_score_v1"
            or score.get("status") != "terminal_native_orthobench_scored"
            or score.get("cell") not in ("p0_c1_r0", "p1_c1_r0")
            or score.get("native_outputs_validated") is not True
            or score.get("prediction_format") != "space_separated_groups"):
        raise ValueError("Require retained terminal-scored candidate-only native cell")
    arm = prep["candidate_arms"][score["cell"][:5]]
    if arm["candidate_partition"] != score["original_prediction"]:
        raise ValueError("Historical partition identity mismatch")
    engine = next(ref for ref in prep["core_sources"] if ref["path"].endswith("/refinement.py"))
    checkpoint = Path(score["prediction"]["path"]).parent / "high_sensitivity_checkpoint"
    checkpoint_refs = [ref for ref in score["evidence"] if Path(ref["path"]).parent == checkpoint]
    if {Path(ref["path"]).name for ref in checkpoint_refs} != {
            "manifest.json", "gene_names.txt", "gene_to_species.npy", "hit_queries.npy", "hit_targets.npy", "hit_scores.npy"}:
        raise ValueError("Require complete native checkpoint bindings")
    sources = [record(Path(__file__).with_name(name)) for name in (
        "probe_native_candidate_order_scores.py", "probe_ob_candidate_order_scores.py",
        "probe_installed_ob_graph.py", "trace_ob_initial_edges.py",
        "diagnose_candidate_trace_variation.py", "link_factorial_scaling_resources.py",
        "prepare_ob_candidate_neighborhood.py", "score_ygob_groups.py")]
    refs = [preparation_ref, score_ref, arm["seed_partition"], prep["cache"], engine,
            arm["candidate_partition"], score["prediction"], *checkpoint_refs, *sources]
    for ref in refs:
        check(ref)
    plan = {"schema": "native_candidate_fixed_seed_factorial_plan_v1", "job_id": score["job_id"],
            "cell": score["cell"], "output": str(output), "labels": list(LABELS),
            "parameters": arm["expansion"]["parameters"], "engine": engine,
            "cache": prep["cache"], "seed": arm["seed_partition"],
            "original_prediction": arm["candidate_partition"], "native_prediction": score["prediction"],
            "checkpoint": str(checkpoint), "checked_records": refs, "runtime": runtime(),
            "attempts": 1, "accuracy_scoring": False, "phylogeny": False,
            "threads": {k: "1" for k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")},
            "limitations": [
                "Fixed original seed boundary; neither new seed inference nor a full end-to-end causal attribution.",
                "Five arms in declared order, one attempt each, exact values without rounding; all failures retained.",
                "Historical and full-fresh controls must match their respective complete saved partitions.",
                "Only candidate expansion is rerun; no benchmark score, phylogeny or comparative timing claim.",
                "Recorded Python/NumPy/source bindings are not a complete dynamic runtime closure."]}
    output.mkdir(parents=True, exist_ok=False)
    save(output / "plan.json", plan)
    return record(output / "plan.json")


def indexed_seed(path, names):
    lookup = {name: i for i, name in enumerate(names)}
    groups = [line.split() for line in Path(path).read_text().splitlines() if line.strip()]
    flat = [g for group in groups for g in group]
    if len(lookup) != len(names) or len(flat) != len(set(flat)) or set(flat) != set(names):
        raise ValueError("Seed does not exactly partition indexed genes")
    return [[lookup[g] for g in group] for group in groups]


def validate_species(old, fresh):
    if old.ndim != 1 or fresh.shape != old.shape or old.dtype.kind not in "iu" or fresh.dtype.kind not in "iu":
        raise ValueError("Malformed species arrays")
    classes = set(zip(old.tolist(), fresh.tolist()))
    if len(classes) != len(set(old.tolist())) or len(classes) != len(set(fresh.tolist())):
        raise ValueError("Species memberships differ")


def load_trusted_hits(cache_ref):
    check(cache_ref)
    with Path(cache_ref["path"]).open("rb") as stream:
        payload = pickle.load(stream)
    if (not isinstance(payload, dict) or not {"all_gene_ids", "gene_to_species", "all_hits"} <= payload.keys()
            or not isinstance(payload["all_hits"], dict) or not isinstance(payload["gene_to_species"], dict)):
        raise ValueError("Malformed retained hit cache")
    names = sorted(payload["all_gene_ids"])
    if not names or len(names) != len(set(names)) or any(not isinstance(g, str) or not g for g in names):
        raise ValueError("Invalid retained gene names")
    lookup = {g: i for i, g in enumerate(names)}
    labels = sorted({str(payload["gene_to_species"][g]) for g in names})
    label_ids = {label: i for i, label in enumerate(labels)}
    species = np.fromiter((label_ids[str(payload["gene_to_species"][g])] for g in names), dtype=np.int32, count=len(names))
    hits = payload["all_hits"]
    q = np.fromiter((lookup[a] for a, b in hits), dtype=np.int32, count=len(hits))
    t = np.fromiter((lookup[b] for a, b in hits), dtype=np.int32, count=len(hits))
    s = np.fromiter(hits.values(), dtype=np.float64, count=len(hits))
    return names, species, q, t, s


def run(plan_ref):
    check(plan_ref)
    plan = json.loads(Path(plan_ref["path"]).read_text())
    if (plan.get("schema") != "native_candidate_fixed_seed_factorial_plan_v1"
            or plan.get("labels") != list(LABELS) or plan.get("attempts") != 1
            or plan.get("accuracy_scoring") is not False or plan.get("phylogeny") is not False):
        raise ValueError("Unexpected factorial protocol")
    output = Path(plan["output"])
    if any((output / name).exists() for name in ("started.json", "report.json", "failure.json")):
        raise FileExistsError("Retained attempt; do not retry or resume")
    if runtime() != plan["runtime"] or any(os.environ.get(k) != v for k, v in plan["threads"].items()):
        raise ValueError("Runtime or thread settings changed")
    for ref in plan["checked_records"]:
        check(ref)
    save(output / "started.json", {"plan": plan_ref, "runtime": runtime(),
        "argv": sys.orig_argv, "affinity": sorted(os.sched_getaffinity(0)), "pid": os.getpid()})
    started = time.monotonic()
    rows, partitions = [], {}
    try:
        # This trusted local pickle is loaded only after its retained checksum check.
        names, species, q, t, s = load_trusted_hits(plan["cache"])
        checkpoint = Path(plan["checkpoint"])
        fresh_names = (checkpoint / "gene_names.txt").read_text().splitlines()
        if names != fresh_names:
            raise ValueError("Historical/fresh gene indexing differs")
        fresh_species, nq, nt, ns = [np.load(checkpoint / filename, mmap_mode="r", allow_pickle=False)
            for filename in ("gene_to_species.npy", "hit_queries.npy", "hit_targets.npy", "hit_scores.npy")]
        validate_species(species, fresh_species)
        manifest = json.loads((checkpoint / "manifest.json").read_text())
        if (manifest.get("schema_version") != 1 or manifest.get("complete") is not True
                or manifest.get("genes") != len(names) or manifest.get("hits") != len(ns)):
            raise ValueError("Checkpoint manifest counts differ")
        arrays, self_hits = factorial_arrays(len(names), (q, t, s), (nq, nt, ns))
        seed = indexed_seed(plan["seed"]["path"], names)
        engine = load_engine(Path(plan["engine"]["path"]))
        for label in LABELS:
            trace = []
            groups, merges, relations, iterations = engine.merge_supported_satellite_candidate_clusters(
                seed, *arrays[label], species, merge_trace=trace, **plan["parameters"])
            path = output / (label + ".txt")
            with path.open("x") as stream:
                for group in groups:
                    stream.write(" ".join(names[g] for g in group) + "\n")
            partitions[label] = partition(path, "space_separated_groups")
            save(output / (label + "_merges.json"), [
                {k: ([names[g] for g in v] if k in ("source_genes", "target_genes") else
                     "positive_infinity" if isinstance(v, float) and v == float("inf") else v)
                 for k, v in row.items()} for row in trace])
            rows.append({"label": label, "partition": record(path), "groups": len(groups),
                         "merges": merges, "relations": relations, "iterations": iterations})
            save(output / (label + "_result.json"), rows[-1])
        controls = {"historical": partitions[LABELS[0]] == partition(Path(plan["original_prediction"]["path"]), "space_separated_groups"),
                    "fresh_full": partitions[LABELS[-1]] == partition(Path(plan["native_prediction"]["path"]), "space_separated_groups")}
        comparisons = {a + "__vs__" + b: compare(partitions[a], partitions[b], set())
                       for i, a in enumerate(LABELS) for b in LABELS[i + 1:]}
        for ref in plan["checked_records"]:
            check(ref)
        save(output / "report.json", {"schema": "native_candidate_fixed_seed_factorial_result_v1",
            "plan": plan_ref, "rows": rows, "controls": controls, "comparisons": comparisons,
            "genes": len(names), "nonself_hits": len(q), "removed_self_hits": self_hits,
            "status": "controls_reproduced" if all(controls.values()) else "control_mismatch",
            "diagnostic_wall_seconds": time.monotonic() - started, "accuracy_scored": False,
            "frozen_method_modified": False, "limitations": plan["limitations"]})
        if not all(controls.values()):
            raise ValueError("Whole-partition historical/fresh control mismatch; retain all arms")
    except BaseException as error:
        save(output / "failure.json", {"error_type": type(error).__name__, "error": str(error),
             "completed_labels": [row["label"] for row in rows], "retry": False})
        raise
    return record(output / "report.json")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preparation", type=Path)
    parser.add_argument("--preparation-sha256")
    parser.add_argument("--score", type=Path)
    parser.add_argument("--score-sha256")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--run", type=Path)
    parser.add_argument("--plan-sha256")
    args = parser.parse_args()
    if args.run:
        ref = record(args.run)
        if ref["sha256"] != args.plan_sha256:
            raise ValueError("Changed plan checksum")
        result = run(ref)
    else:
        if not all((args.preparation, args.score, args.output)):
            parser.error("Require --run or preparation, score and output")
        refs = [record(args.preparation), record(args.score)]
        if [ref["sha256"] for ref in refs] != [args.preparation_sha256, args.score_sha256]:
            raise ValueError("Changed preparation/score checksum")
        result = prepare(*refs, args.output)
    print(json.dumps(result, sort_keys=True))
