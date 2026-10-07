"""Infer label-independent distances from the retained SwissTrees alignments."""

import argparse
import csv
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import re
import signal
import statistics
import subprocess
import sys
import time

import Bio
from Bio import Phylo, SeqIO
from Bio.Phylo.NewickIO import NewickError
import psutil

INPUTS = {
    "alignments": ("corrected_swiss_identity_prepared_22102.json",
                   "7bec0f40c8227f230bbf06cb52622947e256f2edfc2c88b991393fa570603e44"),
    "admission": ("corrected_swiss_identity_admission_22102.json",
                  "cb04162af62fbd58fcfe8f02cbb78bc53bcf49ac20a487aae89911bc3a4e7b2d"),
    "protocol": ("SWISS_MODEL_DIVERGENCE_PROTOCOL_20261007.md",
                 "a596371c443471585b9a277cbdd1466bb3f9870200f3187d20286461d3699bfd"),
}
BINARY = Path("/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/"
              "iqtree-3.0.1-Linux-intel/bin/iqtree3")
BINARY_SHA = "40424ccdb1d79c304641f910cb6c172ebb670214e50958352018e3ff9906ab8f"
MODEL, SEED, TIMEOUT = "WAG+G4", 20261007, 600
UNIT = "model_estimated_expected_amino_acid_substitutions_per_site"
PAIR_FIELDS = ("family", "gene_a", "gene_b", "distance")
SOURCE_FILES = (
    "benchmark_tools/prepare_swiss_model_divergence.py",
    "benchmark_tools/swiss_model_divergence_batch_20261007.sh",
    "benchmark_tools/results/SWISS_MODEL_DIVERGENCE_PROTOCOL_20261007.md",
)
SCOPE = dict(prediction_statistics_evaluated=False, independent_confirmation=False,
             publication_ready=False, scientific_timings_admitted=False,
             new_accuracy_or_resource_admission=False, new_uncertainty=False,
             new_bootstrap_draws=0, raw_scorer_repeated=False, alignment_repeated=False)
LIMITATIONS = [
    "Development-exposed families; descriptive fixed-model distances, not independent confirmation.",
    "Alignment/model/sampling/paralogy/domain/composition effects may affect these distances.",
    "No elapsed biological time, known ancestral history or causal error explanation.",
    "Fixed WAG+G4 model fit and heuristic topology-search uncertainty are not validated.",
    "No support test, model comparison, outcome-selected cutoff, confidence interval or tuning.",
    "All tip pairs, including paralogs, within-species pairs and zero distances, are included.",
    "Only direct selected input/source bindings checked; historical admission is inherited.",
    "Construction cost is shared-host postprocessing, not OrthoHMM inference or isolated speed.",
]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def record(path):
    path = Path(path).resolve(strict=True)
    return dict(path=str(path), bytes=path.stat().st_size,
                sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def check(ref):
    require(record(ref["path"]) == ref, "Changed binding: " + ref["path"])


def save(path, data):
    with Path(path).open("x", encoding="ascii") as stream:
        json.dump(data, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")


def load_inputs(repo):
    base = Path(repo) / "benchmark_tools/results"
    refs, docs = {}, {}
    for key, (name, sha) in INPUTS.items():
        refs[key] = record(base / name)
        require(refs[key]["sha256"] == sha, "Changed prospective input: " + key)
        if name.endswith(".json"):
            docs[key] = json.loads(Path(refs[key]["path"]).read_text())
    inventory, admission = docs["alignments"], docs["admission"]
    require(inventory["status"] == "corrected_swiss_alignments_prepared_unscored"
            and inventory["failed_families"] == []
            and inventory["prediction_statistics_evaluated"] is False
            and admission["status"] == "corrected_swiss_identity_features_verified"
            and admission["prediction_statistics_evaluated"] is False,
            "Incorrect label-independent alignment admission")
    memberships = admission["family_memberships"]
    genes = [g for f in sorted(memberships) for g in memberships[f]]
    require(len(memberships) == 18 and len(genes) == len(set(genes)) == 563
            and all(v == sorted(set(v)) for v in memberships.values()),
            "Changed canonical membership universe")
    runs = inventory["runs"]
    require([r["refog"] for r in runs] == sorted(memberships), "Incomplete selected families")
    for run in runs:
        ref = run["alignment"]
        require(run["status"] == "alignment_validated" and run["exit_code"] == 0
                and ref in admission["records"], "Unadmitted selected alignment")
        check(ref)
        alignment_members(ref["path"], memberships[run["refog"]], run["columns"])
        require(run["genes"] == len(memberships[run["refog"]]), "Incorrect family size")
    return refs, runs, memberships


def alignment_members(path, expected, columns):
    sequences = list(SeqIO.parse(path, "fasta"))
    names = [r.id for r in sequences]
    require(sorted(names) == expected and len(names) == len(set(names)), "Alignment members differ")
    require(type(columns) is int and columns > 0
            and all(len(r.seq) == columns and set(str(r.seq).upper()) <=
                    set("ACDEFGHIKLMNPQRSTVWYX-")
                    and set(str(r.seq).upper()) - {"-", "X"} for r in sequences),
            "Invalid aligned amino-acid data")
    return names


def command(binary, alignment, prefix):
    return [str(binary), "-s", str(alignment), "--seqtype", "AA", "-m", MODEL,
            "--seed", str(SEED), "-T", "1", "--mem", "4G", "-keep-ident",
            "--prefix", str(prefix)]


def tree_features(path, expected):
    try:
        tree = Phylo.read(path, "newick")
    except NewickError as error:
        raise ValueError("Malformed Newick tree") from error
    tips = tree.get_terminals()
    require(sorted(t.name for t in tips) == expected
            and len(tips) == len(set(t.name for t in tips)), "Tree members differ")
    for node in tree.find_clades():
        length = node.branch_length
        require((node is tree.root and length is None)
                or type(length) in (int, float) and math.isfinite(length) and length >= 0,
                "Invalid or missing branch length")
    pairs = []
    for a, b in itertools.combinations(expected, 2):
        distance = tree.distance(a, b)
        require(math.isfinite(distance) and distance >= 0, "Invalid patristic distance")
        pairs.append(dict(gene_a=a, gene_b=b, distance=distance))
    require(bool(pairs), "No unordered pairs")
    values = [p["distance"] for p in pairs]
    features = dict(pairs=len(values), median_pair_distance=statistics.median(values),
                    mean_pair_distance=math.fsum(values) / len(values),
                    minimum_pair_distance=min(values), maximum_pair_distance=max(values),
                    tree_length=math.fsum(n.branch_length for n in tree.find_clades()
                                         if n is not tree.root), members=expected, unit=UNIT)
    return features, pairs


def bins(features):
    require(bool(features) and all(math.isfinite(v["median_pair_distance"])
                                  and v["median_pair_distance"] >= 0 for v in features.values()),
            "Invalid bin features")
    cutoff = statistics.median(v["median_pair_distance"] for v in features.values())
    return cutoff, {
        "all": sorted(features),
        "lower_or_equal_median": sorted(f for f, v in features.items()
                                        if v["median_pair_distance"] <= cutoff),
        "higher_than_median": sorted(f for f, v in features.items()
                                    if v["median_pair_distance"] > cutoff),
    }


def run_family(argv, stdout, stderr, timeout=TIMEOUT):
    with Path(stdout).open("xb") as out, Path(stderr).open("xb") as err:
        process = subprocess.Popen(argv, stdout=out, stderr=err, start_new_session=True)
        try:
            return process.wait(timeout=timeout), False
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
            return process.returncode, True


def committed_sources(repo, commit):
    require(re.fullmatch(r"[a-f0-9]{40}", commit) is not None, "Require explicit full source commit")
    sources = []
    for name in SOURCE_FILES:
        path = Path(repo) / name
        require(path.read_bytes() == subprocess.check_output(["git", "show", commit + ":" + name],
                                                             cwd=repo), "Uncommitted or changed source")
        sources.append(record(path))
    return sources


def execute(repo, output, source_commit, job_id):
    repo, output = Path(repo).resolve(), Path(output).resolve()
    require(not output.exists() and not output.is_symlink(), "Existing attempt; never resume or retry")
    require(output.is_relative_to(repo / "benchmarks/results"), "Output outside new analysis namespace")
    refs, runs, memberships = load_inputs(repo)
    sources = committed_sources(repo, source_commit)
    binary = record(BINARY)
    require(binary["sha256"] == BINARY_SHA and binary["bytes"] == 11333032, "Changed trusted binary")
    require(job_id > 0 and os.environ.get("SLURM_JOB_ID") == str(job_id)
            and os.environ.get("SLURM_CPUS_PER_TASK") == "2"
            and os.environ.get("SLURM_MEM_PER_NODE") == "8192", "Invalid selected Slurm allocation")
    require(all(not os.environ.get(k) for k in ("PYTHONPATH", "PYTHONHOME", "PYTHONUSERBASE",
                "LD_PRELOAD", "LD_LIBRARY_PATH", "LD_AUDIT")), "Unsanitized runtime")
    require(all(os.environ.get(k) == "1" for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                "MKL_NUM_THREADS")), "Unbounded library threading")
    memory = psutil.virtual_memory()
    require(memory.available >= 8 * 1024 ** 3, "Unsafe available RAM; defer only this launch")
    version = subprocess.check_output([str(BINARY), "--version"], text=True)
    require(version.startswith("IQ-TREE version 3.0.1 "), "Unexpected binary runtime")
    output.mkdir()
    start = time.time_ns()
    preflight = dict(inputs=refs, source_commit=source_commit, sources=sources, binary=binary,
                     binary_version=version, python=sys.version, python_executable=record(sys.executable),
                     biopython=Bio.__version__, psutil=psutil.__version__, job_id=job_id,
                     cpus_per_task=2, memory_mib=8192, available_memory_bytes=memory.available,
                     cpu_affinity=psutil.Process().cpu_affinity(), host_load_average=os.getloadavg(),
                     start_ns=start, memberships=memberships, model=MODEL, seed=SEED,
                     family_timeout_seconds=TIMEOUT, unit=UNIT, **SCOPE)
    save(output / "preflight.json", preflight)
    results, features = [], {}
    pairs_path = output / "pairs.tsv"
    with pairs_path.open("x", encoding="ascii", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=PAIR_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for selected in runs:
            family = selected["refog"]
            directory = output / family
            directory.mkdir()
            argv = command(BINARY, selected["alignment"]["path"], directory / "inference")
            result = dict(family=family, alignment=selected["alignment"], columns=selected["columns"],
                          command=argv, start_ns=time.time_ns(), status="failed", features=None)
            try:
                code, timed_out = run_family(argv, directory / "stdout.txt", directory / "stderr.txt")
                result.update(exit_code=code, timed_out=timed_out)
                require(code == 0 and not timed_out, "IQ-TREE failed or timed out")
                report = directory / "inference.iqtree"
                require(re.search(r"^Model of substitution: WAG\+G4\s*$", report.read_text(), re.M)
                        is not None, "Unexpected reported substitution model")
                feature, pairs = tree_features(directory / "inference.treefile", memberships[family])
                result.update(status="feature_constructed", features=feature)
                features[family] = feature
                for pair in pairs:
                    writer.writerow(dict(family=family, **pair))
                stream.flush()
            except (ValueError, OSError) as error:
                result["error"] = str(error)
            result["finish_ns"] = time.time_ns()
            result["outputs"] = [record(p) for p in sorted(directory.iterdir()) if p.is_file()]
            save(directory / "result.json", result)
            results.append(result)
    failed = [r["family"] for r in results if r["status"] != "feature_constructed"]
    cutoff, strata = bins(features) if not failed else (None, None)
    report = dict(schema="swiss_model_divergence_features_v1",
                  status="features_constructed_unverified" if not failed else "selected_family_failure",
                  preflight=record(output / "preflight.json"), source=record(__file__),
                  source_commit=source_commit, inputs=refs, binary=binary, memberships=memberships,
                  model=MODEL, seed=SEED, unit=UNIT, runs=results, failed_families=failed,
                  median_family_distance=cutoff, strata=strata, pairs=record(pairs_path),
                  start_ns=start, finish_ns=time.time_ns(), limitations=LIMITATIONS, **SCOPE)
    save(output / "report.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--job-id", type=int, required=True)
    args = parser.parse_args()
    result = execute(args.repo, args.output, args.source_commit, args.job_id)
    print(json.dumps(dict(status=result["status"], failed_families=result["failed_families"])))
    sys.exit(bool(result["failed_families"]))
