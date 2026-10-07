"""Independent edge-split readback of newly inferred SwissTrees distances."""

import argparse
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import statistics
import subprocess

from Bio import Phylo, SeqIO
from Bio.Phylo.NewickIO import NewickError

UNIT = "model_estimated_expected_amino_acid_substitutions_per_site"
PINS = {
    "alignments": "7bec0f40c8227f230bbf06cb52622947e256f2edfc2c88b991393fa570603e44",
    "admission": "cb04162af62fbd58fcfe8f02cbb78bc53bcf49ac20a487aae89911bc3a4e7b2d",
    "protocol": "a596371c443471585b9a277cbdd1466bb3f9870200f3187d20286461d3699bfd",
}
BINARY_SHA = "40424ccdb1d79c304641f910cb6c172ebb670214e50958352018e3ff9906ab8f"
FLAGS = ("prediction_statistics_evaluated", "independent_confirmation", "publication_ready",
         "scientific_timings_admitted", "new_accuracy_or_resource_admission", "new_uncertainty",
         "raw_scorer_repeated", "alignment_repeated")


def demand(condition, message):
    if not condition:
        raise ValueError(message)


def fingerprint(path):
    path = Path(path).resolve(strict=True)
    return dict(path=str(path), bytes=path.stat().st_size,
                sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def close(actual, expected):
    demand(type(actual) in (int, float) and math.isfinite(actual)
           and math.isclose(actual, expected, rel_tol=0, abs_tol=1e-10), "Numerical readback differs")


def edge_distances(path, members):
    try:
        tree = Phylo.read(path, "newick")
    except NewickError as error:
        raise ValueError("Malformed Newick tree") from error
    leaves = [tip.name for tip in tree.get_terminals()]
    demand(sorted(leaves) == members and len(leaves) == len(set(leaves)), "Invalid full tip set")
    edges = []
    for node in tree.find_clades():
        length = node.branch_length
        demand(node is tree.root and length is None
               or type(length) in (int, float) and math.isfinite(length) and length >= 0,
               "Invalid tree edge length")
        if node is not tree.root:
            edges.append((length, {tip.name for tip in node.get_terminals()}))
    # An edge belongs to a tip-to-tip path iff its split separates the tips.
    pairs = {(a, b): math.fsum(length for length, tips in edges if (a in tips) != (b in tips))
             for a, b in itertools.combinations(members, 2)}
    demand(bool(pairs), "Empty pair set")
    values = sorted(pairs.values())
    middle = len(values) // 2
    median = values[middle] if len(values) % 2 else (values[middle - 1] + values[middle]) / 2
    descriptors = dict(pairs=len(values), median_pair_distance=median,
                       mean_pair_distance=math.fsum(values) / len(values),
                       minimum_pair_distance=values[0], maximum_pair_distance=values[-1],
                       tree_length=math.fsum(length for length, _ in edges), members=members, unit=UNIT)
    return descriptors, pairs


def compare_pairs(path, expected):
    with Path(path).open(newline="", encoding="ascii") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        demand(reader.fieldnames == ["family", "gene_a", "gene_b", "distance"], "Wrong pair TSV fields")
        rows = list(reader)
    keys = [(r["family"], r["gene_a"], r["gene_b"]) for r in rows]
    demand(keys == sorted(expected) and len(keys) == len(set(keys)), "Missing/reordered/duplicate pair rows")
    for row, key in zip(rows, keys):
        close(float(row["distance"]), expected[key])
    return len(keys)


def verify(report_path, repo, source_commit):
    report_path, repo = Path(report_path).resolve(), Path(repo).resolve()
    report_ref = fingerprint(report_path)
    report = json.loads(report_path.read_text())
    checked = [report_ref]

    def bound(ref):
        demand(fingerprint(ref["path"]) == ref, "Changed direct binding")
        if ref not in checked:
            checked.append(ref)
        return Path(ref["path"])

    demand(report["schema"] == "swiss_model_divergence_features_v1"
           and report["status"] == "features_constructed_unverified"
           and report["failed_families"] == [] and report["model"] == "WAG+G4"
           and report["seed"] == 20261007 and report["unit"] == UNIT
           and all(report[k] is False for k in FLAGS)
           and type(report["new_bootstrap_draws"]) is int and report["new_bootstrap_draws"] == 0,
           "Inflated or incomplete feature scope")
    for key, sha in PINS.items():
        demand(report["inputs"][key]["sha256"] == sha, "Changed frozen input")
        bound(report["inputs"][key])
    inventory = json.loads(bound(report["inputs"]["alignments"]).read_text())
    admission = json.loads(bound(report["inputs"]["admission"]).read_text())
    members = admission["family_memberships"]
    genes = [g for f in sorted(members) for g in members[f]]
    demand(report["memberships"] == members and len(members) == 18
           and len(genes) == len(set(genes)) == 563, "Changed canonical universe")
    demand(report["binary"]["sha256"] == BINARY_SHA and report["binary"]["bytes"] == 11333032,
           "Unexpected binary")
    bound(report["binary"])
    preflight = json.loads(bound(report["preflight"]).read_text())
    demand(preflight["memberships"] == members and preflight["inputs"] == report["inputs"]
           and preflight["binary"] == report["binary"] and preflight["unit"] == UNIT
           and preflight["model"] == "WAG+G4" and preflight["seed"] == 20261007
           and preflight["cpus_per_task"] == 2 and preflight["memory_mib"] == 8192
           and preflight["family_timeout_seconds"] == 600
           and preflight["source_commit"] == report["source_commit"]
           and preflight["available_memory_bytes"] >= 8 * 1024 ** 3
           and all(preflight[k] is False for k in FLAGS), "Changed construction preflight")
    demand([Path(r["path"]).relative_to(repo).as_posix() for r in preflight["sources"]] == [
        "benchmark_tools/prepare_swiss_model_divergence.py",
        "benchmark_tools/swiss_model_divergence_batch_20261007.sh",
        "benchmark_tools/results/SWISS_MODEL_DIVERGENCE_PROTOCOL_20261007.md"], "Unexpected source set")
    for ref in preflight["sources"]:
        path = bound(ref)
        demand(path.read_bytes() == subprocess.check_output(["git", "show",
               report["source_commit"] + ":" + path.relative_to(repo).as_posix()], cwd=repo),
               "Source differs from pre-execution commit")
    demand(report["source"] == preflight["sources"][0], "Wrong feature source")
    bound(preflight["python_executable"])
    selected = {r["refog"]: r for r in inventory["runs"]}
    demand([r["family"] for r in report["runs"]] == sorted(members), "Incomplete runs")
    features, pair_index = {}, {}
    for run in report["runs"]:
        family = run["family"]
        demand(run["status"] == "feature_constructed" and run["exit_code"] == 0
               and run["timed_out"] is False and run["alignment"] == selected[family]["alignment"]
               and run["alignment"] in admission["records"]
               and run["columns"] == selected[family]["columns"], "Wrong selected run")
        alignment = list(SeqIO.parse(bound(run["alignment"]), "fasta"))
        demand(sorted(r.id for r in alignment) == members[family]
               and all(len(r.seq) == run["columns"] for r in alignment), "Changed alignment scope")
        directory = report_path.parent / family
        expected_command = [report["binary"]["path"], "-s", run["alignment"]["path"],
                            "--seqtype", "AA", "-m", "WAG+G4", "--seed", "20261007",
                            "-T", "1", "--mem", "4G", "-keep-ident", "--prefix",
                            str(directory / "inference")]
        demand(run["command"] == expected_command and run["finish_ns"] >= run["start_ns"],
               "Changed command or timestamps")
        for ref in run["outputs"]:
            demand(Path(ref["path"]).parent == directory, "Misplaced family output")
            bound(ref)
        demand(directory / "inference.treefile" in [Path(r["path"]) for r in run["outputs"]],
               "Unbound inference tree")
        family_receipt = fingerprint(directory / "result.json")
        demand(json.loads(bound(family_receipt).read_text()) == run, "Changed family receipt")
        descriptor, pairs = edge_distances(directory / "inference.treefile", members[family])
        demand(set(run["features"]) == set(descriptor), "Incorrect feature fields")
        for key, value in descriptor.items():
            if key in ("members", "unit", "pairs"):
                demand(run["features"][key] == value, "Feature descriptor differs")
            else:
                close(run["features"][key], value)
        features[family] = descriptor
        pair_index.update({(family, a, b): distance for (a, b), distance in pairs.items()})
    pair_count = compare_pairs(bound(report["pairs"]), pair_index)
    cutoff = statistics.median(v["median_pair_distance"] for v in features.values())
    close(report["median_family_distance"], cutoff)
    expected_bins = dict(all=sorted(features),
                         lower_or_equal_median=sorted(f for f in features
                             if features[f]["median_pair_distance"] <= cutoff),
                         higher_than_median=sorted(f for f in features
                             if features[f]["median_pair_distance"] > cutoff))
    demand(report["strata"] == expected_bins, "Wrong median/tie bins")
    job = preflight["job_id"]
    accounting = subprocess.check_output(["sacct", "-j", str(job), "--noheader", "--parsable2",
                                          "--format=JobIDRaw,State,ExitCode"], text=True)
    parent = [line.split("|") for line in accounting.splitlines() if line.split("|")[0] == str(job)]
    demand(len(parent) == 1 and parent[0][1:3] == ["COMPLETED", "0:0"],
           "Feature construction job not successfully terminal")
    source = fingerprint(__file__)
    demand(Path(__file__).read_bytes() == subprocess.check_output(["git", "show",
           source_commit + ":benchmark_tools/readback_swiss_model_divergence.py"], cwd=repo),
           "Readback source not committed before execution")
    return dict(schema="swiss_model_divergence_edge_readback_v1", status="features_verified",
                report=report_ref, source=source, source_commit=source_commit, checked_inputs=checked,
                families_checked=18, proteins_checked=563, pairs_checked=pair_count,
                scheduler=dict(job_id=job, state=parent[0][1], exit_code=parent[0][2], raw=accounting),
                features=features, median_family_distance=cutoff, strata=expected_bins,
                limitations=["Independent edge-split arithmetic shares Bio.Phylo parsing.",
                             "Not independent inference, model validation, true history or confirmation."],
                **{k: False for k in FLAGS}, new_bootstrap_draws=0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source-commit", required=True)
    args = parser.parse_args()
    demand(not args.output.exists() and not args.output.is_symlink(), "Existing readback; never overwrite")
    result = verify(args.report, args.repo, args.source_commit)
    with args.output.open("x", encoding="ascii") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps({k: result[k] for k in ("status", "families_checked", "proteins_checked", "pairs_checked")}))
