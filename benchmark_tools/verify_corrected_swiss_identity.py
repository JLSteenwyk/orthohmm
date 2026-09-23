"""Verify frozen SwissTrees alignments and independently recompute identity bins."""

import argparse
import itertools
import json
import math
from pathlib import Path
import statistics
import subprocess

from Bio import SeqIO
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.verify_ygob_validation import require_completed_job

REPORT_SHA = "7bec0f40c8227f230bbf06cb52622947e256f2edfc2c88b991393fa570603e44"
CANONICAL = frozenset("ACDEFGHIKLMNPQRSTVWY")


def fasta(path):
    rows = list(SeqIO.parse(path, "fasta"))
    result = {row.id: str(row.seq).upper() for row in rows}
    if not rows or len(rows) != len(result):
        raise ValueError("Empty or duplicate FASTA identifiers")
    return result


def scalar_identity(aligned):
    if len(aligned) < 2 or len({len(s) for s in aligned.values()}) != 1 or not all(aligned.values()):
        raise ValueError("Invalid aligned sequence dimensions")
    values, missing = [], 0
    for a, b in itertools.combinations(sorted(aligned), 2):
        comparable = [(x, y) for x, y in zip(aligned[a], aligned[b]) if x in CANONICAL and y in CANONICAL]
        if comparable:
            values.append(sum(x == y for x, y in comparable) / len(comparable))
        else:
            missing += 1
    return dict(pairs=len(values) + missing, pairs_without_canonical_overlap=missing,
        mean_pairwise_identity=None if missing else math.fsum(values) / len(values))


def verify_family(run, sequences):
    if run["status"] != "alignment_validated" or run["exit_code"] != 0 or run["accuracy_evaluated"] is not False:
        raise ValueError("Unsuccessful or scored family")
    expected = {}
    for gene, sequence in sequences.items():
        if not sequence or not set(sequence) <= CANONICAL | {"X", "*"}:
            raise ValueError("Unexpected raw residues")
        expected[gene] = sequence.replace("*", "")
        if not expected[gene]:
            raise ValueError("Empty normalized sequence")
    for key in ("input", "alignment"):
        check(run[key])
    if fasta(run["input"]["path"]) != expected:
        raise ValueError("Alignment input differs from corrected sequences")
    aligned = fasta(run["alignment"]["path"])
    if {gene: seq.replace("-", "") for gene, seq in aligned.items()} != expected:
        raise ValueError("Alignment changed sequence residues or membership")
    observed = scalar_identity(aligned)
    for key in ("pairs", "pairs_without_canonical_overlap"):
        if type(run["identity"][key]) is not int or observed[key] != run["identity"][key]:
            raise ValueError("Pair inventory differs")
    old, new = run["identity"]["mean_pairwise_identity"], observed["mean_pairwise_identity"]
    if (old is None) != (new is None) or (new is not None and
            (isinstance(old, bool) or not math.isclose(old, new, rel_tol=0, abs_tol=1e-12))):
        raise ValueError("Independent identity differs")
    removed = {gene: seq.count("*") for gene, seq in sequences.items() if "*" in seq}
    if (run["removed_stop_symbols"] != removed or run["genes"] != len(aligned)
            or run["columns"] != len(next(iter(aligned.values())))):
        raise ValueError("Alignment dimensions or stop inventory differ")
    return observed


def verify(root, output):
    root = root.resolve()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    report_path = root / "benchmark_tools/results/corrected_swiss_identity_prepared_22102.json"
    report = read_frozen(report_path, REPORT_SHA)
    preflight = report["preflight"]
    directory = root / "benchmarks/results/corrected_swiss_identity_v1"
    if (report["status"] != "corrected_swiss_alignments_prepared_unscored" or report["failed_families"]
            or report["prediction_statistics_evaluated"] is not False or report["publication_ready"] is not False
            or preflight["job_id"] != "22102" or preflight["workers"] != 4):
        raise ValueError("Unexpected frozen panel")
    if json.loads((directory / "manifest.json").read_text()) != report or json.loads((directory / "preflight.json").read_text()) != preflight:
        raise ValueError("Native manifests differ")
    accounting = subprocess.check_output(["sacct", "-j", "22102", "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed"], text=True)
    scheduler = require_completed_job(accounting, 22102)
    tracked = [record(report_path), record(directory / "manifest.json"), record(directory / "preflight.json"),
               record(__file__), *preflight["inputs"], *preflight["tools"], *preflight["fasta_inputs"]]
    for item in tracked:
        check(item)
    inventory = json.loads(Path(preflight["inputs"][0]["path"]).read_text())
    families = inventory["family_memberships"]
    if families != preflight["families"] or inventory["fasta_inputs"] != preflight["fasta_inputs"]:
        raise ValueError("Corrected membership or input inventory differs")
    wanted = {g for genes in families.values() for g in genes}
    sequences = {}
    for item in preflight["fasta_inputs"]:
        for row in SeqIO.parse(item["path"], "fasta"):
            gene = row.id.split("|")[1]
            if gene in wanted:
                if gene in sequences:
                    raise ValueError("Duplicate corrected reference accession")
                sequences[gene] = str(row.seq).upper()
    if set(sequences) != wanted or len(sequences) != 563 or len(families) != 18:
        raise ValueError("Incomplete corrected reference universe")
    if sorted(r["refog"] for r in report["runs"]) != sorted(families):
        raise ValueError("Incomplete or duplicate alignment panel")
    environment = json.loads(Path(preflight["inputs"][1]["path"]).read_text())
    mafft = environment["tool_entrypoints"]["mafft"]["absolute_path"]
    features = {}
    for run in report["runs"]:
        name = run["refog"]
        expected_input = directory / name / "input.faa"
        expected_alignment = directory / name / "aligned.faa"
        if (run["command"] != [mafft, "--amino", "--auto", "--thread", "1", str(expected_input)]
                or run["input"]["path"] != str(expected_input) or run["alignment"]["path"] != str(expected_alignment)):
            raise ValueError("Alignment command or paths differ")
        status = directory / name / "status.json"
        if json.loads(status.read_text()) != run:
            raise ValueError("Per-family status differs")
        features[name] = verify_family(run, {g: sequences[g] for g in families[name]})
        tracked.extend([run["input"], run["alignment"], record(status), record(directory / name / "mafft.stderr")])
    # Compare the independent scalar values, but preserve frozen bin membership.
    available = [r["mean_pairwise_identity"] for r in features.values() if r["mean_pairwise_identity"] is not None]
    median = statistics.median(available) if available else None
    bins = dict(lower_identity=[], higher_identity=[], missing_identity=[])
    for name, row in sorted(features.items()):
        value = row["mean_pairwise_identity"]
        bins["missing_identity" if value is None else "lower_identity" if value <= median else "higher_identity"].append(name)
    frozen = report["median_family_identity"]
    if bins != report["strata"] or (median is None) != (frozen is None) or (median is not None and not math.isclose(median, frozen, rel_tol=0, abs_tol=1e-12)):
        raise ValueError("Frozen bins differ from independent identities")
    for item in tracked:
        check(item)
    result = dict(status="corrected_swiss_identity_features_verified", scheduler=scheduler,
        records=tracked, family_memberships=families, identity_features=features, strata=report["strata"],
        median_family_identity=frozen, scalar_median_family_identity=median,
        absolute_tolerance=1e-12, prediction_statistics_evaluated=False, publication_ready=False,
        limitations=["Independent scalar arithmetic shares the Biopython FASTA parser, not an independent alignment method.",
            "Identity is alignment-dependent, not calibrated evolutionary distance or fragment annotation.",
            "No prediction outcomes, subgroup significance or generalization claims are evaluated."])
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    verify(args.root, args.output)
