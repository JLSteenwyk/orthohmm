"""Replay frozen YGOB arithmetic from identifier-free retained pillar counts."""

import argparse
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import sys

import numpy as np


SUMMARY_SHA = "3927d81b1fd851eef3c8ce1343679558ef4f93bcfe2e643435e369587a44a4ef"
FULL_SHA = "5104e000d0c3bb8687103646be8011fbc97f8b05f7eb0a925bbfdd1e76200e83"
METHODS = ("orthohmm_satellite_v2", "orthohmm_high_sensitivity", "orthofinder_full", "orthofinder_sequence_only")
INFERENTIAL = METHODS[:3]
METRICS = ("f1", "precision", "recall")
COLUMNS = ("tp", "fp_twice", "fn", "covered_genes", "exact")
PILLARS, GENES, INPUT_GENES = 10250, 83391, 83404
REPLICATES, SEED = 20000, 20260917
MAX_DECODED = 4 * 1024 ** 2


def identity(data):
    return dict(bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def read_pinned(path, digest, limit):
    path = Path(path)
    if path.is_symlink() or not path.is_file() or path.stat().st_size > limit:
        raise ValueError("Require a bounded direct input file")
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != digest:
        raise ValueError("Input digest differs from its retained identity")
    return data


def require_summary(summary):
    u = summary["uncertainty"]
    if (summary["schema_version"] != 1 or summary["publication_ready"] is not False
            or summary["frozen_evaluation_gates_verified"] is not True
            or summary["completion_gates_verified_by_this_module"] is not False
            or summary["input_gene_count"] != INPUT_GENES
            or summary["protocol"] != "benchmark_tools/results/YGOB_VALIDATION_PROTOCOL_20260916.md"
            or summary["primary_contrast"] != "orthohmm_satellite_v2 versus orthofinder_full"
            or summary["secondary_contrast"] != "orthohmm_high_sensitivity versus orthofinder_full"
            or summary["diagnostic_methods"] != [METHODS[3]] or set(summary["scores"]) != set(METHODS)
            or summary["full_results"]["sha256"] != FULL_SHA or summary["full_results"]["bytes"] != 7945362
            or (u["baseline"], u["replicates"], u["seed"], u["multiplicity_count"], u["rng"], u["numpy_version"])
               != (METHODS[2], REPLICATES, SEED, 6, "PCG64 multinomial", "2.2.6")
            or set(u["comparisons"]) != set(METHODS[:2])):
        raise ValueError("Frozen YGOB study or inferential controls differ")
    for method in METHODS[:2]:
        if set(u["comparisons"][method]) != set(METRICS):
            raise ValueError("Wrong six-endpoint inventory")
        for item in u["comparisons"][method].values():
            if set(item) != {"difference_percentage_points", "paired_95_percent_ci", "bonferroni_ci"}:
                raise ValueError("Wrong endpoint fields")
            for key in ("paired_95_percent_ci", "bonferroni_ci"):
                values = item[key]
                if len(values) != 2 or any(type(v) not in (int, float) or not math.isfinite(v) for v in values) or values[0] > values[1]:
                    raise ValueError("Invalid interval")
            value = item["difference_percentage_points"]
            if type(value) not in (int, float) or not math.isfinite(value):
                raise ValueError("Invalid contrast value")


def statistics(twice_counts):
    tp, fp, fn = np.moveaxis(np.asarray(twice_counts, dtype=np.int64), -1, 0)
    numerators = (2 * tp, tp, tp)
    denominators = (2 * tp + fp + fn, tp + fp, tp + fn)
    return np.stack([np.divide(n, d, out=np.zeros_like(n, dtype=float), where=d > 0)
                     for n, d in zip(numerators, denominators)], axis=-1)


def close(actual, expected):
    if not np.allclose(actual, expected, rtol=0, atol=1e-12):
        raise ValueError("YGOB arithmetic differs from the retained result")


def verify_points(summary, snapshot):
    require_summary(summary)
    if (set(snapshot) != {"schema", "columns", "summary", "full_results", "pillar_sizes", "methods",
                         "ordered_original_pillar_signature_sha256", "raw_reference_or_prediction_memberships_included", "publication_ready"}
            or snapshot["schema"] != "identifier-free-ygob-counts-v1"
            or snapshot["columns"] != list(COLUMNS)
            or snapshot["summary"] != dict(bytes=7978, sha256=SUMMARY_SHA)
            or snapshot["full_results"] != dict(bytes=7945362, sha256=FULL_SHA)
            or snapshot["raw_reference_or_prediction_memberships_included"] is not False
            or snapshot["publication_ready"] is not False
            or set(snapshot["methods"]) != set(METHODS)):
        raise ValueError("Wrong identifier-free observation scope")
    sizes = snapshot["pillar_sizes"]
    if len(sizes) != PILLARS or any(type(n) is not int or n < 1 for n in sizes) or sum(sizes) != GENES:
        raise ValueError("Wrong pillar size inventory")
    signature = snapshot["ordered_original_pillar_signature_sha256"]
    if not isinstance(signature, str) or len(signature) != 64 or any(c not in "0123456789abcdef" for c in signature):
        raise ValueError("Invalid source-order signature")
    arrays = {}
    for method in METHODS:
        rows = snapshot["methods"][method]
        if len(rows) != PILLARS:
            raise ValueError("Missing or extra pillar observations")
        for n, row in zip(sizes, rows):
            if len(row) != 5 or any(type(v) is not int or v < 0 for v in row):
                raise ValueError("Counts must be nonnegative integers, not booleans or rounded half-counts")
            tp, fp2, fn, covered, exact = row
            if (tp + fn != n * (n - 1) // 2 or covered > n or tp > covered * (covered - 1) // 2
                    or fp2 > n * (GENES - n) or exact not in (0, 1)
                    or bool(exact) != (covered == n and fn == 0 and fp2 == 0)):
                raise ValueError("Inconsistent pillar counts, coverage or exact recovery")
        matrix = np.asarray(rows, dtype=np.int64)
        twice = matrix[:, :3] * (2, 1, 2)
        totals = twice.sum(axis=0)
        retained = summary["scores"][method]
        independent = summary["independent_enumerated_counts"][method]
        close(totals / 2, [retained["counts"][k] for k in ("tp", "fp", "fn")])
        close(totals / 2, [independent[k] for k in ("tp", "fp", "fn")])
        close(statistics(totals), [retained["metrics"][m] for m in METRICS])
        covered, exact = int(matrix[:, 3].sum()), int(matrix[:, 4].sum())
        if (retained["reference_genes"] != GENES or retained["reference_groups"] != PILLARS
                or retained["covered_reference_genes"] != covered or retained["exact_reference_groups"] != exact
                or retained["represented_reference_groups"] != int(np.sum(matrix[:, 3] > 0))):
            raise ValueError("Coverage or exact-pillar totals differ")
        close(retained["reference_gene_coverage"], covered / GENES)
        before, projected = retained["predicted_groups_before_projection"], retained["predicted_groups_with_scored_genes"]
        if (type(before) is not int or type(projected) is not int or not 0 <= projected <= before
                or before < 1 or retained["excluded_predicted_genes"] != INPUT_GENES - covered
                or retained["predicted_group_coverage_definition"] != "fraction of supplied groups retaining at least one scored gene"):
            raise ValueError("Projection metadata differs")
        close(retained["predicted_group_coverage"], projected / before)
        arrays[method] = twice
    return arrays


def paired_draws(arrays, replicates, seed, batch_size=128):
    if type(batch_size) is not int or batch_size < 1 or batch_size > 128:
        raise ValueError("Use a bounded batch size from 1 to 128")
    n = len(next(iter(arrays.values())))
    rng = np.random.Generator(np.random.PCG64(seed))
    samples = {m: np.empty((replicates, 3)) for m in arrays}
    for start in range(0, replicates, batch_size):
        end = min(start + batch_size, replicates)
        weights = rng.multinomial(n, np.full(n, 1 / n), size=end - start)
        for method, counts in arrays.items():
            samples[method][start:end] = statistics(weights @ counts)
    return samples


def verify_intervals(summary, arrays, batch_size=128):
    if np.__version__ != summary["uncertainty"]["numpy_version"]:
        raise ValueError("Use the retained NumPy 2.2.6 numerical implementation")
    samples = paired_draws({m: arrays[m] for m in INFERENTIAL}, REPLICATES, SEED, batch_size)
    baseline = METHODS[2]
    for method in METHODS[:2]:
        delta = 100 * (samples[method] - samples[baseline])
        point = 100 * (statistics(arrays[method].sum(axis=0)) - statistics(arrays[baseline].sum(axis=0)))
        for i, metric in enumerate(METRICS):
            retained = summary["uncertainty"]["comparisons"][method][metric]
            close(retained["difference_percentage_points"], point[i])
            close(retained["paired_95_percent_ci"], np.quantile(delta[:, i], [.025, .975], method="linear"))
            close(retained["bonferroni_ci"], np.quantile(delta[:, i], [.025 / 6, 1 - .025 / 6], method="linear"))
    return 6


def export_counts(full_path, summary_path, output):
    output = Path(output)
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    full_bytes = read_pinned(full_path, FULL_SHA, 8 * 1024 ** 2)
    summary_bytes = read_pinned(summary_path, SUMMARY_SHA, 8192)
    full, summary = json.loads(full_bytes), json.loads(summary_bytes)
    require_summary(summary)
    if set(full) != set(summary) - {"full_results"} or any(full[k] != summary[k] for k in full if k != "scores"):
        raise ValueError("Full result and admitted summary differ")
    base = full["scores"][METHODS[2]]["records"]
    signature = [(r["pillar"], r["genes"]) for r in base]
    if len({name for name, _ in signature}) != PILLARS or signature != sorted(signature):
        raise ValueError("Source pillars must be unique and in original sorted order")
    snapshot = dict(schema="identifier-free-ygob-counts-v1", columns=list(COLUMNS),
                    summary=identity(summary_bytes), full_results=identity(full_bytes),
                    pillar_sizes=[r["genes"] for r in base], methods={},
                    ordered_original_pillar_signature_sha256=hashlib.sha256(json.dumps(signature, separators=(",", ":")).encode()).hexdigest(),
                    raw_reference_or_prediction_memberships_included=False, publication_ready=False)
    for method in METHODS:
        score = full["scores"][method]
        if {k: v for k, v in score.items() if k != "records"} != summary["scores"][method]:
            raise ValueError("Full result point metadata differs")
        if [(r["pillar"], r["genes"]) for r in score["records"]] != signature:
            raise ValueError("Paired pillar identities/order differ")
        rows = []
        for r in score["records"]:
            fp2 = 2 * r["fp"]
            if not math.isfinite(fp2) or fp2 != int(fp2) or type(r["exact"]) is not bool:
                raise ValueError("Invalid retained half-count or exact flag")
            rows.append([r["tp"], int(fp2), r["fn"], r["covered_genes"], int(r["exact"])])
        snapshot["methods"][method] = rows
    verify_points(summary, snapshot)
    raw = (json.dumps(snapshot, separators=(",", ":"), sort_keys=True, allow_nan=False) + "\n").encode()
    if len(raw) > MAX_DECODED:
        raise ValueError("Observation payload exceeds its bounded size")
    payload = gzip.compress(raw, compresslevel=6, mtime=0)
    if read_pinned(full_path, FULL_SHA, 8 * 1024 ** 2) != full_bytes or read_pinned(summary_path, SUMMARY_SHA, 8192) != summary_bytes:
        raise ValueError("Source changed during projection")
    with output.open("xb") as stream:
        stream.write(payload)
    return dict(status="identifier_free_ygob_counts_exported", compressed=identity(payload), decoded=identity(raw),
                pillars=PILLARS, methods=4, integer_cells=PILLARS * 5 * 4, raw_memberships_included=False,
                native_inference_repeated=False, intervals_recomputed=False, publication_ready=False)


def reproduce(snapshot_path, snapshot_sha, summary_path, output, batch_size=128):
    output = Path(output)
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    result = dict(status="validation_failed", publication_ready=False, native_inference_repeated=False,
                  raw_scoring_repeated=False, historical_evidence_paths_accessed=False)
    try:
        compressed = read_pinned(snapshot_path, snapshot_sha, 1024 ** 2)
        summary_bytes = read_pinned(summary_path, SUMMARY_SHA, 8192)
        with gzip.GzipFile(fileobj=io.BytesIO(compressed)) as stream:
            raw = stream.read(MAX_DECODED + 1)
        if len(raw) > MAX_DECODED:
            raise ValueError("Decoded count payload exceeds its bounded size")
        summary, snapshot = json.loads(summary_bytes), json.loads(raw)
        arrays = verify_points(summary, snapshot)
        endpoints = verify_intervals(summary, arrays, batch_size)
        if read_pinned(snapshot_path, snapshot_sha, 1024 ** 2) != compressed or read_pinned(summary_path, SUMMARY_SHA, 8192) != summary_bytes:
            raise ValueError("Input changed during numerical replay")
        result.update(status="frozen_ygob_arithmetic_reproduced", snapshot=identity(compressed),
                      decoded=identity(raw), summary=identity(summary_bytes), source=identity(Path(__file__).read_bytes()),
                      methods_checked=4, pillars_checked=PILLARS, point_metrics_checked=12,
                      interval_endpoints_checked=endpoints, nominal_and_adjusted_bounds_checked=24,
                      replicates=REPLICATES, seed=SEED, batch_size=batch_size, numpy_version=np.__version__,
                      python=sys.version, absolute_tolerance=1e-12,
                      limitations=["Replays retained sufficient statistics, not raw reference/prediction or native admission.",
                          "Uses independent integer count arithmetic, but the same NumPy PCG64 and quantiles, not a new statistical engine.",
                          "Pillar exchangeability and half-allocation of cross-pillar errors remain assumptions.",
                          "Novel-taxon transfer retains homolog overlap; no family-disjoint or unrestricted generalization claim.",
                          "Projection-group metadata is checked internally, not reconstructed from omitted memberships.",
                          "No controlled resources, complete release, redistribution clearance or archival deposition."])
    except Exception as error:
        result.update(error_type=type(error).__name__, error=str(error))
        raise
    finally:
        with output.open("x") as stream:
            json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    exporter = commands.add_parser("export")
    exporter.add_argument("--full-results", type=Path, required=True)
    exporter.add_argument("--summary", type=Path, required=True)
    exporter.add_argument("--output", type=Path, required=True)
    checker = commands.add_parser("reproduce")
    checker.add_argument("--snapshot", type=Path, required=True)
    checker.add_argument("--snapshot-sha256", required=True)
    checker.add_argument("--summary", type=Path, required=True)
    checker.add_argument("--output", type=Path, required=True)
    checker.add_argument("--batch-size", type=int, default=128)
    args = parser.parse_args()
    result = (export_counts(args.full_results, args.summary, args.output) if args.command == "export"
              else reproduce(args.snapshot, args.snapshot_sha256, args.summary, args.output, args.batch_size))
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
