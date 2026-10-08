"""Recompute planned paired intervals from available admitted native counts."""

import argparse
from copy import deepcopy
import json
from pathlib import Path

import numpy as np

from benchmark_tools import bind_composed_native_qfo_swiss_uncertainty as binding
from benchmark_tools.bootstrap_qfo_factorial import CELLS, contrasts, validated_values
from benchmark_tools.bootstrap_qfo_swiss_stages import METRICS, aggregate
from benchmark_tools.export_native_factorial_progress import load, require
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

REPLICATES = 100000
SEED = 20260922
MULTIPLICITY = 42


def analyze(native, retained):
    require(native and set(native).issubset(CELLS), "Require known native cells")
    validated_values({**retained, "status": "qfo_factorial_swiss_counts_verified"})
    validation = deepcopy(retained)
    validation["status"] = "qfo_factorial_swiss_counts_verified"
    for row in validation["cells"]:
        if row["cell"] in native:
            actual = native[row["cell"]]
            require(actual["cell"] == row["cell"], "Native cell key differs")
            require(len(actual["families"]) == len(row["families"]), "Native family count differs")
            for fresh, previous in zip(actual["families"], row["families"]):
                fresh_counts, old_counts = fresh["counts_without_prior"], previous["counts_without_prior"]
                require(fresh["family"] == previous["family"]
                    and sorted(fresh["represented_genes"]) == sorted(previous["represented_genes"])
                    and fresh_counts["TP"] + fresh_counts["FN"] == old_counts["TP"] + old_counts["FN"]
                    and fresh_counts["FP"] + fresh_counts["TN"] == old_counts["FP"] + old_counts["TN"],
                    "Native reference members or truth totals differ")
            row.update(families=actual["families"], aggregate=actual["aggregate"])
    # Retained cells establish the fixed reference universe for validation only.
    # Discard their values before resampling or computing any native contrast.
    checked = validated_values(validation)
    selected = [cell for cell in CELLS if cell in native]
    values = checked[[CELLS.index(cell) for cell in selected]]
    definitions = contrasts()
    needed = {row["name"]: [cell for cell, weight in zip(CELLS, row["weights"]) if weight]
        for row in definitions}
    require(any(set(cells).issubset(native) for cells in needed.values()),
        "No estimable native contrast; do not generate unneeded draws")
    families = retained["families"]
    n = len(families)
    multiplicities = np.random.Generator(np.random.PCG64(SEED)).multinomial(
        n, np.full(n, 1 / n), size=REPLICATES)
    draws = np.asarray([aggregate(multiplicities @ cell / n) for cell in values])
    points, per_family = aggregate(values.mean(axis=1)), aggregate(values)
    rows = []
    for definition in definitions:
        absent = [cell for cell in needed[definition["name"]] if cell not in native]
        row = dict(definition, status="native_counts_unavailable" if absent else "native_counts_resampled",
            missing_cells=absent, metrics=None, family_differences=None)
        if not absent:
            weights = np.asarray([definition["weights"][CELLS.index(cell)] for cell in selected])
            delta = np.tensordot(weights, draws, axes=1)
            point = weights @ points
            differences = np.tensordot(weights, per_family, axes=1)
            row["metrics"] = {metric: dict(difference=float(point[j]),
                paired_percentile_ci=np.quantile(delta[:, j], [.025, .975], method="linear").tolist(),
                bonferroni_percentile_ci=np.quantile(delta[:, j],
                    [.05 / (2 * MULTIPLICITY), 1 - .05 / (2 * MULTIPLICITY)], method="linear").tolist(),
                family_wins=int(np.sum(differences[:, j] > 1e-10)),
                family_ties=int(np.sum(np.abs(differences[:, j]) <= 1e-10)),
                family_losses=int(np.sum(differences[:, j] < -1e-10)))
                for j, metric in enumerate(METRICS)}
            row["family_differences"] = [dict(family=family, **dict(zip(METRICS, value.tolist())))
                for family, value in zip(families, differences)]
        rows.append(row)
    return dict(status="paired_composed_native_qfo_swiss_intervals", families=families,
        observed_cells=selected, point_estimates={cell: dict(zip(METRICS, point.tolist()))
            for cell, point in zip(selected, points)}, comparisons=rows,
        replicates=REPLICATES, seed=SEED, alpha=.05, multiplicity_endpoints=MULTIPLICITY,
        quantile_method="linear", rng="numpy.PCG64 multinomial; shared draws across observed native cells",
        numpy_version=np.__version__, units="raw 0-to-1 metric units",
        new_bootstrap_draws=REPLICATES, retained_intervals_reused=False,
        unobserved_cells_imputed=False, independent_confirmation=False,
        new_accuracy_or_resource_admission=False, publication_ready=False)


def run(snapshot_path, snapshot_sha, audits, retained_path, bootstrap_path):
    reviewed = binding.bind(snapshot_path, snapshot_sha, audits, retained_path, bootstrap_path)
    evidence = [*reviewed["evidence"], reviewed["snapshot"], reviewed["retained_counts"],
        reviewed["bootstrap"], reviewed["source"], *reviewed["helpers"]]
    retained, _ = load(retained_path, reviewed["retained_counts"]["sha256"], evidence)
    native = {}
    for path, digest in audits:
        audit, ref = load(path, digest, evidence)
        for row in audit["cells"]:
            require(reviewed["bound_cells"][row["cell"]]["count_audit"] == ref,
                "Bootstrap input differs from checked count audit")
            native[row["cell"]] = row
    result = analyze(native, retained)
    helpers = [record(module.__file__) for module in (
        binding, binding.composed, binding.original)]
    helpers.extend(record(Path(__file__).with_name(name)) for name in (
        "bootstrap_qfo_factorial.py", "bootstrap_qfo_swiss_stages.py", "audit_qfo_swiss_counts.py"))
    source = record(__file__)
    for ref in [*evidence, *helpers, source]:
        check(ref)
    result.update(schema="composed_native_qfo_swiss_bootstrap_v1", source=source,
        helpers=helpers, evidence=evidence, snapshot=reviewed["snapshot"],
        retained_counts=reviewed["retained_counts"], retained_bootstrap=reviewed["bootstrap"],
        bound_cells=reviewed["bound_cells"], limitations=[
            "New shared family draws use actual admitted native counts, not old intervals or missing cached cells.",
            "Retained counts establish reference membership/truth and existing provenance, not unobserved native predictions.",
            "All 42 planned endpoints remain in adjustment; missing contrasts retain null metrics.",
            "Harmonic mean of family-mean PPV and TPR is recomputed in every replicate, not mean family F1.",
            "Only 18 development-exposed families; exchangeability and approximate percentile coverage limits remain.",
            "SwissTrees only; no uncertainty for other challenges or the secondary six-metric mean.",
            "New conditional uncertainty is not independent confirmation, raw readmission or isolated efficiency evidence."])
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("snapshot", "retained-counts", "bootstrap", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--snapshot-sha256", required=True)
    parser.add_argument("--counts-audit", nargs=2, action="append", required=True)
    args = parser.parse_args()
    require(not args.output.exists() and not args.output.is_symlink(), "Output already exists")
    result = run(args.snapshot, args.snapshot_sha256, args.counts_audit, args.retained_counts, args.bootstrap)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps(dict(observed_cells=len(result["observed_cells"]), new_draws=REPLICATES,
        estimable_contrasts=sum(row["metrics"] is not None for row in result["comparisons"]))))


if __name__ == "__main__":
    main()
