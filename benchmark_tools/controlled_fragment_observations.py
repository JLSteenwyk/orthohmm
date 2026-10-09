"""Frozen synthetic truncations, pair strata and whole-seed uncertainty."""

from copy import deepcopy
import hashlib
from pathlib import Path

from benchmark_tools.simulation_conditions import canonical_pairs, ranked_names, score_pairs

CONDITION = "fragment20_center60_v1"
SEEDS = tuple(range(20261101, 20261111))
METHODS = ("orthohmm_high_sensitivity", "orthohmm_satellite_v2", "orthofinder_full",
           "orthofinder_sequence_only")
INFERENCE_METHODS = METHODS[:3]
METRICS = ("f1", "precision", "recall")
COMPARISONS = tuple(("fragment", method, "baseline", method) for method in INFERENCE_METHODS) + (
    ("fragment", METHODS[0], "fragment", METHODS[2]),
    ("fragment", METHODS[1], "fragment", METHODS[2]),
)
BOOTSTRAP_REPLICATES = 20000
BOOTSTRAP_SEED = 20261011
PLANNED_ENDPOINTS = 15


def sequence_hash(sequence):
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def observe(sequences, owners, seed):
    """Select by ID only; preserve the full gene and species universe."""
    if type(seed) is not int or seed <= 0:
        raise ValueError("Positive integer seed required")
    if not sequences or set(sequences) != set(owners):
        raise ValueError("Sequence and owner universes differ or are empty")
    if any(not isinstance(owner, str) or not owner for owner in owners.values()):
        raise ValueError("Nonempty species owners required")
    selected = set(ranked_names(list(sequences), seed, CONDITION)[:len(sequences) // 5])
    observed, coordinates = {}, []
    for gene in sorted(sequences):
        sequence = sequences[gene]
        if not isinstance(sequence, str) or not sequence or any(c.isspace() for c in sequence):
            raise ValueError("Nonempty unspaced sequence required")
        fragment = gene in selected
        length = len(sequence)
        retained = 3 * length // 5 if fragment else length
        if fragment and not 0 < retained < length:
            raise ValueError("Unsupported parent length for truncation: " + gene)
        start = (length - retained) // 2 if fragment else 0
        stop = start + retained
        observed[gene] = sequence[start:stop]
        coordinates.append({"gene": gene, "species": owners[gene], "fragment": fragment,
                            "parent_length": length, "observed_length": retained,
                            "start": start, "stop": stop,
                            "parent_sha256": sequence_hash(sequence),
                            "observed_sha256": sequence_hash(observed[gene])})
    return observed, coordinates


def verify_observation(parent, observed, owners, coordinates, seed):
    """Read back exact coordinates and strings without calling the transform."""
    if type(seed) is not int or seed <= 0 or not parent:
        raise ValueError("Nonempty parent and positive integer seed required")
    if set(parent) != set(observed) or set(parent) != set(owners):
        raise ValueError("Observation changed the gene universe")
    expected_order = sorted(parent, key=lambda gene: (
        hashlib.sha256(f"{seed}:{CONDITION}:{gene}".encode("utf-8")).digest(), gene))
    selected = set(expected_order[:len(parent) // 5])
    if len(coordinates) != len(parent) or [row["gene"] for row in coordinates] != sorted(parent):
        raise ValueError("Incomplete, duplicate or unordered coordinate records")
    for row in coordinates:
        gene = row["gene"]
        sequence = parent[gene]
        length = len(sequence)
        fragment = gene in selected
        retained = (3 * length) // 5 if fragment else length
        start = (length - retained) // 2 if fragment else 0
        if not sequence or not 0 < retained <= length or (fragment and retained == length):
            raise ValueError("Unsupported coordinate length")
        expected = {"gene": gene, "species": owners[gene], "fragment": fragment,
                    "parent_length": length, "observed_length": retained,
                    "start": start, "stop": start + retained,
                    "parent_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
                    "observed_sha256": hashlib.sha256(observed[gene].encode("ascii")).hexdigest()}
        if type(row["fragment"]) is not bool or row != expected or observed[gene] != sequence[start:start + retained]:
            raise ValueError("Observation, coordinates or species differ from frozen rule")
    return {"genes": len(parent), "fragment_genes": len(selected), "status": "exact_substrings_verified"}


def score_strata(predictions, truth, owners, flags):
    """Partition ALL cross-species pairs, not just within-family predictions."""
    if set(flags) != set(owners) or any(type(value) is not bool for value in flags.values()):
        raise ValueError("Boolean fragment flags must cover the exact gene universe")
    predictions, truth = list(predictions), list(truth)
    predicted, _ = canonical_pairs(predictions, owners)
    reference, duplicates = canonical_pairs(truth, owners)
    if duplicates:
        raise ValueError("Duplicate truth pairs")
    # Keep the audited kernel's coverage/counts, but expose undefined ratios as null.
    total = score_pairs(predictions, truth, owners)
    for metric in total["undefined_ratios"]:
        total[metric] = None
    strata = []
    for count in range(3):
        p = [pair for pair in predicted if sum(flags[g] for g in pair) == count]
        t = [pair for pair in reference if sum(flags[g] for g in pair) == count]
        scored = score_pairs(p, t, owners)
        for metric in scored["undefined_ratios"]:
            scored[metric] = None
        strata.append({"fragment_endpoints": count, "score": scored})
    for name in ("tp", "fp", "fn", "predicted_pairs", "eligible_true_pairs"):
        if total[name] != sum(row["score"][name] for row in strata):
            raise ValueError("Fragment strata do not partition counts")
    return {"score": total, "strata": strata}


def paired_intervals(records):
    """Fixed 15-endpoint family; failures never become zero-valued observations."""
    import math
    import numpy as np

    expected = {(arm, seed, method) for arm in ("baseline", "fragment") for seed in SEEDS for method in METHODS}
    indexed = {}
    for row in records:
        key = row["arm"], row["seed"], row["method"]
        if key not in expected or key in indexed:
            raise ValueError("Unknown or duplicate fragment comparison record")
        if row["status"] not in {"complete", "failed", "unavailable"}:
            raise ValueError("Nonterminal fragment comparison record")
        if row["status"] != "complete" and "score" in row:
            raise ValueError("Failed outcome must not carry a score")
        if row["status"] == "complete":
            for metric in METRICS:
                value = row["score"][metric]
                if value is not None and (isinstance(value, bool) or not math.isfinite(value) or not 0 <= value <= 1):
                    raise ValueError("Invalid fragment metric")
        indexed[key] = row
    if set(indexed) != expected:
        raise ValueError("Incomplete 80-record baseline/fragment inventory")
    result = []
    for target_arm, target_method, reference_arm, reference_method in COMPARISONS:
        for metric in METRICS:
            eligible, excluded, values = [], [], []
            for seed in SEEDS:
                target = indexed[target_arm, seed, target_method]
                reference = indexed[reference_arm, seed, reference_method]
                reasons = []
                for role, row in (("target", target), ("reference", reference)):
                    if row["status"] != "complete":
                        reasons.append({"role": role, "status": row["status"], "reason": row.get("reason", row["status"])})
                    elif row["score"][metric] is None:
                        reasons.append({"role": role, "status": "undefined_metric", "reason": metric})
                if reasons:
                    excluded.append({"seed": seed, "reasons": reasons})
                else:
                    eligible.append(seed)
                    values.append(target["score"][metric] - reference["score"][metric])
            row = {"target_arm": target_arm, "target_method": target_method,
                   "reference_arm": reference_arm, "reference_method": reference_method,
                   "metric": metric, "eligible_seeds": eligible, "excluded_seeds": excluded,
                   "paired_differences": values, "planned_endpoints": PLANNED_ENDPOINTS,
                   "replicates": BOOTSTRAP_REPLICATES, "rng": "PCG64", "rng_seed": BOOTSTRAP_SEED,
                   "statistic": "mean paired seed-level metric difference",
                   "quantile_method": "linear"}
            if values:
                array = np.asarray(values, dtype=float)
                rng = np.random.Generator(np.random.PCG64(BOOTSTRAP_SEED))
                draws = rng.integers(0, len(values), size=(BOOTSTRAP_REPLICATES, len(values)))
                replicates = array[draws].mean(axis=1)
                row.update(status="conditional_approximation", estimate=float(array.mean()),
                           nominal_interval=np.quantile(replicates, [.025, .975], method="linear").tolist(),
                           adjusted_interval=np.quantile(replicates, [.05 / 30, 1 - .05 / 30], method="linear").tolist())
            else:
                row.update(status="unavailable", estimate=None, nominal_interval=None, adjusted_interval=None)
            result.append(row)
    return result


def overlaps(first, second):
    first, second = Path(first).resolve(), Path(second).resolve()
    return first == second or first in second.parents or second in first.parents


def fresh_methods(baseline, inputs, destination):
    """Copy admitted argv; change only FASTA, output and metrics paths."""
    inputs, destination = Path(inputs).resolve(), Path(destination).resolve()
    old_input = Path(baseline["input"]).resolve()
    if overlaps(inputs, old_input) or overlaps(inputs, destination):
        raise ValueError("New inputs overlap baseline inputs or new inference outputs")
    if destination.exists():
        raise FileExistsError("Existing inference destination; no automatic restart")
    for original in baseline["methods"].values():
        if overlaps(destination, original["output"]) or overlaps(inputs, original["output"]):
            raise ValueError("New paths overlap original inference artifacts")
        if original.get("metrics") and (overlaps(destination, original["metrics"]) or overlaps(inputs, original["metrics"])):
            raise ValueError("New paths overlap original metrics")
    if overlaps(destination, old_input):
        raise ValueError("New inference destination overlaps baseline inputs")
    result = deepcopy(baseline["methods"])
    for method in INFERENCE_METHODS:
        original = baseline["methods"][method]
        config = result[method]
        output = destination / method
        argv = config["argv"]
        if method.startswith("orthohmm_"):
            if len(argv) < 5 or argv[2:5] != [baseline["input"], original["output"], original["metrics"]]:
                raise ValueError("Unexpected baseline HMM positional paths")
            metrics = destination / (method + ".json")
            argv[2:5] = [str(inputs), str(output), str(metrics)]
            config["metrics"] = str(metrics)
        else:
            if argv.count("-f") != 1 or original["copy_inputs_from"] != baseline["input"]:
                raise ValueError("Unexpected baseline OrthoFinder input command")
            index = argv.index("-f") + 1
            if index >= len(argv) or argv[index] != original["copy_inputs_to"]:
                raise ValueError("Baseline copied FASTA path differs from command")
            copy = output / "input"
            argv[index] = str(copy)
            config.update(copy_inputs_from=str(inputs), copy_inputs_to=str(copy))
        config["output"] = str(output)
    result["orthofinder_sequence_only"]["output"] = result["orthofinder_full"]["output"]
    return result
