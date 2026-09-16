"""Label-free coverage diagnostics for frozen HMM and DIAMOND search controls."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_historical_profile_ablation import verify_file
from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.recover_factorial_postflight import PREPARED_HASH
from benchmark_tools.replay_high_sensitivity import load_replay_input
from benchmark_tools.run_sequence_search_control import MANIFEST_SHA, verify_plan
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.verify_ygob_validation import require_completed_job


def same_species_partition(first, second):
    if len(first) != len(second):
        raise ValueError("Different species universe lengths")
    pairs = np.unique(np.column_stack((first, second)), axis=0)
    if len(pairs) != len(np.unique(first)) or len(pairs) != len(np.unique(second)):
        raise ValueError("Species ownership differs between hit inputs")


def summarize(queries, targets, scores, species):
    queries, targets, scores, species = map(np.asarray, (queries, targets, scores, species))
    n = len(species)
    if (not n or any(a.ndim != 1 for a in (queries, targets, scores, species))
            or len(queries) != len(targets) or len(queries) != len(scores)
            or queries.dtype.kind not in "iu" or targets.dtype.kind not in "iu"
            or species.dtype.kind not in "iu" or np.any(species < 0)
            or not np.isfinite(scores).all() or np.any(scores <= 0)
            or np.any(queries < 0) or np.any(targets < 0)
            or np.any(queries >= n) or np.any(targets >= n)):
        raise ValueError("Invalid hit arrays or species universe")
    codes = np.sort(queries.astype(np.int64) * n + targets)
    if len(codes) > 1 and np.any(codes[1:] == codes[:-1]):
        raise ValueError("Duplicate directed query-target pair")
    self_mask = queries == targets
    cross = species[queries] != species[targets]
    nonself_codes = codes[codes // n != codes % n]
    reverse = (nonself_codes % n) * n + nonself_codes // n
    reciprocal = len(np.intersect1d(nonself_codes, reverse, assume_unique=True))
    outgoing = np.bincount(queries, minlength=n)
    cross_outgoing = np.bincount(queries[cross], minlength=n)
    probabilities = [0, .25, .5, .75, .95, .99, 1]
    report = {"genes": n, "directed_hits": len(codes), "self_hits": int(self_mask.sum()),
              "nonself_hits": len(nonself_codes), "cross_species_hits": int(cross.sum()),
              "queries_without_hits": int(np.count_nonzero(outgoing == 0)),
              "queries_without_nonself_hits": int(n - len(np.unique(queries[~self_mask]))),
              "queries_without_cross_species_hits": int(np.count_nonzero(cross_outgoing == 0)),
              "targets_without_hits": int(n - len(np.unique(targets))),
              "reciprocal_nonself_directed_hits": reciprocal,
              "reciprocal_nonself_unordered_pairs": reciprocal // 2,
              "quantile_probabilities": probabilities,
              "outgoing_hit_count_quantiles": np.quantile(outgoing, probabilities).tolist(),
              "normalized_score_quantiles": np.quantile(scores, probabilities).tolist() if len(scores) else None,
              "species_directions": []}
    for source in np.unique(species):
        query_mask = species[queries] == source
        for target in np.unique(species):
            selected = query_mask & (species[targets] == target)
            report["species_directions"].append({"query_species": int(source), "target_species": int(target),
                "directed_hits": int(selected.sum()), "queries_with_hits": len(np.unique(queries[selected])),
                "query_species_genes": int(np.count_nonzero(species == source))})
    return report, codes


def overlap(first, second, n):
    result = {}
    for label, a, b in (("all", first, second),
                        ("nonself", first[first // n != first % n], second[second // n != second % n])):
        common = len(np.intersect1d(a, b, assume_unique=True))
        union = len(a) + len(b) - common
        result[label] = {"intersection": common, "first_only": len(a) - common,
                         "second_only": len(b) - common, "union": union,
                         "jaccard": common / union if union else None,
                         "fraction_of_first_recovered": common / len(a) if len(a) else None,
                         "fraction_of_second_recovered": common / len(b) if len(b) else None}
    return result


def assemble(root, output):
    root, output = root.resolve(), output.resolve()
    if output.exists():
        raise FileExistsError(output)
    accounting = subprocess.check_output(["sacct", "-j", "21292", "--parsable2",
                                         "--format=JobIDRaw,State,ExitCode,Elapsed"], text=True)
    scheduler = require_completed_job(accounting, 21292)
    results = root / "benchmark_tools/results"
    prepared = read_frozen(results / "orthobench_factorial_prepared_20260916.json", PREPARED_HASH)
    plan = read_frozen(results / "ob_sequence_search_prepared_20260916.json", MANIFEST_SHA)
    verify_plan(plan)
    conversion_path = root / "benchmarks/results/ob_sequence_numeric_v1/manifest.json"
    conversion_record = file_provenance(conversion_path)
    converted = json.loads(conversion_path.read_text())
    if (converted["status"] != "numeric_checkpoints_verified" or converted["accuracy_evaluated"] is not False
            or converted["scheduler"]["JobIDRaw"] != "21291"
            or set(converted["variants"]) != {"all_hits", "top100"}):
        raise ValueError("Unadmitted conversion or wrong variants")
    expected_source = root / "benchmarks/work/publication_ob_sequence_conversion_v1/benchmark_tools/convert_sequence_search_control.py"
    if converted["source"] != file_provenance(expected_source):
        raise ValueError("Conversion executor differs")
    verify_file(Path(converted["execution"]["path"]), converted["execution"])
    verify_file(Path(prepared["cache"]["path"]), prepared["cache"])
    names, species, q, t, s, evidence = load_replay_input(pickle_path=Path(prepared["cache"]["path"]))
    if evidence != prepared["cache"]:
        raise ValueError("HMM cache changed during loading")
    metadata = json.loads(Path(plan["gene_metadata"]["path"]).read_text())
    if names != sorted(metadata) or len(names) != 251378:
        raise ValueError("Search gene universes differ")
    species_names = sorted({item["species"] for item in metadata.values()})
    ids = {name: i for i, name in enumerate(species_names)}
    expected_species = np.array([ids[metadata[name]["species"]] for name in names], dtype=np.int32)
    same_species_partition(species, expected_species)
    summary, hmm_codes = summarize(q, t, s, expected_species)
    del q, t, s
    reports, intersections, admissions = {"hmm": summary}, {}, {}
    all_codes = None
    for label in ("all_hits", "top100"):
        variant = converted["variants"][label]
        if variant["cap"] != (None if label == "all_hits" else 100):
            raise ValueError("Unexpected reporting cap")
        other_names, other_species, q, t, s, admission = load_replay_input(
            checkpoint=Path(variant["checkpoint"]), checkpoint_sha256=variant["manifest"]["sha256"])
        if list(other_names) != names:
            raise ValueError("Numeric checkpoint gene order differs")
        if not np.array_equal(other_species, expected_species):
            raise ValueError("Numeric checkpoint species encoding differs")
        reports[label], codes = summarize(q, t, s, expected_species)
        intersections["hmm_vs_" + label] = overlap(hmm_codes, codes, len(names))
        admissions[label] = admission
        if label == "all_hits":
            all_codes = codes
        else:
            intersections["all_hits_vs_top100"] = overlap(all_codes, codes, len(names))
            if intersections["all_hits_vs_top100"]["all"]["second_only"]:
                raise ValueError("Top100 contains hits absent from all-hits control")
        del q, t, s
    verify_plan(plan)
    verify_file(Path(prepared["cache"]["path"]), prepared["cache"])
    verify_file(conversion_path, conversion_record)
    report = {"schema_version": 1, "accuracy_evaluated": False, "scheduler": scheduler,
              "conversion": conversion_record, "hmm_cache": evidence, "numeric_admissions": admissions,
              "species_labels": species_names, "searches": reports, "overlaps": intersections,
              "source": file_provenance(Path(__file__)), "numpy_version": np.__version__,
              "limitations": ["Hit-set overlap is not sensitivity against biological ground truth.",
                  "Equal E-value cutoffs and normalization formulas do not equate score calibration or sensitivity.",
                  "Normalized score quantiles are descriptive, not cross-engine significance thresholds.",
                  "Top100 is a post-search reporting diagnostic, not an HMM prefilter emulation.",
                  "No benchmark reference labels or accuracy outcomes were used."]}
    with output.open("x") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    assemble(args.root, args.output)


if __name__ == "__main__":
    main()
