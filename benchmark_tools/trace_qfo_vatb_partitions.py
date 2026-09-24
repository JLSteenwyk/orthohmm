"""Post-hoc VATB membership trace; co-grouping is not native orthology scoring."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_simulation_methods import read_frozen

COUNTS_SHA = "1d6bd4f68a9d6728e0180bdde31f55ca0da7b4fc22e09bcc797dcb8db6cec628"
BASELINE_SHA = "d8385c50426e690afd6d32f3c5302e678de6977c9841451d013201b0f75b564a"


def trace(path, accessions):
    wanted = set(accessions)
    if not wanted or len(wanted) != len(accessions):
        raise ValueError("Invalid reference accession inventory")
    seen, groups = set(), []
    with path.open() as stream:
        for number, line in enumerate(stream, 1):
            tokens = line.split()
            selected = []
            for token in tokens:
                parts = token.split("|")
                accession = parts[1] if len(parts) == 3 and parts[0] in {"sp", "tr"} else token
                if accession in wanted:
                    if accession in seen:
                        raise ValueError("Duplicate or ambiguous reference accession")
                    seen.add(accession)
                    selected.append(dict(accession=accession, gene=token))
            if selected:
                groups.append(dict(line=number, total_genes=len(tokens), reference_genes=selected))
    if seen != wanted:
        raise ValueError("Reference members missing from partition")
    return dict(reference_genes=len(seen), represented_groups=len(groups),
        singleton_reference_genes=sum(len(g["reference_genes"]) for g in groups if g["total_genes"] == 1),
        co_grouped_reference_pairs=sum(len(g["reference_genes"]) * (len(g["reference_genes"]) - 1) // 2 for g in groups),
        groups=groups)


def run(root, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    counts_path = root / "benchmark_tools/results/qfo_parameter_uncertainty_partial_20260923.json"
    counts = read_frozen(counts_path, COUNTS_SHA)
    families = {}
    for arm in counts["reconstructed_counts"]["arms"]:
        if arm["arm"] in {"control", "cpm_low"}:
            families[arm["arm"]] = next(row for row in arm["families"] if row["family"] == "VATB")
    accessions = families["control"]["represented_genes"]
    if len(accessions) != 28 or accessions != families["cpm_low"]["represented_genes"]:
        raise ValueError("Wrong or inconsistent VATB reference inventory")
    seed_path = root / "benchmarks/work/qfo_cpm_variant_admission_22082_0.json"
    candidate_path = root / "benchmarks/work/qfo_cpm_candidates_admission_22086_0.json"
    baseline_path = root / "benchmarks/results/qfo_corrected_factorial_v1/manifest.json"
    records = [record(__file__), record(counts_path), record(seed_path), record(candidate_path), record(baseline_path)]
    seed = read_frozen(seed_path, records[2]["sha256"])
    candidate = read_frozen(candidate_path, records[3]["sha256"])
    baseline = read_frozen(baseline_path, BASELINE_SHA)
    if (seed["status"] != "cpm_variant_replay_admitted_unscored" or seed["arm"] != "cpm_low"
            or candidate["status"] != "cpm_candidates_admitted_unscored" or candidate["arm"] != "cpm_low"):
        raise ValueError("Wrong low-CPM stage admission")
    stages = [("cpm_low", row["label"], row["output"]) for row in seed["coverage"]]
    stages.append(("cpm_low", "candidate", candidate["candidate_arm"]["candidate_partition"]))
    stages.extend(("control", label, baseline["candidate_arms"]["p1_c1"][key]) for label, key in (
        ("strict_profiles_refined", "seed_partition"), ("candidate", "candidate_partition")))
    rows = []
    for arm, stage, item in stages:
        check(item)
        rows.append(dict(arm=arm, stage=stage, partition=item, **trace(Path(item["path"]), accessions)))
        records.append(item)
    for item in records:
        check(item)
    result = dict(status="vatb_partition_membership_traced", family="VATB", selection="post_hoc_largest_low_cpm_family_decline",
        reference_counts=families, stages=rows, checked_records=records, publication_ready=False,
        limitations=["Post-hoc development-exposed error localization, not an independent confirmation or causal intervention.",
                     "Co-grouped pair counts include all reference pairs and are not native ortholog predictions or true positives.",
                     "Reference-accession uniqueness is checked here; full-partition integrity relies on retained stage admissions.",
                     "No fresh search, clustering, refinement, phylogeny or scoring was executed."])
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.root.resolve(), args.output.absolute())
