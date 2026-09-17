"""Recalculate WGD diagnostics from admitted native groups without the scorer."""

import argparse
from collections import Counter
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.assemble_wgd_application import ADMISSIONS, artifact, unchanged
from benchmark_tools.read_wgd_native_groups import read_orthohmm, read_orthofinder_root_ids, read_species_table
from benchmark_tools.report_ygob_validation import read_checkpoint
from benchmark_tools.run_wgd_application import pinned
from benchmark_tools.snapshot_orthohmm_input_order import record
from benchmark_tools.validate_scaling_outputs import input_universe

METHODS = ("orthohmm_high_sensitivity", "orthohmm_satellite_v2", "orthofinder_full", "sonicparanoid")
DIAGNOSTIC = "orthofinder_mcl_checkpoint"
ENDPOINTS = ("separation_rate", "supported_separation_rate", "mean_non_scer_coverage")


def recalculate(cohort, groups, reference, owners):
    # Deliberately separate arithmetic from score_wgd_application and its helpers.
    assignments = {}
    for group, members in groups.items():
        for gene in members:
            if gene in assignments or gene not in owners:
                raise ValueError("Invalid native partition")
            assignments[gene] = group
    curated = {gene: pillar for pillar, genes in reference.items() for gene in genes}
    if len(curated) != sum(map(len, reference.values())):
        raise ValueError("Ambiguous reference membership")
    output = []
    for source in cohort:
        group_ids = [assignments.get(g) for g in source["orf_pair"]]
        members = [set(groups[g]) if g is not None else set() for g in group_ids]
        union = members[0] | members[1]
        if not source["split_eligible"]:
            state = "input_excluded"
        elif any(g is None for g in group_ids):
            state = "incomplete_assignment"
        elif group_ids[0] == group_ids[1]:
            state = "merged"
        else:
            state = "separated"
        row = dict(source, assignment_state=state, anchor_groups=group_ids,
                   anchor_group_sizes=[len(g) for g in members], union_size=len(union),
                   separation_rate=None if state == "input_excluded" else int(state == "separated"))
        for key in ("supported_separation_rate", "mean_non_scer_coverage", "homolog_support_by_anchor",
                    "coverage_numerator", "coverage_denominator", "foreign_pillar_members", "unmapped_members",
                    "pillar_native_group_count", "unassigned_pillar_members"):
            row[key] = None
        if source["reference_eligible"]:
            pillar = source["reference_pillar"]
            homologs = [g for g in reference[pillar] if owners[g] != "Scerevisiae"]
            support = [sum(g in group for g in homologs) for group in members]
            recovered = sum(g in union for g in homologs)
            row.update(homolog_support_by_anchor=support, coverage_numerator=recovered,
                       coverage_denominator=len(homologs),
                       mean_non_scer_coverage=recovered / len(homologs) if homologs else None,
                       supported_separation_rate=int(state == "separated" and support[0] > 0 and support[1] > 0),
                       foreign_pillar_members=sorted(g for g in union if g in curated and curated[g] != pillar),
                       unmapped_members=sorted(g for g in union if g not in curated),
                       pillar_native_group_count=len({assignments[g] for g in reference[pillar] if g in assignments}),
                       unassigned_pillar_members=sorted(g for g in reference[pillar] if g not in assignments))
        output.append(row)
    return output


def verify_entry(entry, recalculated):
    expected = {tuple(r["orf_pair"]): r for r in recalculated}
    observed = {tuple(r["orf_pair"]): r for r in entry["rows"]}
    if len(observed) != len(entry["rows"]) or expected != observed:
        raise ValueError("Native-membership rescore differs from reported rows")
    for label, summary in [(None, entry["summary"]), *entry["strata"].items()]:
        rows = [r for r in recalculated if label is None or r["experimental_class"] == label]
        eligible = [r for r in rows if r["reference_eligible"]]
        check = {"pairs": len(rows), "input_eligible": sum(r["split_eligible"] for r in rows),
                 "reference_eligible": len(eligible), "assignment_states": dict(Counter(r["assignment_state"] for r in rows)),
                 "endpoints": {}}
        for endpoint in ENDPOINTS:
            values = [r[endpoint] for r in eligible if r[endpoint] is not None]
            check["endpoints"][endpoint] = {"pairs": len(values), "mean": sum(values) / len(values) if values else None}
        if check != summary:
            raise ValueError("Recalculated summary or stratum differs")


def audit(repo, report_path):
    report = json.loads(report_path.read_text())
    for item in report["sources"] + report["admissions"] + [report["protocol"], report["contrast_protocol"]]:
        unchanged(item)
    prepared = pinned(report["input_manifest"])
    owners, species = input_universe({**prepared, "proteomes": 4})
    reference = pinned(prepared["reference"])
    cohort = prepared["cohort_pairs"]
    groups, admissions = {}, []
    for method, (filename, digest) in zip(METHODS, ADMISSIONS):
        path = repo / "benchmark_tools/results" / filename
        native = pinned({"path": str(path), "sha256": digest})
        admissions.append(record(path))
        if native["method"] != method:
            raise ValueError("Method admission mismatch")
        for item in native["native"]["checked_files"]:
            unchanged(item)
        if method == METHODS[0]:
            groups[method] = read_orthohmm(artifact(native, "orthohmm_orthogroups.txt"), "named_groups", owners)
        elif method == METHODS[1]:
            groups[method] = read_orthohmm(artifact(native, "orthohmm_root_hogs.tsv"), "root_hogs", owners)
        elif method == METHODS[2]:
            ids = artifact(native, "SequenceIDs.txt")
            groups[method] = read_orthofinder_root_ids(artifact(native, "N0.ids.tsv"), ids, owners, {s: s for s in species})
            checkpoints = {Path(r["path"]).name for r in native["native"]["checked_files"]
                           if Path(r["path"]).name.startswith("clusters_OrthoFinder_I") and r["path"].endswith("_id_pairs.txt")}
            if len(checkpoints) != 1:
                raise ValueError("Ambiguous MCL checkpoint")
            groups[DIAGNOSTIC] = read_checkpoint(artifact(native, checkpoints.pop()), ids, owners)
        else:
            columns = {Path(r["path"]).name: Path(r["path"]).stem for r in prepared["inputs"]}
            groups[method] = read_species_table(artifact(native, "ortholog_groups.tsv"), "sonicparanoid", owners, columns)
    if admissions != report["admissions"]:
        raise ValueError("Unexpected admission provenance")
    rows = {}
    for method, partition in groups.items():
        rows[method] = recalculate(cohort, partition, reference, owners)
        verify_entry(report["methods"][method], rows[method])
    examples = report["prespecified_examples"]
    if len(examples) != len(prepared["prespecified_examples"]):
        raise ValueError("Changed example count")
    for source, example in zip(prepared["prespecified_examples"], examples):
        if {k: v for k, v in example.items() if k != "methods"} != source:
            raise ValueError("Changed prespecified example")
        for method in groups:
            expected = next(r for r in rows[method] if r["orf_pair"] == source["orf_pair"])
            if example["methods"][method] != expected:
                raise ValueError("Changed example result")
    population = sorted((r for r in cohort if r["reference_eligible"]), key=lambda r: r["reference_pillar"])
    if len(population) != 231 or len({r["reference_pillar"] for r in population}) != 231:
        raise ValueError("Direct resampling audit requires the frozen 231 distinct pillars")
    order = [tuple(r["orf_pair"]) for r in population]
    draws = np.random.Generator(np.random.PCG64(20260920)).integers(231, size=(20000, 231))
    comparisons = report["uncertainty"]["comparisons"]
    expected_contrasts = [(METHODS[1], METHODS[0]), (METHODS[1], METHODS[2]),
                          (METHODS[1], METHODS[3]), (METHODS[0], METHODS[2])]
    if [(c["first"], c["second"], c["endpoint"]) for c in comparisons] != [
            (a, b, e) for a, b in expected_contrasts for e in ENDPOINTS]:
        raise ValueError("Changed contrast inventory")
    for comparison in comparisons:
        columns = []
        for method in (comparison["first"], comparison["second"]):
            lookup = {tuple(r["orf_pair"]): r for r in rows[method]}
            columns.append(np.array([lookup[k][comparison["endpoint"]] for k in order]))
        delta = columns[0] - columns[1]
        samples = delta[draws].mean(axis=1) * 100
        checks = {"difference_pp": delta.mean() * 100,
                  "nominal95_pp": np.quantile(samples, [.025, .975]),
                  "bonferroni12_pp": np.quantile(samples, [.05 / 24, 1 - .05 / 24]),
                  "wins": (delta > 0).sum(), "losses": (delta < 0).sum(), "ties": (delta == 0).sum(),
                  "pairs": 231, "pillars": 231, "undefined_replicates": 0}
        for key, value in checks.items():
            if not np.allclose(comparison[key], value, rtol=0, atol=1e-10):
                raise ValueError("Direct bootstrap differs: " + key)
    return {"status": "native_membership_and_interval_audit_passed", "report": record(report_path),
            "auditor": record(__file__), "methods": len(groups), "pairs_per_method": len(cohort),
            "contrasts": len(comparisons), "prespecified_examples": len(examples),
            "admissions": admissions,
            "limitations": ["Native format readers and identity checks are shared with the report pipeline.",
                            "Separate arithmetic and direct resampling do not validate biological reference truth.",
                            "Direct resampling check is restricted to this frozen one-pair-per-pillar population."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.repo.resolve(), args.report.resolve())
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
