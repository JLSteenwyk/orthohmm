"""Run the frozen four-rule diagnostic without changing trees or method defaults."""

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import dendropy

from benchmark_tools.audit_wgd_results import recalculate
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.reconstruct_reconciliation_trace import apply_logged_constraints
from benchmark_tools.run_publication_pipeline import save
from benchmark_tools.run_wgd_application import pinned
from benchmark_tools.trace_wgd_cases import run as trace, partition_index
from benchmark_tools.validate_scaling_outputs import input_universe


RULES = ("species_overlap", "supported_children", "confidence", "mapped_event")
PROTOCOL_SHA = "1659a135e2584c1b1b947447c5517d69ea52dcee1d804f8d9bd30b1569a27764"
TRACE_SHA = "5b346d99d1ac063fdf169ecf76fde69390d876073e9c27d43968469936ac07cc"
SOURCE_SHA = "216d97608e6dede8960f3b06fe7036bd23da8fae5677df88583ba9369ec3c1bf"


def partition(groups):
    canonical = sorted(tuple(sorted(g)) for g in groups)
    partition_index({str(i): list(g) for i, g in enumerate(canonical)})
    return canonical


def biological_row(row):
    return {k: v for k, v in row.items() if k != "anchor_groups"}


def frozen_module(path):
    if record(path)["sha256"] != SOURCE_SHA:
        raise ValueError("Frozen phylogeny source changed")
    spec = importlib.util.spec_from_file_location("_wgd_fixed_tree_phylogeny", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def run(repo):
    results = repo / "benchmark_tools/results"
    protocol = record(results / "WGD_FIXED_TREE_RULE_PROTOCOL_20260928.md")
    if protocol["sha256"] != PROTOCOL_SHA:
        raise ValueError("Diagnostic protocol changed")
    trace_path = results / "biological_wgd_case_trace_20260917.json"
    trace_ref = record(trace_path)
    if trace_ref["sha256"] != TRACE_SHA:
        raise ValueError("Retained WGD trace changed")
    retained = json.loads(trace_path.read_text())
    if trace(repo) != retained:
        raise ValueError("Retained native trace did not reproduce")
    refs = [protocol, trace_ref, *retained["checked_artifacts"]]
    paths = {Path(r["path"]).name: Path(r["path"]) for r in retained["checked_artifacts"]}
    module = frozen_module(paths["phylogeny.py"])
    audit = pinned(retained["application_audit"])
    application = pinned(audit["report"])
    prepared = pinned(application["input_manifest"])
    reference = pinned(prepared["reference"])
    refs.extend([retained["application_audit"], audit["report"], application["input_manifest"], prepared["reference"]])
    owners, _ = input_universe({**prepared, "proteomes": 4})
    keys = [tuple(case["orf_pair"]) for case in retained["cases"]]
    cohort = [next(r for r in prepared["cohort_pairs"] if tuple(r["orf_pair"]) == key) for key in keys]
    expected_rows = {tuple(c["orf_pair"]): c["stages"]["final"] for c in retained["cases"]}
    pair_evidence, arms = {}, []
    for rule in RULES:
        families, final_groups = {}, {}
        for family, original in retained["families"].items():
            genes = set(original["candidate_members"])
            if original["reconciled"]:
                gene_tree = dendropy.Tree.get(path=str(paths[family + ".rooted.nwk"]),
                    schema="newick", preserve_underscores=True, rooting="force-rooted")
                species_tree = dendropy.Tree.get(path=str(paths["species_tree.rooted.nwk"]),
                    schema="newick", preserve_underscores=True, rooting="force-rooted")
                value = module.reconcile_gene_tree(gene_tree, species_tree, owners,
                    family_id=family, root_duplication_rule=rule, pair_orthology_rule="positive_paralogy")
                pairs = (value.ortholog_pairs, value.paralog_pairs, value.ortholog_pair_confidence)
                if rule == RULES[0]:
                    pair_evidence[family] = pairs
                elif pairs != pair_evidence[family]:
                    raise ValueError("Root-rule intervention changed pair evidence")
                groups = value.root_groups
                high = {tuple(sorted((a, b))) for a, b, confidence in value.ortholog_pair_confidence
                        if confidence == "high"}
            else:
                groups, high = [genes], set()
                if original["constraints"]:
                    raise ValueError("Unexpected constraint in bypass family")
            reconstruction = dict(root_groups=list(map(set, groups)), high_confidence_pairs=high)
            events = [(e["event_index"], e) for e in original["constraints"]]
            final, constraints = apply_logged_constraints(reconstruction, events, genes)
            if set().union(*final) != genes:
                raise ValueError("Alternative lost candidate genes")
            if rule == RULES[0] and (partition(groups) != partition(original["preconstraint_groups"])
                                    or partition(final) != partition(original["final_groups"])):
                raise ValueError("Fixed-tree baseline groups do not reproduce native output")
            families[family] = dict(preconstraint_groups=partition(groups), final_groups=partition(final),
                                    constraints=constraints)
            for i, group in enumerate(final):
                final_groups[f"{family}:{i}"] = sorted(group)
        scored = recalculate(cohort, final_groups, reference, owners)
        if rule == RULES[0] and any(biological_row(r) != biological_row(expected_rows[tuple(r["orf_pair"])])
                                  for r in scored):
            raise ValueError("Baseline biological scoring fields differ")
        assignments = partition_index(final_groups)
        destinations = []
        for case in retained["cases"]:
            for lost in case["homolog_destinations"]:
                if lost["stage"] == "root_lineage_split":
                    gene, anchors = lost["gene"], case["orf_pair"]
                    destinations.append(dict(gene=gene, orf_pair=anchors, group=assignments[gene],
                        with_anchor=[a for a in anchors if assignments[a] == assignments[gene]]))
        arms.append(dict(rule=rule, families=families, cases=scored, focal_homologs=destinations))
    for ref in refs:
        check(ref)
    return dict(status="fixed_tree_rule_diagnostic_complete", protocol=protocol, trace=trace_ref,
        inputs=refs, source=record(__file__), dendropy_version=dendropy.__version__,
        arms=arms, all_pair_evidence_unchanged=True,
        pair_evidence_sha256=hashlib.sha256(json.dumps(pair_evidence, sort_keys=True).encode()).hexdigest(),
        defaults_changed=False, independent_validation=False, publication_ready=False,
        limitations=["Post hoc selected-case mechanism diagnostic, not population accuracy or a tuning recommendation.",
            "Tree topology, rooting, candidate membership, pair rules and constraints held fixed.",
            "Rule effects conditional on these inputs do not prove ancestral-copy truth or correct topology.",
            "No new uncertainty interval or global method ranking."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    save(args.output, run(args.repo.resolve()))
