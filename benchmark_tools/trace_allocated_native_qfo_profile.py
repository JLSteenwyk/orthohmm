"""Localize every observed profile-contrast pair change without rerunning inference."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

from benchmark_tools import trace_native_qfo_swiss_reconciliation as topology
from benchmark_tools import trace_native_qfo_swiss_transitions as transitions
from benchmark_tools.export_native_factorial_progress import load, require
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

CELLS = ("p0_c0_r1", "p1_c0_r1")


def member_hash(genes):
    return hashlib.sha256(("\n".join(sorted(genes)) + "\n").encode()).hexdigest()


def describe(a, b, selected, members, nodes, predictions):
    left, right = selected[a], selected[b]
    family_a, family_b = left["source_family"], right["source_family"]
    same = family_a == family_b
    pair = tuple(sorted((left["gene"], right["gene"])))
    node = topology.pair_lca(nodes[family_a], *pair) if same else None
    predicted = pair in predictions
    require(predicted is (same and node["pair_event"] != "duplication"),
        "Native pair disagrees with observed grouping/positive-paralogy rule")
    return dict(gene_a=left["gene"], gene_b=right["gene"], source_family_a=family_a,
        source_family_b=family_b, same_source_family=same, predicted=predicted,
        root_hog_a=left["root_hog"], root_hog_b=right["root_hog"],
        candidate_members_a_sha256=member_hash(members[family_a]),
        candidate_members_b_sha256=member_hash(members[family_b]),
        candidate_size_a=len(members[family_a]), candidate_size_b=len(members[family_b]),
        lca=None if node is None else {k: node[k] for k in ("node_id", "event", "pair_event",
            "event_confidence", "species_overlap_count", "mapping_conflict", "branch_support")})


def run(readback_path, readback_sha, output):
    require(not output.exists() and not output.is_symlink(), "Require fresh profile trace output")
    evidence = []
    readback, readback_ref = load(readback_path, readback_sha, evidence)
    require(readback["schema"] == "allocated_native_qfo_profile_swiss_rational_readback_v1"
        and readback["source"] == record(Path(__file__).with_name("readback_allocated_native_qfo_profile_swiss.py"))
        and readback["cells"] == list(CELLS) and readback["profile_pair_labels_matched"] == 10765
        and readback["families_checked"] == 18 and readback["prior_matched_contrasts_checked"] == 2
        and readback["prior_matched_contrasts_unchanged"] is True
        and all(readback[k] is False for k in ("new_accuracy_or_resource_admission",
            "scientific_timings_admitted", "independent_confirmation", "publication_ready")),
        "Require independently read-back profile contrast")
    evidence.extend([*readback["checked_inputs"], readback["source"]])
    binding, _ = load(readback["binding"]["path"], readback["binding"]["sha256"], evidence)
    require(binding["schema"] == "allocated_native_qfo_retained_swiss_uncertainty_binding_v1",
        "Changed interval binding route")
    count_rows = []
    for cell in CELLS:
        ref = binding["bound_cells"][cell]["count_audit"]
        audit, _ = load(ref["path"], ref["sha256"], evidence)
        selected = [r for r in audit["cells"] if r["cell"] == cell]
        require(len(selected) == 1 and binding["bound_cells"][cell]["status"] == "native_records_matched",
            "Missing matched native counts")
        count_rows.append(selected[0])
    families = binding["families"]
    labels = [transitions.read_labels(Path(r["raw_file"]["path"]), families) for r in count_rows]
    comparison, changes = transitions.compare_labels(*labels, families, *count_rows)
    require(changes, "No changed profile pairs to localize")
    targets = {p for row in changes for p in row[1:3]}
    states, observed = [], []
    for index, cell in zip((7, 10), count_rows):
        admission, _ = load(cell["admission"]["path"], cell["admission"]["sha256"], evidence)
        require(admission["accuracy_admitted"] is True and admission["native_index"] == index
            and admission["cell"] == cell["cell"] and admission["native_job_id"] == cell["native_job_id"],
            "Changed native admission identity")
        stage = admission["conversion"]
        review_ref = stage["scientific_recovery"] if index == 7 else stage["terminal_review"]
        review, _ = load(review_ref["path"], review_ref["sha256"], evidence)
        output_ref = review["outputs"] if index == 7 else review["reviews"]["outputs_or_failure"]
        reviewed, _ = load(output_ref["path"], output_ref["sha256"], evidence)
        prediction = stage["native_input"]
        require(reviewed["native_outputs_validated"] is True and reviewed["cell"] == cell["cell"]
            and prediction in reviewed["checked_files"] and stage["conversion_kind"] == "native",
            "Unadmitted native profile prediction")
        if index == 7:
            require(admission["resources"] is None and admission["eligible_for_timing_comparison"] is False,
                "Profile trace repairs failed reference timing")
        directory = Path(prediction["path"]).parent
        records = {}
        for name in ("orthohmm_root_hogs.tsv", "provenance_manifest.json"):
            refs = [r for r in reviewed["checked_files"] if r["path"] == str(directory / name)]
            require(len(refs) == 1, "Missing admitted profile " + name)
            records[name] = refs[0]
            evidence.append(refs[0])
            check(refs[0])
        manifest = json.loads(Path(records["provenance_manifest.json"]["path"]).read_text())
        require(manifest["membership_reconciliation"] is None and manifest["pair_orthology_rule"] == "positive_paralogy",
            "Different family/pair reconciliation rule")
        selected, members, reconstruction = topology.source_groups(
            Path(records["orthohmm_root_hogs.tsv"]["path"]), targets, manifest["input_cluster_sha256"])
        require(reconstruction["input_genes"] == reviewed["input_genes"]
            and reconstruction["root_hogs"] == reviewed["phylogeny"]["root_hogs"], "Changed complete root universe")
        # Node tables were not inventoried by the original native-output validator.
        node_ref = record(directory / "orthohmm_reconciliation_nodes.tsv")
        nodes, node_rows = topology.read_nodes(Path(node_ref["path"]), members)
        evidence.extend([prediction, node_ref, review["source"], reviewed["source"]])
        queries = {tuple(sorted((selected[row[1]]["gene"], selected[row[2]]["gene"]))) for row in changes}
        found, pair_rows = topology.selected_predictions(Path(prediction["path"]), queries)
        require(pair_rows == reviewed["phylogeny"]["native_pair_rows"], "Changed admitted native pair inventory")
        points = [describe(row[1], row[2], selected, members, nodes, found) for row in changes]
        require(all(p["predicted"] is (row[3 if index == 7 else 4] in ("TP", "FP"))
            for p, row in zip(points, changes)), "Scored relation differs from actual native prediction")
        states.append(points)
        observed.append(dict(cell=cell["cell"], admission=cell["admission"], native_output_review=output_ref,
            reconstruction=reconstruction, root_hogs=records["orthohmm_root_hogs.tsv"],
            manifest=records["provenance_manifest.json"], observed_node_annotations=node_ref,
            annotation_rows_read=node_rows, selected_source_families=len(members), native_pair_rows_checked=pair_rows,
            species_tree_sha256=manifest["species_tree_sha256"]))
    rows, summary = [], Counter()
    for changed, before, after in zip(changes, *states):
        require((before["gene_a"], before["gene_b"]) == (after["gene_a"], after["gene_b"]),
            "Changed native target gene identity")
        category = "candidate_separation" if not after["same_source_family"] else "positive_paralogy_exclusion"
        require(before["predicted"] and not after["predicted"], "Unexpected non-removal transition")
        summary[category] += 1
        rows.append(dict(family=changed[0], protein_a=changed[1], protein_b=changed[2],
            before_label=changed[3], after_label=changed[4], before=before, after=after,
            localization=category, candidate_members_a_unchanged=before["candidate_members_a_sha256"] == after["candidate_members_a_sha256"],
            candidate_members_b_unchanged=before["candidate_members_b_sha256"] == after["candidate_members_b_sha256"]))
    helpers = [record(topology.__file__), record(transitions.__file__)]
    source = record(__file__)
    for ref in [*evidence, *helpers, source]:
        check(ref)
    result = dict(schema="allocated_native_qfo_profile_pair_localization_v1", source=source,
        readback=readback_ref, evidence=evidence, helpers=helpers, comparison=comparison, states=observed,
        changed_pairs=rows, changed_pairs_traced=len(rows), summary=dict(summary),
        species_tree_bytes_identical=observed[0]["species_tree_sha256"] == observed[1]["species_tree_sha256"],
        node_annotations_previously_inventoried=False, new_scoring_or_admission=False,
        uncertainty_admitted=False, scientific_timings_admitted=False, independent_confirmation=False,
        publication_ready=False, limitations=[
            "Every changed assessed SwissTrees profile pair, not all input or submitted relations.",
            "Observed grouping/LCA rule explains software inclusion, not true duplication history or tree correctness.",
            "Root/pair/manifest outputs admitted; node annotations newly observed, topology/rule checked without rerun.",
            "Does not identify individual profile scores/edges as a cause or isolate species-tree changes.",
            "Group hashes refer to complete candidate members, not earlier search/profile histories.",
            "Development exposure and failed reference timing remain; no bootstrap, scoring or inference rerun."])
    with output.open("x") as stream:
        json.dump(result, stream, sort_keys=True, indent=2, allow_nan=False)
        stream.write("\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("readback", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--readback-sha256", required=True)
    args = parser.parse_args()
    result = run(args.readback, args.readback_sha256, args.output)
    print(json.dumps({k: result[k] for k in ("changed_pairs_traced", "summary", "species_tree_bytes_identical")}))


if __name__ == "__main__":
    main()
