"""Complete grouping replay, separate graph reader and provenance failures."""

import ast
from collections import Counter
import copy
import json
from pathlib import Path

import pytest

from benchmark_tools import readback_native_qfo_candidate_groups as reader
from benchmark_tools import trace_native_qfo_candidate_groups as primary
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_diagnose_candidate_trace_variation import row

ROOT = Path(__file__).resolve().parents[2]


def example():
    initial = {frozenset(g) for g in (("s",), ("t",), ("x",), ("a", "b"), ("u",))}
    final = {frozenset(("a", "b", "s", "t", "x")), frozenset(("u",))}
    trace = [row(), row("t", source_cluster=2),
             row("x", target=["a", "b", "s", "t"], source_cluster=3, iteration=1)]
    changes = [("a", "s", "FN", "TP"), ("s", "t", "not_scored", "FP"),
               ("s", "x", "not_scored", "FP"), ("t", "x", "not_scored", "FP")]
    return initial, final, trace, changes


def keyed(groups):
    return {min(g): g for g in groups}


def test_two_algorithms_reconstruct_all_groups_and_localizations():
    initial, final, trace, changes = example()
    ledger, summary, genes = primary.localize(initial, final, trace, changes)
    assert (ledger, summary) == reader.localize(keyed(initial), keyed(final), trace, changes)
    assert genes == 6
    assert [(r[7], r[8], r[9]) for r in ledger] == [
        (0, 0, "direct_cross_endpoint"), (0, "", "transitive_union"),
        (1, 2, "direct_cross_endpoint"), (1, 2, "direct_cross_endpoint")]
    assert sum(s["pairs"] for s in summary) == 4
    rebuilt, histories = reader.traverse(keyed(initial), trace, {"a", "s", "t", "x"})
    assert rebuilt == keyed(final)
    assert histories[0]["x"] != histories[0]["a"]


@pytest.mark.parametrize("which", ["primary", "reader"])
@pytest.mark.parametrize("mode", ["partial", "cycle", "label", "iteration", "duplicate_event",
    "overlap", "missing_gene", "collision", "unknown_pair", "unordered_pair", "duplicate_pair",
    "unsupported_state", "initial_cogroup", "no_changes"])
def test_localization_rejects_inconsistent_inputs(which, mode):
    initial, final, trace, changes = example()
    if mode == "partial":
        trace[2]["target_genes"] = ["a", "b", "s"]; trace[2]["target_size"] = 3
    elif mode == "cycle":
        trace.insert(2, row("s", target=["t"], target_size=1, target_cluster=2))
    elif mode == "label":
        trace[1]["source_cluster"] = 1
    elif mode == "iteration":
        trace[2]["iteration"] = True
    elif mode == "duplicate_event":
        trace.insert(0, copy.deepcopy(trace[0]))
    elif mode == "overlap":
        trace[0]["source_genes"] = ["a"]
    elif mode == "missing_gene":
        final.remove(frozenset(("u",)))
    elif mode == "collision":
        initial.update((frozenset(("sp|z|A",)), frozenset(("tr|z|B",))))
        final.update((frozenset(("sp|z|A",)), frozenset(("tr|z|B",))))
    elif mode == "unknown_pair":
        changes[0] = ("a", "unknown", "FN", "TP")
    elif mode == "unordered_pair":
        changes[0] = ("s", "a", "FN", "TP")
    elif mode == "duplicate_pair":
        changes.append(changes[0])
    elif mode == "unsupported_state":
        changes[0] = ("a", "s", "TP", "FN")
    elif mode == "initial_cogroup":
        changes[0] = ("a", "b", "FN", "TP")
    elif mode == "no_changes":
        changes = []
    fn = primary.localize if which == "primary" else reader.localize
    with pytest.raises(ValueError):
        fn(initial if which == "primary" else keyed(initial), final if which == "primary" else keyed(final), trace, changes)


@pytest.mark.parametrize("content", ["", "\n", "a a\n", "a b\nb c\n"])
def test_independent_partition_parser_rejects_bad_membership(tmp_path, content):
    path = tmp_path / "groups.txt"; path.write_text(content)
    with pytest.raises(ValueError):
        reader.partitions(path)


def test_accession_normalization_and_gene_keys_are_distinct():
    initial = {frozenset(("sp|a|A", "tr|b|B")), frozenset(("sp|s|S",))}
    final = {frozenset().union(*initial)}
    trace = [row("sp|s|S", target=["sp|a|A", "tr|b|B"])]
    changes = [("a", "s", "FN", "TP")]
    actual = primary.localize(initial, final, trace, changes)
    assert actual[:2] == reader.localize(keyed(initial), keyed(final), trace, changes)
    assert actual[0][0][4:7] == ["sp|a|A", "sp|s|S", "sp|a|A"]


@pytest.fixture
def retained(tmp_path):
    def write(name, value):
        path = tmp_path / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(value) if isinstance(value, (dict, list)) else value)
        return record(path)

    def build(mode=None):
        _, _, trace, changes = example()
        baseline = write("baseline/orthohmm_edges_clustered.txt", "s\nt\nx\na b\nu\n")
        partition_text = "a b s t x\nu\n" if mode != "partition" else "a b s t x\n"
        candidate = write("candidate/orthohmm_edges_clustered.txt", partition_text)
        checkpoint = write("candidate/phylogeny_candidate_superfamilies.txt",
                           partition_text if mode != "checkpoint" else "a b s t\nx\nu\n")
        if mode == "endpoint":
            trace[2]["target_genes"] = ["a", "b", "s"]; trace[2]["target_size"] = 3
        trace_ref = write("candidate/phylogeny_candidate_merges.json", trace)
        profile = dict(profile="satellite_v2", seed_families=5, candidate_families=2, merges=3)
        if mode == "seed_count": profile["seed_families"] = 4
        if mode == "merge_count": profile["merges"] = 4
        if mode == "profile": profile["profile"] = "tuned"
        metrics = write("candidate/metrics.json", dict(status="complete", metadata=dict(phylogeny_candidate_profile=profile)))
        methods = []
        for index, cell, participant, partition in ((6, "p0_c0_r0", "baseline", baseline), (8, "p0_c1_r0", "candidate", candidate)):
            job = 22435 if index == 6 else 22444
            output = dict(schema="native_factorial_output_review_v1", status="native_outputs_validated",
                native_outputs_validated=True, index=index, cell=cell, job_id=job, input_genes=6,
                checked_files=[partition] + ([trace_ref, metrics, checkpoint] if index == 8 else []))
            terminal = dict(schema="native_factorial_terminal_review_v1", status="native_success",
                native_outputs_validated=True, index=index, cell=cell, job_id=job)
            conversion = dict(schema="full_native_factorial_qfo_conversion_v1",
                status="full_native_factorial_qfo_pairs_prepared_unscored", conversion_kind="group",
                native_index=index, cell=cell, native_job_id=job, native_input=partition)
            admission = dict(schema="full_native_factorial_qfo_admission_v1",
                status="full_native_factorial_qfo_assessment_admitted", accuracy_admitted=True, publication_ready=False,
                native_index=index, cell=cell, native_job_id=job, participant=participant)
            if index == 8:
                if mode == "output_schema": output["schema"] = "unreviewed"
                if mode == "output_status": output["status"] = "unvalidated"
                if mode == "unbound_partition": output["checked_files"].remove(partition)
                if mode == "unbound_trace": output["checked_files"].remove(trace_ref)
                if mode == "ambiguous_trace":
                    output["checked_files"].append(write("other/phylogeny_candidate_merges.json", trace))
                if mode == "identical_inventory_duplicate": output["checked_files"].append(trace_ref)
                if mode == "input_genes": output["input_genes"] = 5
                if mode == "terminal": terminal["status"] = "native_failed"
                if mode == "conversion_kind": conversion["conversion_kind"] = "pairwise"
                if mode == "job": conversion["native_job_id"] = 0
                if mode == "admission_schema": admission["schema"] = "unadmitted"
                if mode == "admission_status": admission["status"] = "pending"
                if mode == "accuracy": admission["accuracy_admitted"] = False
                if mode == "cell": admission["cell"] = "p1_c1_r0"
            output_ref = write(participant + "/outputs_or_failure.json", output)
            terminal["reviews"] = dict(outputs_or_failure=output_ref)
            terminal_ref = write(participant + "/review.json", terminal)
            conversion["terminal_review"] = terminal_ref
            conversion_ref = write(participant + "/conversion.json", conversion)
            admission["pairs_manifest"] = conversion_ref
            admission["checked_records"] = [conversion_ref] if mode != "unbound_conversion" or index != 8 else []
            methods.append(dict(participant=participant, admission=write(participant + "/admission.json", admission)))
        rows = [("a", "b", "TP", "TP"), *changes]
        if mode == "duplicate_transition": rows.append(rows[0])
        if mode == "unsupported_transition": rows[1] = ("a", "s", "TP", "FN")
        text = "protein_left\tprotein_right\tblock_left\tblock_right\tp0_c0_r0\tp0_c1_r0\n"
        text += "".join("\t".join([a, b, "block1", "block2", l, r]) + "\n" for a, b, l, r in rows)
        pairs_ref = write("transitions.tsv", text)
        counts = Counter((r[2], r[3]) for r in rows)
        decomposition = dict(schema="native_qfo_candidate_vgnc_decomposition_v1", status="new_candidate_scored_rows_decomposed",
            source=record(ROOT / "benchmark_tools/export_native_qfo_candidate_vgnc.py"),
            baseline=methods[0], candidate=methods[1], transition_table=pairs_ref, union_scored_pairs=len(rows),
            transition_counts=[dict(baseline=a, candidate=b, pairs=n) for (a, b), n in sorted(counts.items())],
            uncertainty_admitted=False, new_scoring_or_admission=False, publication_ready=False,
            failed_r1_timing_remains_ineligible=True)
        if mode == "pair_count": decomposition["union_scored_pairs"] -= 1
        if mode == "scope": decomposition["publication_ready"] = True
        decomp_ref = write("decomposition.json", decomposition)
        readback = dict(schema="native_qfo_candidate_vgnc_readback_v1", status="candidate_decomposition_independently_verified",
            source=record(ROOT / "benchmark_tools/readback_native_qfo_candidate_vgnc.py"), report=decomp_ref,
            uncertainty_admitted=False, new_scoring_or_admission=False, publication_ready=False,
            failed_r1_timing_remains_ineligible=True)
        if mode == "readback_link": readback["report"] = pairs_ref
        return decomp_ref, write("decomposition_readback.json", readback)
    return write, build


def test_full_provenance_export_and_independent_reader(retained, tmp_path):
    _, build = retained
    refs = build()
    result = primary.execute(*refs, tmp_path / "result")
    path = tmp_path / "result/report.json"
    checked = reader.review(path, record(path)["sha256"])
    assert checked["changed_pairs"] == result["changed_pairs"] == 4
    assert checked["localized_summary"] == result["localized_summary"]
    assert checked["whole_candidate_partition_reconstructed"]
    assert not checked["primary_or_accepted_replay_imported"]
    for key in ("accuracy_rescored", "native_inference_reexecuted", "uncertainty_admitted", "publication_ready"):
        assert result[key] is checked[key] is False
    assert checked["failed_r1_timing_remains_ineligible"]
    with pytest.raises(ValueError, match="fresh direct"):
        primary.execute(*refs, tmp_path / "result")


def test_identical_inventory_duplicates_are_not_conflicts(retained, tmp_path):
    _, build = retained
    primary.execute(*build("identical_inventory_duplicate"), tmp_path / "result")
    path = tmp_path / "result/report.json"
    assert reader.review(path, record(path)["sha256"])["changed_pairs"] == 4


@pytest.mark.parametrize("mode", ["partition", "checkpoint", "endpoint", "seed_count", "merge_count", "profile",
    "output_schema", "output_status", "unbound_partition", "unbound_trace", "ambiguous_trace", "input_genes",
    "terminal", "conversion_kind", "job", "admission_schema", "admission_status", "accuracy", "cell",
    "unbound_conversion", "duplicate_transition", "unsupported_transition", "pair_count", "scope", "readback_link"])
def test_provenance_or_scientific_mismatch_retains_failure(retained, tmp_path, mode):
    _, build = retained
    refs = build(mode)
    with pytest.raises(ValueError):
        primary.execute(*refs, tmp_path / "result")
    failure = json.loads((tmp_path / "result/failure.json").read_text())
    assert failure["error_type"] == "ValueError" and failure["error"]
    assert failure["automatic_retry"] is failure["publication_ready"] is False
    assert not (tmp_path / "result/report.json").exists()


@pytest.mark.parametrize("mode", ["scope", "count", "summary", "unbound_source", "missing_method", "checkpoint", "ledger"])
def test_independent_reader_rejects_modified_output(retained, tmp_path, mode):
    write, build = retained
    report = primary.execute(*build(), tmp_path / "result")
    if mode == "scope": report["uncertainty_admitted"] = True
    elif mode == "count": report["genes"] -= 1
    elif mode == "summary": report["localized_summary"] = []
    elif mode == "unbound_source": report["checked_records"].remove(report["source"])
    elif mode == "missing_method": report["methods"].pop()
    elif mode == "checkpoint": report["candidate_checkpoint"] = report["methods"][0]["partition"]
    elif mode == "ledger":
        path = Path(report["ledger"]["path"])
        path.write_text("\n".join(path.read_text().splitlines()[:-1]) + "\n")
        report["ledger"] = record(path)
    ref = write("tampered_report.json", report)
    with pytest.raises(ValueError):
        reader.review(ref["path"], ref["sha256"])


def test_changed_input_hash_is_rejected(retained, tmp_path):
    _, build = retained
    refs = build()
    Path(refs[0]["path"]).write_text("{}")
    with pytest.raises(ValueError):
        primary.execute(*refs, tmp_path / "result")
    assert (tmp_path / "result/failure.json").is_file()


def test_reader_does_not_import_primary_scientific_kernels():
    tree = ast.parse(Path(reader.__file__).read_text())
    imports = [n.module for n in ast.walk(tree) if isinstance(n, ast.ImportFrom)]
    assert "benchmark_tools" in imports
    names = [alias.name for n in ast.walk(tree) if isinstance(n, (ast.Import, ast.ImportFrom)) for alias in n.names]
    assert "readback_native_qfo_vgnc_blocks" in names
    assert not set(names) & {"trace_native_qfo_candidate_groups", "replay_native_candidate_trace", "partition", "numpy"}
