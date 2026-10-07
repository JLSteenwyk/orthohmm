"""Proven original-protein joins must preserve all scored IDs and paths."""

import ast
from collections import Counter
import copy
import csv
import gzip
import json
from pathlib import Path
import sqlite3

import pytest

from benchmark_tools import diagnose_native_qfo_candidate_group_ids as diagnosis
from benchmark_tools import join_native_qfo_candidate_aliases as primary
from benchmark_tools import readback_native_qfo_candidate_aliases as reader
from benchmark_tools import trace_native_qfo_candidate_groups as previous
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_native_qfo_candidate_groups import example, keyed, retained

ROOT = Path(__file__).resolve().parents[2]


def bridge_fixture():
    return (dict(mapping={"x": 4, "z": 4}, Goff=[0, 100], species=["TEST"]),
            [[1, 4, "x", "TEST"], [2, 4, "z", "TEST"]], {"x": "x"}, ["z"])


def test_original_qfo_identity_and_db_rows_prove_bridge():
    args = bridge_fixture()
    assert primary.bridge(*args) == reader.bridges(*args) == [dict(scored_accession="z", native_accession="x",
        native_gene="x", prot_nr=4, species="TEST", original_proteome_rows=[
            dict(rowid=1, prot_nr=4, uniprot_id="x", species="TEST"), dict(rowid=2, prot_nr=4, uniprot_id="z", species="TEST")])]


@pytest.mark.parametrize("which", ["primary", "reader"])
@pytest.mark.parametrize("mode", ["no_alias", "invalid_number", "invalid_offsets", "bad_species", "missing_alias_row",
    "missing_native_row", "ambiguous_native", "different_protein", "multiple_targets", "no_missing", "present_alias", "bool_native"])
def test_bridge_rejects_guessing_or_conflicting_proof(which, mode):
    mapping, rows, native, missing = copy.deepcopy(bridge_fixture())
    if mode == "no_alias": del mapping["mapping"]["z"]
    elif mode == "invalid_number": mapping["mapping"]["z"] = True
    elif mode == "invalid_offsets": mapping["Goff"] = [100, 0]
    elif mode == "bad_species": rows[1][3] = "OTHER"
    elif mode == "missing_alias_row": rows.pop()
    elif mode == "missing_native_row": rows.pop(0)
    elif mode == "ambiguous_native":
        mapping["mapping"]["x2"] = 4; native["x2"] = "x2"; rows.append([3, 4, "x2", "TEST"])
    elif mode == "different_protein": mapping["mapping"]["x"] = 5
    elif mode == "multiple_targets":
        mapping["mapping"]["z2"] = 4; missing.append("z2"); rows.append([3, 4, "z2", "TEST"])
    elif mode == "no_missing": missing = []
    elif mode == "present_alias": native["z"] = "z"
    elif mode == "bool_native":
        mapping["mapping"] = {"x": True, "z": 1}
        for r in rows: r[1] = 1
    fn = primary.bridge if which == "primary" else reader.bridges
    with pytest.raises(ValueError):
        fn(mapping, rows, native, missing)


def test_join_preserves_original_pair_orientation_and_all_paths():
    initial, final, trace, changes = example()
    changes[-1] = ("t", "z", "not_scored", "FP")
    proof = primary.bridge(*bridge_fixture())
    rows, summary, genes = primary.joined_localization(initial, final, trace, changes, proof)
    assert (rows, summary) == reader.localize(keyed(initial), keyed(final), trace, changes, proof)
    assert rows[-1][:6] == ["t", "z", "t", "x", "not_scored", "FP"]
    assert genes == 6 and len(rows) == 4
    # A canonical ID that moves to the other side must not reverse raw orientation.
    changes[-1] = ("q", "t", "not_scored", "FP")
    proof = [dict(scored_accession="q", native_accession="x")]
    rows, summary, _ = primary.joined_localization(initial, final, trace, changes, proof)
    assert (rows, summary) == reader.localize(keyed(initial), keyed(final), trace, changes, proof)
    assert rows[-1][:8] == ["q", "t", "x", "t", "not_scored", "FP", "x", "t"]


@pytest.mark.parametrize("which", ["primary", "reader"])
@pytest.mark.parametrize("mode", ["self", "alias_collision", "missing", "duplicate_bridge"])
def test_join_rejects_collapsed_or_unmapped_pairs(which, mode):
    initial, final, trace, changes = example()
    proof = [dict(scored_accession="z", native_accession="x")]
    if mode == "self": changes.append(("x", "z", "not_scored", "FP"))
    elif mode == "alias_collision": changes.append(("t", "z", "not_scored", "FP"))
    elif mode == "missing": changes[-1] = ("t", "unknown", "not_scored", "FP")
    elif mode == "duplicate_bridge": proof *= 2
    fn = primary.joined_localization if which == "primary" else reader.localize
    with pytest.raises(ValueError):
        fn(initial if which == "primary" else keyed(initial), final if which == "primary" else keyed(final), trace, changes, proof)


@pytest.fixture
def bound_aliases(retained, tmp_path):
    write, build = retained
    def prepare(mode=None):
        refs = build("identifier")
        decomposition = json.loads(Path(refs[0]["path"]).read_text())
        readback = json.loads(Path(refs[1]["path"]).read_text())
        mapping_path = tmp_path / "mapping.json.gz"
        mapping = dict(mapping={"a": 1, "b": 2, "s": 3, "x": 4, "t": 5, "u": 6, "z": 4}, Goff=[0, 100], species=["TEST"])
        if mode == "bridge": mapping["mapping"]["z"] = 99
        with gzip.open(mapping_path, "wt") as stream: json.dump(mapping, stream)
        mapping_ref = record(mapping_path)
        db_path = tmp_path / "candidate.db"
        with sqlite3.connect(db_path) as connection:
            connection.execute("CREATE TABLE proteomes (prot_nr INTEGER, uniprot_id TEXT, species TEXT)")
            connection.executemany("INSERT INTO proteomes VALUES (?, ?, ?)", [(4, "x", "TEST"), (4, "z", "TEST")])
        database_ref = record(db_path)
        execution = dict(status="process_succeeded_pending_independent_admission", native_index=8, cell="p0_c1_r0",
                         exit_code=0, outputs=[database_ref] if mode != "unbound_db" else [])
        execution_ref = write("execution.json", execution)
        for name in ("baseline", "candidate"):
            admission = json.loads(Path(decomposition[name]["admission"]["path"]).read_text())
            conversion = json.loads(Path(admission["pairs_manifest"]["path"]).read_text())
            conversion["mapping"] = mapping_ref
            conversion["checked_records"] = [mapping_ref] if mode != "unbound_mapping" else []
            conversion_ref = write(name + "/conversion.json", conversion)
            admission["pairs_manifest"] = conversion_ref; admission["checked_records"] = [conversion_ref]
            if name == "candidate":
                admission["execution_report"] = execution_ref
                admission["checked_records"].extend([execution_ref, database_ref,
                    record(ROOT / "qfo_benchmark/benchmark-webservice/map_relations.py"),
                    record(ROOT / "qfo_benchmark/benchmark-webservice/vgnc_benchmark.py")])
                if mode == "unbound_scorer": admission["checked_records"].pop()
                decomposition[name]["execution"] = execution_ref; decomposition[name]["database"] = database_ref
            decomposition[name]["admission"] = write(name + "/admission.json", admission)
        decomp_ref = write("decomposition.json", decomposition)
        readback["report"] = decomp_ref
        rb_ref = write("decomposition_readback.json", readback)
        with pytest.raises(ValueError, match="unmapped changed pair"):
            previous.execute(decomp_ref, rb_ref, tmp_path / "original")
        failure_ref = record(tmp_path / "original/failure.json")
        diagnosis.execute(failure_ref, tmp_path / "prior_diagnosis")
        return record(tmp_path / "prior_diagnosis/report.json")
    return write, prepare


def test_complete_original_provenance_and_both_algorithms(bound_aliases, tmp_path):
    _, prepare = bound_aliases
    prior = prepare()
    result = primary.execute(prior, tmp_path / "new_join")
    path = tmp_path / "new_join/report.json"
    check = reader.review(path, record(path)["sha256"])
    assert check["bridges"] == result["bridges"] and check["changed_pairs"] == 4
    assert check["localized_summary"] == result["localized_summary"]
    assert check["complete_pair_localization"] and check["original_scored_identifiers_preserved"]
    assert not check["exporter_or_primary_replay_imported"]
    assert not check["original_failed_export_retried"] and not check["publication_ready"]
    assert not (tmp_path / "original/report.json").exists()
    with pytest.raises(ValueError, match="fresh direct"):
        primary.execute(prior, tmp_path / "new_join")


@pytest.mark.parametrize("mode", ["bridge", "unbound_db", "unbound_mapping", "unbound_scorer"])
def test_provenance_failure_retained_before_new_join(bound_aliases, tmp_path, mode):
    _, prepare = bound_aliases
    with pytest.raises(ValueError): primary.execute(prepare(mode), tmp_path / "new_join")
    result = json.loads((tmp_path / "new_join/failure.json").read_text())
    assert result["automatic_retry"] is result["original_failed_export_retried"] is False
    assert not (tmp_path / "new_join/report.json").exists()


@pytest.mark.parametrize("mode", ["bridge", "sql_rows", "scope", "count", "summary", "source", "kernel_source", "ledger", "original_ids"])
def test_independent_reader_rejects_changed_report(bound_aliases, tmp_path, mode):
    write, prepare = bound_aliases
    result = primary.execute(prepare(), tmp_path / "new_join")
    if mode == "bridge": result["bridges"][0]["native_accession"] = "t"
    elif mode == "sql_rows": result["original_database_rows"] = []
    elif mode == "scope": result["accuracy_rescored"] = True
    elif mode == "count": result["changed_pairs"] -= 1
    elif mode == "summary": result["localized_summary"] = []
    elif mode == "source": result["checked_records"].remove(result["source"])
    elif mode == "kernel_source":
        result["checked_records"].remove(record(ROOT / "benchmark_tools/replay_native_candidate_trace.py"))
    elif mode == "original_ids": result["original_scored_identifiers_preserved"] = False
    elif mode == "ledger":
        path = Path(result["ledger"]["path"])
        path.write_text("\n".join(path.read_text().splitlines()[:-1]) + "\n")
        result["ledger"] = record(path)
    ref = write("changed_report.json", result)
    with pytest.raises(ValueError): reader.review(ref["path"], ref["sha256"])


def test_reader_import_boundary():
    tree = ast.parse(Path(reader.__file__).read_text())
    names = [a.name for n in ast.walk(tree) if isinstance(n, (ast.Import, ast.ImportFrom)) for a in n.names]
    assert "readback_native_qfo_candidate_groups" in names
    assert not set(names) & {"join_native_qfo_candidate_aliases", "trace_native_qfo_candidate_groups", "replay_native_candidate_trace", "numpy"}


def test_actual_full_ledger_readback_and_failure_identities_without_rerun():
    directory = ROOT / "benchmark_tools/results"
    path = directory / "native_qfo_candidate_alias_group_20261006_v1/report.json"
    if not path.is_file(): pytest.skip("Selected alias result not installed")
    assert record(path)["sha256"] == "3b3ab8f1132e70362e1450ac990caf7a9dad536e528bc54d6af39db16183cd53"
    report = json.loads(path.read_text())
    rb_path = directory / "native_qfo_candidate_alias_group_readback_20261006_v1.json"
    assert record(rb_path)["sha256"] == "4afc4aa6853a69027b9dc710602309f44923315b1e68dbae4cd1e61866f4b92d"
    readback = json.loads(rb_path.read_text())
    assert readback["report"] == record(path) and readback["bridges"] == report["bridges"]
    assert readback["localized_summary"] == report["localized_summary"]
    assert (report["genes"], report["baseline_groups"], report["candidate_groups"], report["accepted_merges"],
            report["changed_pairs"]) == (984137, 394328, 353638, 40690, 2295)
    assert [(b["scored_accession"], b["native_accession"], b["prot_nr"], b["species"]) for b in report["bridges"]] == [
        ("Q17QN5_BOVIN", "Q17QN5", 594577, "BOVIN"), ("Q1RMT5_BOVIN", "Q1RMT5", 594851, "BOVIN")]
    assert record(report["original_failure"]["path"])["sha256"] == "bfad1eefd8f084adc6ab66d5da385c9ce8a1f17440948c59330cf1a761efb695"
    assert record(report["diagnosis"]["path"])["sha256"] == "3bc6646e1dc7cb6488ceb6aaaf763e3a4e3cdcbd41268f2066ef2d8008f5be74"
    assert record(report["ledger"]["path"])["sha256"] == "3f0598f7e00cfc0942146604ab1863ad9f03f2d5017efd5244c7e76a29d43642"
    for ref in report["checked_records"]: assert record(ref["path"]) == ref
    with Path(report["ledger"]["path"]).open(newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    assert len(rows) == 2295 and len({(r["protein_left"], r["protein_right"]) for r in rows}) == 2295
    observed = Counter((r["candidate_state"], int(r["first_connected_iteration"]), r["connection_path"]) for r in rows)
    assert [dict(candidate_state=s, iteration=i, connection_path=p, pairs=n) for (s, i, p), n in sorted(observed.items())] == report["localized_summary"]
    assert sum(r["candidate_state"] == "TP" for r in rows) == 162
    assert sum(r["candidate_state"] == "FP" for r in rows) == 2133
    assert sum(r["candidate_state"] == "FP" and r["connection_path"] == "direct_cross_endpoint" for r in rows) == 2033
    assert sum(any("_BOVIN" in r[k] for k in ("protein_left", "protein_right")) for r in rows) == 8
    for key in ("whole_candidate_partition_reconstructed", "complete_pair_localization", "original_scored_identifiers_preserved"):
        assert report[key] is readback[key] is True
    for key in ("original_failed_export_retried", "accuracy_rescored", "native_inference_reexecuted",
                "new_scoring_or_admission", "uncertainty_admitted", "publication_ready"):
        assert report[key] is readback[key] is False


def test_manuscript_result_table_matches_exact_retained_summary():
    directory = ROOT / "benchmark_tools/results"
    report = json.loads((directory / "native_qfo_candidate_alias_group_20261006_v1/report.json").read_text())
    counts = {(r["candidate_state"], r["iteration"], r["connection_path"]): r["pairs"] for r in report["localized_summary"]}
    for name in ("PUBLICATION_MANUSCRIPT_DRAFT_20260916.md", "NATIVE_QFO_CANDIDATE_ALIAS_GROUP_RESULT_20261006.md"):
        text = (directory / name).read_text()
        for iteration in (0, 1):
            for path, label in (("direct_cross_endpoint", "Direct cross-endpoint"), ("transitive_union", "Transitive union")):
                tp, fp = counts["TP", iteration, path], counts["FP", iteration, path]
                assert f"| {iteration} | {label} | {tp:,} | {fp:,} |" in text
        assert "| Total | All paths | 162 | 2,133 |" in text
    manuscript = (directory / "PUBLICATION_MANUSCRIPT_DRAFT_20260916.md").read_text()
    assert "The initial diagnosis does not admit" in manuscript
    assert "ruling out an all-transitive explanation" in manuscript
    assert "not suffix stripping or an independent biological validation" in manuscript
    claims = (directory / "PUBLICATION_CLAIMS_20260916.md").read_text()
    assert "localizes ALL2,295changed pair paths" in claims
    assert "not a causal support" in claims
