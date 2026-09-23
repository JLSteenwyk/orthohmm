import pytest
import hashlib
import json

from benchmark_tools.audit_unisave_fragment_source import choose_version, inspect_entry
from benchmark_tools import verify_swiss_historical_fragments as verifier
from benchmark_tools.verify_swiss_historical_fragments import selected_history, family_bins
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def row(version, sv=1, first="12-Aug-2020", last="12-Aug-2020"):
    return dict(accession="P12345", sequenceVersion=sv, entryVersion=version,
                firstReleaseDate=first, lastReleaseDate=last)


@pytest.mark.parametrize("rows,sv", [([row(1)], 1),
    ([row(2, first="02-Dec-2020", last="02-Dec-2020"), row(1)], 1),
    ([row(3, 2, "10-Feb-2021", "10-Feb-2021"), row(2, 2, "02-Dec-2020", "02-Dec-2020"), row(1)], 2)])
def test_independent_selection_agrees_on_valid_histories(rows, sv):
    history = dict(results=rows)
    assert selected_history(history, "P12345", sv) == choose_version(history, "P12345", sv)


@pytest.mark.parametrize("rows", [[], [row(1), row(1)], [row(1), row(2)],
    [row(1, first="02-Dec-2020")], [row(1, first="17-Jun-2020", last="17-Jun-2020")],
    [dict(row(1), accession="OTHER")]])
def test_bad_history_rejected(rows):
    with pytest.raises(ValueError):
        selected_history(dict(results=rows), "P12345", 1)


def annotation(positive=False, feature=False, later=False):
    return dict(fragment_flag=positive, incomplete_sequence_features=[{}] if feature else [],
                selection_class="later_sequence_version" if later else "baseline_release")


def test_bins_preserve_missing_and_later_sensitivity():
    families = {"positive_missing": ["a", "b"], "unflagged": ["c"],
                "feature": ["d"], "later": ["e"], "missing": ["f", "c"]}
    annotations = dict(a=annotation(True), b=None, c=annotation(),
                       d=annotation(feature=True), e=annotation(True, later=True), f=None)
    bins = family_bins(families, annotations)
    assert bins == dict(annotation_positive=["feature", "later", "positive_missing"],
                        all_matched_unflagged=["unflagged"], missing_without_positive=["missing"])
    baseline = family_bins(families, annotations, True)
    assert baseline["annotation_positive"] == ["feature", "positive_missing"]
    assert baseline["missing_without_positive"] == ["later", "missing"]


@pytest.mark.parametrize("families,annotations", [({"f": []}, {}), ({"f": ["a", "a"]}, {"a": None}),
    ({"f": ["a"]}, {}), ({"f": ["a"]}, {"a": None, "extra": None})])
def test_incomplete_or_ambiguous_bins_rejected(families, annotations):
    with pytest.raises(ValueError):
        family_bins(families, annotations)


@pytest.mark.parametrize("problem", [None, "altered_raw", "wrong_gene_digest", "incomplete", "wrong_count"])
def test_full_panel_validation(tmp_path, monkeypatch, problem):
    directory = tmp_path / "benchmarks/work/swiss_historical_fragment_panel_20260923"
    directory.mkdir(parents=True)
    results = tmp_path / "benchmark_tools/results"
    results.mkdir(parents=True)
    genes, descriptors, text = {}, {}, []
    for i in range(563):
        gene = f"P{i:05d}"
        description = f"sp|{gene}|TEST Protein OX=9606 SV=1"
        text.append(f">{description}\nAAA\n")
        genes[gene] = dict(accession=gene, sequence_version=1, taxid="9606", length=3,
            sequence_sha256=hashlib.sha256(b"AAA").hexdigest())
        descriptors[gene] = dict(description=description)
    families = {f"family{i}": list(genes)[i::18] for i in range(18)}
    fastas = []
    for i in range(78):
        path = tmp_path / f"input{i}.fa"
        path.write_text("".join(text) if i == 0 else ">unused|OUTSIDE|TEST\nAAA\n")
        fastas.append(record(path))
    inventory_path = results / "corrected_swiss_sequence_strata_20260918.json"
    inventory_path.write_text(json.dumps(dict(genes=descriptors, family_memberships=families, fasta_inputs=fastas)))
    inputs = [record(inventory_path)]
    monkeypatch.setattr(verifier, "INVENTORY_SHA", inputs[0]["sha256"])
    for name in ("PROTOCOL_SHA", "HELPER_SHA", "COLLECTOR_SHA"):
        path = tmp_path / name
        path.write_text(name)
        inputs.append(record(path))
        monkeypatch.setattr(verifier, name, inputs[-1]["sha256"])
    entries = {g: dict(status="missing", error_type="HTTPError", error="Synthetic missing source") for g in genes}
    gene = "P00000"
    source = directory / gene
    source.mkdir()
    selected = dict(row(2), accession=gene)
    (source / "history.json").write_text(json.dumps(dict(results=[selected])))
    raw = ("ID   TEST_HUMAN              Reviewed;         3 AA.\nAC   P00000;\n"
        "DT   01-JAN-2000, integrated into UniProtKB/Swiss-Prot.\n"
        "DT   01-JAN-2000, sequence version 1.\nDT   12-AUG-2020, entry version 2.\n"
        "DE   RecName: Full=Test;\nDE   Flags: Fragment;\nOS   Homo sapiens (Human).\n"
        "OC   Eukaryota.\nOX   NCBI_TaxID=9606;\n"
        "SQ   SEQUENCE   3 AA;  100 MW;  0000000000000000 CRC64;\n     AAA\n//\n")
    (source / "entry.txt").write_text(raw)
    audit = inspect_entry(raw, gene, 2, 1, genes[gene]["sequence_sha256"], "9606")
    audit.update(selection=selected, selection_class="baseline_release", outcome_analysis_performed=False,
        source=inputs[2], records=[record(source / n) for n in ("history.json", "entry.txt")])
    (source / "audit.json").write_text(json.dumps(audit))
    entries[gene] = dict(status="sequence_matched", selection_class="baseline_release", audit=record(source / "audit.json"))
    report = dict(status="collection_complete_pending_independent_validation", job_id="22116",
        prediction_statistics_evaluated=False, annotation_panel_admitted=False, publication_ready=False,
        preflight=dict(genes=genes, families=families, records=inputs+fastas), entries=entries, matched=1, missing=562)
    if problem == "altered_raw":
        (source / "entry.txt").write_text(raw.replace("AAA", "BBB"))
    elif problem == "wrong_gene_digest":
        genes[gene]["sequence_sha256"] = "0"*64
    elif problem == "incomplete":
        report["status"] = "collecting"
    elif problem == "wrong_count":
        report["missing"] = 0
    (directory / "status.json").write_text(json.dumps(report))
    monkeypatch.setattr(verifier.subprocess, "check_output", lambda *a, **k:
        "JobIDRaw|State|ExitCode|Elapsed\n22116|COMPLETED|0:0|00:01:00\n")
    output = tmp_path / "verified.json"
    if problem:
        with pytest.raises(ValueError):
            verifier.verify(tmp_path, output)
        assert not output.exists()
    else:
        result = verifier.verify(tmp_path, output)
        assert result["matched"] == 1 and result["missing"] == 562
        assert result["strata"]["annotation_positive"] == ["family0"]
        assert len(result["strata"]["missing_without_positive"]) == 17
        assert result["publication_ready"] is False
        with pytest.raises(FileExistsError):
            verifier.verify(tmp_path, output)
