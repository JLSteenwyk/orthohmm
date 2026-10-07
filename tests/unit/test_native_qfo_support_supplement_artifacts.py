"""Check the actual companion and render bindings without rerunning inference."""

import csv
import hashlib
import json
from pathlib import Path
import subprocess

import fitz
from PIL import Image

from benchmark_tools import prepare_native_qfo_support_supplement as generator

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "benchmark_tools/results/native_qfo_support_supplement_20261006_v1"
PRINT = OUT.parent / "native_qfo_support_supplement_print_20261006_v1"
REVIEW = OUT.parent / "native_qfo_support_supplement_pdf_review_20261006_v1"


def load(path):
    return json.loads(path.read_text())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def table(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def test_actual_generation_outputs_and_copied_evidence_exact():
    receipt = load(OUT / "generation.json")
    assert len(receipt["outputs"]) == 14
    for ref in receipt["outputs"] + receipt["evidence"]:
        path = Path(ref["path"])
        assert path.stat().st_size == ref["bytes"] and digest(path) == ref["sha256"]
    for name, ref in receipt["copied_outputs"].items():
        assert Path(ref["path"]).read_bytes() == (OUT / name).read_bytes()
    assert (OUT / "report.json").read_bytes() == Path(receipt["evidence"][0]["path"]).read_bytes()
    assert (OUT / "readback.json").read_bytes() == Path(receipt["evidence"][1]["path"]).read_bytes()


def test_committed_generator_produces_same_manuscript_from_verified_snapshot():
    generation = load(OUT / "generation.json")
    report, reader = load(OUT / "report.json"), load(OUT / "readback.json")
    groups = generator.snapshot(report, reader, generation["evidence"][0])
    assert generator.manuscript(report, groups) == (OUT / "supplement.md").read_text()
    committed = subprocess.run(["git", "show", "f721207e:benchmark_tools/prepare_native_qfo_support_supplement.py"],
        cwd=ROOT, check=True, capture_output=True).stdout
    assert committed == Path(generator.__file__).read_bytes()
    assert hashlib.sha256(committed).hexdigest() == generation["source"]["sha256"]


def test_complete_tables_match_actual_summary_fields():
    report = load(OUT / "report.json")
    cohorts = table(OUT / "event_cohorts.tsv")
    features = table(OUT / "feature_summaries.tsv")
    paths = table(OUT / "pair_paths.tsv")
    assert len(cohorts) == 12 and len(features) == 204 and len(paths) == 8
    expected_cohorts, expected_features = [], []
    for row in report["cohort_summaries"]:
        iteration = "all" if row["iteration"] is None else str(row["iteration"])
        expected_cohorts.append(dict(iteration=iteration, cohort=row["cohort"], events=str(row["events"])))
        for feature in generator.FEATURES:
            expected_features.append(dict(iteration=iteration, cohort=row["cohort"], events=str(row["events"]), feature=feature,
                **{k: "" if v is None else str(v) for k, v in row["features"][feature].items()}))
    assert cohorts == expected_cohorts and features == expected_features
    assert paths == [{k: str(v) for k, v in row.items()} for row in report["localized_summary"]]


def test_claims_do_not_upgrade_scientific_evidence():
    claims = load(OUT / "claims.json")
    assert len(claims["claims"]) == 4
    assert all(c["status"] == "independently_verified_descriptive" for c in claims["claims"])
    assert claims["not_established"] == ["biological correctness of unlabeled events", "direct protein-pair HMM hits", "calibrated confidence",
        "causal error mechanisms", "valid VGNC confidence intervals", "independent confirmation", "optimal thresholds", "OrthoFinder superiority", "publication readiness"]
    assert all((OUT / c["evidence"]).is_file() for c in claims["claims"])
    generation = load(OUT / "generation.json")
    assert generation["new_bootstrap_draws"] == 0 and generation["original_trace_or_partition_replayed"] is False
    assert generation["new_accuracy_or_resource_admission"] is False and generation["defaults_changed"] is False


def test_actual_figure_pixels_and_pdf_text_nonblank():
    with Image.open(OUT / "accepted_event_support.png") as image:
        assert image.size == (2560, 1760)
        colors = image.convert("RGB").getcolors(maxcolors=2560 * 1760)
        assert len(colors) > 100 and sum(n for n, color in colors if color != (255, 255, 255)) > 50000
    with fitz.open(OUT / "accepted_event_support.pdf") as doc:
        assert len(doc) == 1
        text = " ".join(doc[0].get_text().split())
        assert "not confidence intervals" in text and "Unlabeled does not mean correct" in text
        assert "40,241" in text and "33,699" in text and "6,542" in text


def test_actual_local_asset_bindings_and_three_page_pdf():
    render = load(OUT / "render.json")
    printing = load(PRINT / "print.json")
    review = load(REVIEW / "report.json")
    assert render["local_occurrences"] == render["unique_targets"] == 11 and render["citation_ids"] == []
    for ref in [render["html"], *render["sources"], *render["targets"], printing["pdf"], *review["rendered_pages"]]:
        assert digest(Path(ref["path"])) == ref["sha256"]
    assert printing["status"] == "verified_html_printed" and printing["returncode"] == 0
    assert review["page_count"] == printing["page_count"] == 3 and review["bounds_violations"] == []
    assert review["phrase_matches"] == {"Methods": [1], "Interpretation": [2], "Machine-readable": [3]}
    with fitz.open(PRINT / "document.pdf") as doc:
        all_text = " ".join(" ".join(page.get_text().split()) for page in doc)
        assert "40,690 accepted events" in all_text and "2,295 changed pair paths" in all_text
        assert "46" in all_text and "Failed R-on timing remains ineligible" in all_text
        assert "Claim-To-Evidence Checklist" in all_text


def test_scoped_manual_review_anchors_and_limits():
    manual = load(OUT.parent / "native_qfo_support_supplement_visual_review_20261006_v1.json")
    targets = dict(generation=OUT / "generation.json", manuscript=OUT / "supplement.md", claims=OUT / "claims.json",
        figure_pdf=OUT / "accepted_event_support.pdf", render=OUT / "render.json", printed_pdf=PRINT / "document.pdf",
        pdf_bounds_review=REVIEW / "report.json")
    assert all(manual[key + "_sha256"] == digest(path) for key, path in targets.items())
    assert manual["pdf_pages_viewed"] == [1, 2, 3] and manual["full_resolution_figure_viewed"] is True
    assert manual["all_companion_pdf_pages_viewed"] is True and manual["full_study_or_main_manuscript_visual_review"] is False
    assert manual["independent_scientific_reproduction"] is False and manual["publication_ready"] is False
