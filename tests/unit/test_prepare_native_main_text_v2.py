"""Check presentation arithmetic, preservation and refusal behavior, not new inference."""

import copy
import hashlib
import json
from pathlib import Path
import subprocess

import pytest

from benchmark_tools import prepare_native_main_text_v2 as current
from benchmark_tools.render_manuscript_review import citation_ids, local_assets


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / current.BASE


@pytest.fixture
def values():
    result = {}
    for key in ("original", "scores", "swiss", "vgnc", "aliases", "aliases_readback"):
        path = BASE / current.PINS[key][0]
        data = path.read_bytes()
        assert hashlib.sha256(data).hexdigest() == current.PINS[key][1]
        result[key] = data.decode("utf-8") if key == "original" else json.loads(data)
    return result


def compose(values):
    return current.compose(*(values[k] for k in (
        "original", "scores", "swiss", "vgnc", "aliases", "aliases_readback")))


def normalized(values):
    return " ".join(compose(values).split())


def test_all_eighteen_native_scores_and_coverage_are_from_admitted_values(values):
    text = compose(values)
    rows = values["scores"]["rows"][:3]
    for endpoint in current.ENDPOINTS:
        statistic = "F1" if endpoint in current.ENDPOINTS[:3] else (
            "Sample mean" if endpoint == "FAS" else "Similarity")
        expected = f"| {endpoint} | {statistic} | " + " | ".join(
            f'{row["scores"][endpoint]:.8f}' for row in rows) + " |"
        assert expected in text
    for row in rows:
        assert f'{row["secondary_mean"]:.8f}' in text
        assert f'{row["submitted_pairs"]:,}' in text
        assert f'{100 * row["relation_coverage"]:.4f}%' in text
    for phrase in ("Three of seven", "four remain unavailable", "failed 22437 timing",
                   "not official QfO F1", "not the selected-default all-tool",
                   "initial sensitive HMM search on", "not a complete fresh factorial"):
        assert phrase in normalized(values)


def test_all_retained_contrast_metrics_and_zero_crossing_limits(values):
    text = normalized(values)
    matched = [c for c in values["swiss"]["contrasts"] if c["metrics"] is not None]
    assert len(matched) == 2
    for contrast in matched:
        for metric in contrast["metrics"].values():
            assert f'{100 * metric["difference"]:+.4f}' in text
            low, high = metric["bonferroni_percentile_ci"]
            assert f'[{100 * low:.4f}, {100 * high:.4f}]' in text
            assert "/".join(str(metric[k]) for k in (
                "family_wins", "family_ties", "family_losses")) in text
    assert "Both adjusted F1 intervals include zero" in text
    assert "other 12 remain unavailable" in text
    assert "no new draws" in text
    assert "neither equivalence nor independent confirmation" in text


def test_complete_path_table_matches_independent_group_summary(values):
    text = compose(values)
    summary = {(r["iteration"], r["connection_path"], r["candidate_state"]): r["pairs"]
               for r in values["aliases"]["localized_summary"]}
    for iteration in (0, 1):
        for path, label in (("direct_cross_endpoint", "Direct group endpoints"),
                            ("transitive_union", "Transitive union")):
            assert (f"| {iteration} | {label} | {summary[(iteration, path, 'TP')]:,} | "
                    f"{summary[(iteration, path, 'FP')]:,} |") in text
    assert "| Total | All localized paths | 162 | 2,133 |" in text
    for phrase in ("2,295 changed scored pairs", "without suffix guessing",
                   "all-transitive software-path explanation", "does not prove a direct protein-pair HMM hit",
                   "not validated independent uncertainty units", "not a true negative"):
        assert phrase in normalized(values)
    for number in values["vgnc"]["differences"].values():
        assert f"{100 * number:+.6f}" in text


def section(text, start, end):
    return current.bounded(text, start, end)[1]


def test_unrelated_science_and_functional_diagnostic_are_preserved_exactly(values):
    old, new = values["original"], compose(values)
    for start, end in (
        ("## Abstract\n", "### Shared-Host Resource Measurement\n"),
        ("### Shared-Host Resource Measurement\n", current.NATIVE_START),
        (current.NATIVE_END, "## Reproducibility And Availability\n"),
    ):
        if start == "## Abstract\n":
            old_part = section(old, start, end)
            assert section(new, start, end).startswith(old_part)
        else:
            assert section(old, start, end) == section(new, start, end)
    assert section(old, current.FUNCTIONAL_START, current.NATIVE_END) == section(
        new, current.FUNCTIONAL_START, current.NATIVE_END)
    assert section(old, "## Reproducibility And Availability\n", current.AVAILABILITY_START) == section(
        new, "## Reproducibility And Availability\n", current.AVAILABILITY_START)
    assert new.split("## References\n", 1)[1] == old.split("## References\n", 1)[1]


def test_existing_local_targets_and_eighteen_citations_remain_with_new_evidence(values):
    def parse(text):
        return json.loads(subprocess.check_output([
            "pandoc", "--from=markdown", "--to=json"], input=text, text=True))
    old, new = parse(values["original"]), parse(compose(values))
    main = BASE / "PUBLICATION_MAIN_TEXT_20261006_v2.md"
    _, old_targets = local_assets(old, main, ROOT)
    _, new_targets = local_assets(new, main, ROOT)
    assert old_targets.keys() <= new_targets.keys()
    assert citation_ids(old) == citation_ids(new) and len(citation_ids(new)) == 18
    for key in ("scores", "swiss", "swiss_readback", "vgnc", "vgnc_readback",
                "aliases", "aliases_readback", "ledger", "figure_pdf"):
        assert str(current.BASE / current.PINS[key][0]) in new_targets


@pytest.mark.parametrize("mutation", [
    "cohort", "unadmitted_score", "mean", "coverage", "timing", "draws", "contrasts",
    "family_counts", "interval", "path_total", "independent_paths", "incomplete", "partition",
])
def test_inconsistent_scientific_scope_or_summary_is_rejected(values, mutation):
    changed = copy.deepcopy(values)
    s, u, a, r = (changed[k] for k in ("scores", "swiss", "aliases", "aliases_readback"))
    if mutation == "cohort": s["rows"][2]["accuracy_admitted"] = False
    elif mutation == "unadmitted_score": s["rows"][3]["scores"]["FAS"] = 0
    elif mutation == "mean": s["rows"][0]["secondary_mean"] = 0
    elif mutation == "coverage": s["rows"][0]["relation_coverage"] = 1
    elif mutation == "timing": s["rows"][1]["timing_eligible"] = True
    elif mutation == "draws": u["new_bootstrap_draws"] = 1
    elif mutation == "contrasts": u["contrasts"].pop()
    elif mutation == "family_counts":
        next(c for c in u["contrasts"] if c["metrics"])["metrics"]["F1"]["family_wins"] += 1
    elif mutation == "interval":
        next(c for c in u["contrasts"] if c["metrics"])["metrics"]["F1"]["bonferroni_percentile_ci"] = [1, 2]
    elif mutation == "path_total":
        a["localized_summary"][0]["pairs"] += 1
        r["localized_summary"] = copy.deepcopy(a["localized_summary"])
    elif mutation == "independent_paths": r["localized_summary"][0]["pairs"] += 1
    elif mutation == "incomplete": a["complete_pair_localization"] = False
    elif mutation == "partition": a["candidate_groups"] += 1
    with pytest.raises(ValueError): compose(changed)


@pytest.mark.parametrize("start,end", [("missing", current.NATIVE_END),
                                      (current.NATIVE_START, "missing")])
def test_missing_or_duplicate_boundaries_rejected(values, start, end):
    with pytest.raises(ValueError): current.bounded(values["original"], start, end)
    with pytest.raises(ValueError):
        current.bounded(values["original"] + current.NATIVE_START,
                        current.NATIVE_START, current.NATIVE_END)


def copied_inputs(tmp_path):
    base = tmp_path / current.BASE
    base.mkdir(parents=True)
    for name, _ in current.PINS.values():
        destination = base / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes((BASE / name).read_bytes())
    return base


def test_actual_generation_is_reproducible_and_does_not_overwrite(tmp_path, values):
    base = copied_inputs(tmp_path)
    output, receipt = base / "new.md", base / "receipt.json"
    result = current.build(output, receipt, tmp_path)
    assert output.read_text() == compose(values)
    assert json.loads(receipt.read_text()) == result
    assert result["publication_ready"] is False and result["new_bootstrap_draws"] == 0
    assert result["original_evidence_modified"] is False
    assert result["unavailable_cells"] == 4 and result["unavailable_swiss_contrasts"] == 12
    for name, sha in current.PINS.values():
        assert hashlib.sha256((base / name).read_bytes()).hexdigest() == sha
    with pytest.raises(FileExistsError): current.build(output, receipt, tmp_path)


@pytest.mark.parametrize("mode", ["input", "location", "same_path", "receipt_exists", "symlink"])
def test_invalid_input_or_destination_fails_before_writing(tmp_path, mode):
    base = copied_inputs(tmp_path)
    output, receipt = base / "new.md", base / "receipt.json"
    if mode == "input": (base / current.PINS["original"][0]).write_text("changed")
    elif mode == "location": output = tmp_path / "outside.md"
    elif mode == "same_path": receipt = output
    elif mode == "receipt_exists": receipt.write_text("existing")
    elif mode == "symlink": output.symlink_to(base / current.PINS["original"][0])
    with pytest.raises(FileExistsError if mode in ("receipt_exists", "symlink") else ValueError):
        current.build(output, receipt, tmp_path)
    if mode != "symlink": assert not output.exists()
