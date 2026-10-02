from copy import deepcopy
from fractions import Fraction
import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from benchmark_tools import check_swiss_descriptive_tables as checker

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "benchmark_tools/results"
GUARDED_REPLAY = """
import json, os, runpy, sys
from pathlib import Path
blocked, script, results, output = [Path(p).resolve() for p in sys.argv[1:]]
events = []
def audit(event, args):
    if event == 'subprocess.Popen':
        raise PermissionError('subprocess forbidden in derived-table replay')
    if event in ('open', 'os.chdir') and args and isinstance(args[0], (str, bytes, os.PathLike)):
        path = Path(os.fsdecode(args[0])).resolve()
        if path == blocked or path.is_relative_to(blocked):
            events.append(event)
            raise PermissionError('original checkout access rejected')
sys.addaudithook(audit)
try:
    (blocked / 'benchmark_tools/results/qfo_recovered_swiss_uncertainty_22178.json').read_bytes()
except PermissionError:
    assert events == ['open']
else:
    raise RuntimeError('canary was not rejected')
events.clear()
sys.argv = [str(script), '--results', str(results), '--output', str(output)]
runpy.run_path(str(script), run_name='__main__')
report = json.loads(output.read_text())
print(json.dumps(dict(canary_rejected=True, original_path_events_after_canary=len(events),
    rows=report['rows'], cells=report['score_and_difference_cells'],
    checked_records_within_copy=all(Path(r['path']).is_relative_to(results) for r in report['checked_records']),
    checker_within_copy=Path(report['checker']['path']) == script)))
"""


def payload_paths():
    paths = [Path(checker.COUNTS)]
    for kind, (feature, _) in checker.PANELS.items():
        paths.append(Path(feature))
        directory = Path(f"swiss_{kind}_strata_20260926")
        paths.extend(directory / name for name in ("manifest.json", "scores.tsv", "scores.md"))
    return paths


@pytest.fixture
def copied(tmp_path):
    root = tmp_path / "table review with spaces"
    for name in payload_paths():
        target = root / "results" / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(RESULTS / name, target)
    shutil.copyfile(ROOT / "benchmark_tools/check_swiss_descriptive_tables.py", root / "checker.py")
    return root


def fixture_counts():
    zero = dict(TP=0, FP=0, FN=0, TN=0)
    other = dict(TP=8, FP=0, FN=12, TN=0)
    methods = []
    for name, raw in ((checker.REFERENCE, [zero, zero]), ("candidate", [zero, other])):
        families = []
        for family, row in zip(("A", "B"), raw):
            p = Fraction(row["TP"]+2, row["TP"]+row["FP"]+4)
            r = Fraction(row["TP"]+2, row["TP"]+row["FN"]+4)
            families.append(dict(family=family, counts_without_prior=dict(row),
                statistics_with_prior=dict(F1=float(2*p*r/(p+r)), PPV=float(p), TPR=float(r))))
        methods.append(dict(method=name, status="counts_verified", families=families,
                            prediction_semantics="fixture native pairs"))
    return dict(status="corrected_swiss_comparison_intervals_audited", scientific_inputs_admitted=True,
        families=["A", "B"], reconstructed_counts=dict(methods=methods),
        point_estimates={checker.REFERENCE: dict(F1=.5, PPV=.5, TPR=.5),
                        "candidate": dict(F1=float(Fraction(44, 81)), PPV=float(Fraction(2, 3)),
                                          TPR=float(Fraction(11, 24)))})


def test_exact_macro_harmonic_statistic_and_missingness():
    rows = checker.rational_rows(fixture_counts(), {"all": ["A", "B"], "empty": []}, True)
    candidate = next(r for r in rows if r["method"] == "candidate" and r["stratum"] == "all")
    assert candidate["F1"] == Fraction(44, 81)
    assert candidate["F1"] != Fraction(19, 36)  # Mean family F1 is a different endpoint.
    assert candidate["delta_F1"] == Fraction(7, 162)
    for row in rows:
        if row["stratum"] == "empty":
            assert row["status"] == "empty_bin"
            assert all(row[k] is None for k in (*checker.METRICS, "delta_F1", "delta_PPV", "delta_TPR"))


def test_not_admitted_method_remains_missing():
    counts = fixture_counts()
    counts["reconstructed_counts"]["methods"][1].update(status="not_admitted", families=[])
    counts["point_estimates"]["candidate"] = None
    row = checker.rational_rows(counts, {"all": ["A", "B"]}, True)[1]
    assert row["status"] == "method_not_admitted" and row["F1"] is row["delta_F1"] is None


@pytest.mark.parametrize("fault", ["status", "admission", "family_duplicate", "family_order", "method_duplicate",
                                  "point_inventory", "raw_negative", "raw_bool", "raw_float", "raw_fields",
                                  "family_score", "point_score", "missing_imputation", "unknown_status"])
def test_malformed_counts_rejected(fault):
    counts = fixture_counts()
    method = counts["reconstructed_counts"]["methods"][1]
    family = method["families"][0]
    if fault == "status":
        counts["status"] = "unadmitted"
    elif fault == "admission":
        counts["scientific_inputs_admitted"] = False
    elif fault == "family_duplicate":
        counts["families"] = ["A", "A"]
    elif fault == "family_order":
        method["families"].reverse()
    elif fault == "method_duplicate":
        counts["reconstructed_counts"]["methods"].append(deepcopy(method))
    elif fault == "point_inventory":
        del counts["point_estimates"]["candidate"]
    elif fault.startswith("raw_"):
        if fault == "raw_fields":
            del family["counts_without_prior"]["TN"]
        else:
            family["counts_without_prior"]["TP"] = {"raw_negative": -1, "raw_bool": False, "raw_float": 0.}[fault]
    elif fault == "family_score":
        family["statistics_with_prior"]["F1"] = .6
    elif fault == "point_score":
        counts["point_estimates"]["candidate"]["F1"] = .6
    elif fault == "missing_imputation":
        method.update(status="not_admitted", families=[])
    else:
        method["status"] = "pending"
    with pytest.raises(ValueError):
        checker.rational_rows(counts, {"all": ["A", "B"]}, True)


@pytest.mark.parametrize("members", [["A", "A"], ["X"], ["A", "B", "X"]])
def test_bad_stratum_membership_rejected(members):
    with pytest.raises(ValueError, match="stratum membership"):
        checker.rational_rows(fixture_counts(), {"all": ["A", "B"], "bad": members}, False)


@pytest.mark.parametrize("value", [True, None, float("nan"), float("inf"), "NA", "nan", "inf", .7])
def test_nonfinite_or_wrong_numeric_values_rejected(value):
    with pytest.raises(ValueError):
        checker.number_matches(value, Fraction(1, 2))


def test_absolute_tolerance_is_fixed_and_missing_is_not_zero():
    checker.number_matches("0.500000000001", Fraction(1, 2))
    with pytest.raises(ValueError, match="Rational score differs"):
        checker.number_matches("0.5000000000011", Fraction(1, 2))
    with pytest.raises(ValueError, match="Missing score was imputed"):
        checker.number_matches(0., None)


@pytest.mark.parametrize("fault", ["duplicate", "missing", "extra_field", "semantics", "members", "families",
                                  "status", "score", "difference", "missing_zero", "family_bool", "family_float"])
def test_table_inventory_metadata_and_values_rejected(fault):
    expected = checker.rational_rows(fixture_counts(), {"all": ["A", "B"], "empty": []}, True)
    retained = [{k: float(v) if isinstance(v, Fraction) else deepcopy(v) for k, v in r.items()} for r in expected]
    if fault == "duplicate":
        retained.append(deepcopy(retained[0]))
    elif fault == "missing":
        retained.pop()
    elif fault == "missing_zero":
        retained[1]["F1"] = 0.
    elif fault in ("family_bool", "family_float"):
        retained[1]["families"] = False if fault == "family_bool" else 0.
    else:
        key, value = {"extra_field": ("extra", True), "semantics": ("prediction_semantics", "group cliques"),
            "members": ("family_members", ["A"]), "families": ("families", 1), "status": ("status", "pending"),
            "score": ("F1", .7), "difference": ("delta_F1", .2)}[fault]
        retained[0][key] = value
    with pytest.raises(ValueError):
        checker.compare_rows(retained, expected)


def test_retained_all_four_tables_verify(copied):
    result = checker.verify(copied / "results")
    assert result["rows"] == 208 and result["score_and_difference_cells"] == 984
    assert [p["rows"] for p in result["panels"]] == [88, 32, 56, 32]
    assert result["raw_inputs_revalidated"] is result["bootstrap_intervals_recomputed"] is False
    assert result["publication_ready"] is False
    assert len({r["path"] for r in result["checked_records"]}) == 17
    assert all(Path(r["path"]).is_relative_to(copied) for r in result["checked_records"])


@pytest.mark.parametrize("target", [checker.COUNTS, checker.PANELS["duplication"][0],
    "swiss_identity_strata_20260926/manifest.json", "swiss_fragment_strata_20260926/scores.tsv",
    "swiss_descriptive_strata_20260926/scores.md"])
def test_changed_component_bytes_rejected(copied, target):
    path = copied / "results" / target
    path.write_bytes(path.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="Changed retained bytes"):
        checker.verify(copied / "results")


def test_standalone_copy_and_no_overwrite(copied):
    output = copied / "report.json"
    argv = [sys.executable, "-I", "-B", str(copied / "checker.py"), "--results", str(copied / "results"),
            "--output", str(output)]
    subprocess.run(argv, cwd=copied, capture_output=True, text=True, check=True)
    original = output.read_bytes()
    assert json.loads(original)["rows"] == 208
    failed = subprocess.run(argv, cwd=copied, capture_output=True, text=True)
    assert failed.returncode != 0 and "FileExistsError" in failed.stderr
    assert output.read_bytes() == original


def test_copied_child_rejects_original_checkout_reads(copied):
    run = subprocess.run([sys.executable, "-I", "-B", "-c", GUARDED_REPLAY, str(ROOT),
        str(copied / "checker.py"), str(copied / "results"), str(copied / "guarded.json")],
        cwd=copied, capture_output=True, text=True, check=True)
    assert json.loads(run.stdout) == dict(canary_rejected=True, original_path_events_after_canary=0,
        rows=208, cells=984, checked_records_within_copy=True, checker_within_copy=True)


def test_payload_change_during_final_panel_is_detected(copied, monkeypatch):
    original = checker.compare_rows

    def change_after_last_comparison(retained, expected, *, tsv=False):
        original(retained, expected, tsv=tsv)
        if tsv and any(r["stratum"] == "lower_duplication_fraction" for r in expected):
            path = copied / "results" / checker.COUNTS
            path.write_bytes(path.read_bytes() + b"\n")

    monkeypatch.setattr(checker, "compare_rows", change_after_last_comparison)
    with pytest.raises(ValueError, match="Changed retained bytes"):
        checker.verify(copied / "results")
