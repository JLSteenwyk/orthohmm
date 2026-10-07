"""Fixtures only before selected projection; actual artifacts are tested separately."""

import ast
from copy import deepcopy
from fractions import Fraction
import json
from pathlib import Path

import pytest

from benchmark_tools import export_native_qfo_three_cell_strata as export
from benchmark_tools import readback_native_qfo_three_cell_strata as failed_reader
from benchmark_tools import readback_native_qfo_three_cell_strata_v2 as readback

ROOT = Path(__file__).resolve().parents[2]


def rational(counts):
    p = Fraction(counts["TP"] + 2, counts["TP"] + counts["FP"] + 4)
    r = Fraction(counts["TP"] + 2, counts["TP"] + counts["FN"] + 4)
    return dict(F1=2 * p * r / (p + r), PPV=p, TPR=r)


def mean_points(values):
    if not values:
        return dict.fromkeys(export.METRICS)
    p, r = (sum(v[k] for v in values) / len(values) for k in ("PPV", "TPR"))
    return dict(F1=2 * p * r / (p + r), PPV=p, TPR=r)


def fixture():
    fs = [f"F{i:02}" for i in range(18)]
    members = {f: [f + f"_g{j:02}" for j in range(36 if i == 0 else 31)] for i, f in enumerate(fs)}
    family_rows, points = [], {}
    for c in export.CELLS[:2]:
        for i, f in enumerate(fs):
            counts = (dict(TP=10, FP=20, FN=30, TN=531 + 127 * (i == 0)) if c == export.CELLS[0]
                      else dict(TP=8, FP=10, FN=32, TN=541 + 127 * (i == 0)))
            value = {k: float(v) for k, v in rational(counts).items()}
            points[c, f] = value
            family_rows.append(dict(cell=c, family=f, counts_without_prior=counts, **value))
    primary = {"lower_entropy": fs[:9], "higher_entropy": fs[9:], "missing_entropy": []}
    secondary = {"concentrated": [], "explicit_fragment": [], "missing": [],
                 "no_explicit_fragment": fs, "not_concentrated": fs,
                 "not_short_relative": fs[7:], "short_relative": fs[:7]}
    bins = dict(sequence=dict(all=fs, **primary, **secondary),
                domain=dict(all=fs, low_types=fs[:12], high_types=fs[12:], low_repeat=fs[:15], high_repeat=fs[15:]),
                duplication=dict(all=fs, lower=fs[:9], upper=fs[9:], missing=[]))
    docs = dict(features=dict(family_memberships=members, primary_strata=primary,
                              secondary_strata=secondary, prediction_statistics_evaluated=False))
    for suite in export.SUITES:
        rows = [dict(cell=c, stratum=n, family_members=vs, families=len(vs),
                     prediction_semantics="resolved_native_pairs" if c == export.CELLS[1] else "group_clique",
                     **mean_points([points[c, f] for f in vs]))
                for c in export.CELLS[:2] for n, vs in bins[suite].items()]
        docs[suite] = dict(schema=f"native_qfo_swiss_{suite}_strata_v1", family_rows=deepcopy(family_rows),
                           memberships=members, bins=bins[suite], rows=rows,
                           **{k: False for k in readback.FLAGS[:4]})
        docs[suite + "_reader"] = dict(family_rows_checked=36, raw_rows_checked=21530)
    docs["domain"]["cells"] = [dict(cell=c, native_job_id=22435 + 2 * i, raw_file={}, admission={},
                                   **(dict(timing_admitted=False, timing_eligible=False) if i else {}))
                                for i, c in enumerate(export.CELLS[:2])]
    candidate_rows = []
    for i, f in enumerate(fs):
        counts = dict(TP=12, FP=30, FN=28, TN=521 + 127 * (i == 0))
        value = {k: float(v) for k, v in rational(counts).items()}
        points[export.CELLS[2], f] = value
        candidate_rows.append(dict(family=f, counts_without_prior=counts, represented_genes=members[f],
                                   statistics_with_prior=value))
    candidate = dict(cell=export.CELLS[2], index=8, native_job_id=22444, raw_file={}, admission={},
                     retained_family_records_identical=True, retained_aggregate_identical=True,
                     families=candidate_rows, aggregate=mean_points([points[export.CELLS[2], f] for f in fs]))
    docs["candidate"] = dict(schema="native_qfo_swiss_family_count_audit_v1",
        status="supplied_native_swiss_family_counts_verified", families=fs, reference_relation_count=10765,
        cells=[candidate], new_bootstrap_draws=0, historical_intervals_attached=False,
        **{k: False for k in readback.FLAGS[1:4]})
    docs["candidate_reader"] = dict(cells=[export.CELLS[0], export.CELLS[2]], families_checked=18,
        native_family_records_checked=36, candidate_pair_labels_matched=10765,
        rational_macro_points={c: mean_points([points[c, f] for f in fs]) for c in (export.CELLS[0], export.CELLS[2])})
    return docs


def projected(docs):
    return dict(schema="native_qfo_three_cell_strata_v1", limitations=export.LIMITATIONS,
                **export.SCOPE, **export.project(docs))


def test_complete_fixture_projection_and_independent_rational_reader():
    docs = fixture()
    report = projected(docs)
    result = readback.validate(docs, report)
    assert result == dict(families_checked=18, proteins_checked=563, family_rows_checked=54,
                          score_rows_checked=60, differences_checked=40, inherited_score_rows_reproduced=40)
    assert len(report["rows"]) == 60 and len(report["differences"]) == 40
    assert all(r["F1"] is None for r in report["rows"] if r["status"] == "empty_bin")
    assert all(r["TPR"] is None for r in report["differences"] if r["status"] == "empty_bin")
    assert report["cells"][1]["timing_admitted"] is False


@pytest.mark.parametrize("counts", [dict(TP=0, FP=0, FN=1, TN=0), dict(TP=3, FP=7, FN=5, TN=9),
                                    dict(TP=0, FP=100, FN=1000, TN=40)])
def test_qfo_halving_and_prior(counts):
    for k, v in rational(counts).items():
        assert export.statistics(counts)[k] == pytest.approx(float(v), abs=1e-12)
        assert readback.points(counts)[k] == v


def test_macro_is_not_pair_pooling_or_mean_family_f1():
    counts = [dict(TP=10, FP=100, FN=0, TN=1), dict(TP=1000, FP=0, FN=10000, TN=1)]
    values = [export.statistics(c) for c in counts]
    macro = export.aggregate(values)
    pooled = export.statistics({k: sum(c[k] for c in counts) for k in export.COUNTS})
    assert abs(macro["F1"] - pooled["F1"]) > .01
    assert abs(macro["F1"] - sum(v["F1"] for v in values) / 2) > .01
    assert macro["F1"] == pytest.approx(float(readback.macro([readback.points(c) for c in counts])["F1"]))


@pytest.mark.parametrize("bad", [dict(TP=True, FP=0, FN=0, TN=1), dict(TP=-1, FP=0, FN=0, TN=1),
                                 dict(TP=1.0, FP=0, FN=0, TN=1), dict(TP=0, FP=0, FN=0, TN=0),
                                 dict(TP=1, FP=2, FN=3)])
def test_invalid_counts_rejected_by_both_implementations(bad):
    for fn in (export.statistics, readback.points):
        with pytest.raises(ValueError):
            fn(bad)


@pytest.mark.parametrize("change", ["member", "bin", "old_row", "old_counts", "candidate_count",
                                    "candidate_prior", "reference", "source_scope", "timing", "duplicate_family"])
def test_input_tampering_rejected(change):
    docs = fixture()
    if change == "member":
        docs["candidate"]["cells"][0]["families"][0]["represented_genes"] = ["wrong"]
    elif change == "bin":
        docs["domain"]["bins"] = deepcopy(docs["domain"]["bins"])
        docs["domain"]["bins"]["high_types"] = []
    elif change == "old_row":
        docs["sequence"]["rows"][0]["F1"] += .01
    elif change == "old_counts":
        docs["domain"]["family_rows"][0]["counts_without_prior"]["TN"] += 1
    elif change == "candidate_count":
        docs["candidate"]["cells"][0]["families"][0]["counts_without_prior"]["TN"] += 1
    elif change == "candidate_prior":
        docs["candidate"]["cells"][0]["families"][0]["statistics_with_prior"]["F1"] += .01
    elif change == "reference":
        docs["candidate"]["reference_relation_count"] += 1
    elif change == "source_scope":
        docs["sequence"]["new_uncertainty"] = True
    elif change == "timing":
        docs["domain"]["cells"][1]["timing_admitted"] = True
    else:
        docs["candidate"]["cells"][0]["families"][0]["family"] = "F01"
    with pytest.raises(ValueError):
        export.project(docs)
    with pytest.raises(ValueError):
        readback.validate(docs, projected(fixture()))


@pytest.mark.parametrize("change", ["score", "empty", "delta", "missing", "reorder", "bin", "family", "member", "timing"])
def test_output_tampering_rejected(change):
    docs = fixture()
    report = projected(docs)
    if change == "score":
        report["rows"][0]["PPV"] += .01
    elif change == "empty":
        next(r for r in report["rows"] if r["status"] == "empty_bin")["F1"] = 0
    elif change == "delta":
        report["differences"][0]["F1"] += .01
    elif change == "missing":
        report["rows"].pop()
    elif change == "reorder":
        report["rows"].reverse()
    elif change == "bin":
        report["bins"] = deepcopy(report["bins"])
        report["bins"]["duplication"]["upper"] = []
    elif change == "family":
        report["family_rows"][0]["counts_without_prior"] = dict(TP=999, FP=0, FN=0, TN=0)
    elif change == "member":
        report["memberships"] = deepcopy(report["memberships"])
        report["memberships"]["F00"].pop()
    else:
        report["cells"][1]["timing_eligible"] = True
    with pytest.raises(ValueError):
        readback.validate(docs, report)


@pytest.mark.parametrize("flag", [*readback.FLAGS, "new_bootstrap_draws"])
def test_scope_inflation_rejected(flag):
    docs = fixture()
    report = projected(docs)
    report[flag] = 1 if flag == "new_bootstrap_draws" else True
    with pytest.raises(ValueError):
        readback.validate(docs, report)


def test_tables_roundtrip_and_tampering(tmp_path):
    report = projected(fixture())
    # JSON reload reproduces alphabetical bin-key order of the actual report.
    report = json.loads(json.dumps(report, sort_keys=True))
    score = tmp_path / "scores.tsv"
    export.write_tsv(score, report["rows"], export.SCORE_FIELDS)
    readback.tsv(str(score), report["rows"], readback.SCORE_FIELDS)
    table = tmp_path / "TABLE.md"
    table.write_text(export.table(report))
    readback.human_table(table, report)
    text = table.read_text()
    table.write_text(text.replace("+", "-"))
    with pytest.raises(ValueError):
        readback.human_table(table, report)
    score.write_text(score.read_text().replace("descriptive", "invented", 1))
    with pytest.raises(ValueError):
        readback.tsv(score, report["rows"], readback.SCORE_FIELDS)


@pytest.mark.parametrize("numeric", ["nan", "inf", "garbage"])
def test_invalid_tsv_metric(tmp_path, numeric):
    score = tmp_path / "scores.tsv"
    rows = projected(fixture())["rows"]
    export.write_tsv(score, rows, export.SCORE_FIELDS)
    score.write_text(score.read_text().replace(str(rows[0]["F1"]), numeric, 1))
    with pytest.raises(ValueError):
        readback.tsv(score, rows, readback.SCORE_FIELDS)


def test_refuse_existing_output_before_loading_inputs(tmp_path):
    with pytest.raises(ValueError, match="already exists"):
        export.export(tmp_path / "missing_repo", tmp_path)


def test_sources_are_stdlib_only_and_reader_does_not_import_exporter():
    allowed = {"argparse", "csv", "fractions", "hashlib", "json", "math", "pathlib"}
    for module in (export, readback):
        tree = ast.parse(Path(module.__file__).read_text())
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                assert all(n.name in allowed for n in node.names)
            elif isinstance(node, ast.ImportFrom):
                assert node.module in allowed
    assert readback.PINS == export.PINS


def test_frozen_direct_inputs_and_sources_match_without_selected_projection():
    docs, refs, checked = export.load_inputs(ROOT)
    assert len(refs) == 10 and len(checked) == 19
    assert docs["candidate_reader"]["audit"] == refs["candidate"]
    for suite in export.SUITES:
        assert docs[suite + "_reader"]["report"] == refs[suite]


def test_original_failed_reader_preserved_and_string_path_bug_is_explicit(tmp_path):
    assert export.record(failed_reader.__file__)["sha256"] == (
        "f59bd8cc924e5cbd0a6ef253a2ed47c243e92d1aaf2155dd3bc1ac286390583f"
    )
    score = tmp_path / "scores.tsv"
    rows = projected(fixture())["rows"]
    export.write_tsv(score, rows, export.SCORE_FIELDS)
    with pytest.raises(AttributeError, match="open"):
        failed_reader.tsv(str(score), rows, failed_reader.SCORE_FIELDS)
    readback.tsv(str(score), rows, readback.SCORE_FIELDS)


def test_full_fixture_export_and_json_string_path_readback(tmp_path, monkeypatch):
    docs = fixture()
    root = tmp_path / "benchmark_tools/results"
    root.mkdir(parents=True)
    refs, sources = {}, []
    for key in (*export.SUITES, *(s + "_reader" for s in export.SUITES),
                "features", "candidate", "candidate_reader"):
        path = tmp_path / (key + ".py")
        path.write_text("# Fixture bound source: " + key + "\n")
        docs[key]["source"] = export.record(path)
        sources.append(docs[key]["source"])
    names = {k: name for k, (name, _) in export.PINS.items()}

    def save(key):
        path = root / names[key]
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(docs[key], sort_keys=True))
        refs[key] = export.record(path)

    save("features")
    docs["sequence"]["strata"] = refs["features"]
    for suite in export.SUITES:
        save(suite)
        docs[suite + "_reader"]["report"] = refs[suite]
        save(suite + "_reader")
    save("candidate")
    docs["candidate_reader"]["audit"] = refs["candidate"]
    save("candidate_reader")
    path = root / names["protocol"]
    path.write_text("# Fixture protocol\n")
    refs["protocol"] = export.record(path)
    pins = {k: (names[k], refs[k]["sha256"]) for k in export.PINS}
    monkeypatch.setattr(export, "PINS", pins)
    monkeypatch.setattr(readback, "PINS", pins)
    out = tmp_path / "fresh_projection"
    result = export.export(tmp_path, out)
    assert result["source"] == export.record(export.__file__)
    report = export.record(out / "report.json")
    failure_path = tmp_path / "fixture_failure.json"
    failure_path.write_text(json.dumps(dict(report=report, execution_exit_code=1, output_written=False,
        report_or_export_overwritten=False, source=export.record(failed_reader.__file__))))
    amendment = tmp_path / "fixture_amendment.md"
    amendment.write_text("# Fixture repair amendment\n")
    monkeypatch.setattr(readback, "REPAIR_PINS", {
        "failure": (str(failure_path), export.record(failure_path)["sha256"]),
        "amendment": (str(amendment), export.record(amendment)["sha256"]),
    })
    checked = readback.verify(str(out / "report.json"), report["sha256"])
    assert checked["score_rows_checked"] == 60
    assert checked["differences_checked"] == 40
    with pytest.raises(ValueError, match="already exists"):
        export.export(tmp_path, out)
    (out / "scores.tsv").write_text((out / "scores.tsv").read_text() + "tampered\n")
    with pytest.raises(ValueError, match="output"):
        readback.verify(out / "report.json", report["sha256"])
