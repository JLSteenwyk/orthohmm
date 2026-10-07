"""Wholly invented domain coordinates and counts; no selected result evaluation."""

import copy
import inspect
import itertools
import json
from pathlib import Path
import random
import subprocess
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pytest

from benchmark_tools import prepare_swiss_ordered_pfam as primary
from benchmark_tools import readback_swiss_ordered_pfam as independent


def annotation(instances=(), length=100):
    rows = [dict(start=a, end=b, domain="pfam_" + d) for a, b, d in instances]
    types = sorted({r["domain"] for r in rows})
    return dict(length=length, ordered_pfam_instances=rows, pfam_types=types,
                pfam_type_count=len(types), pfam_instance_count=len(rows), has_repeated_pfam_type=len(rows) > len(types))


@pytest.mark.parametrize("instances,state,signature", [
    ([], "zero_pfam", None),
    ([(5, 8, "A")], "usable", ["pfam_A"]),
    ([(25, 30, "A"), (5, 8, "A"), (12, 20, "B")], "usable", ["pfam_A", "pfam_B", "pfam_A"]),
    ([(5, 8, "A"), (8, 11, "B")], "order_ambiguous", None),
    ([(5, 50, "A"), (9, 11, "B"), (20, 30, "C")], "order_ambiguous", None),
    ([(5, 8, "A"), (5, 7, "B")], "order_ambiguous", None),
    ([(5, 5, "A")], "order_ambiguous", None),
])
def test_primary_adjacent_and_independent_all_pairs_agree(instances, state, signature):
    row = annotation(instances)
    feature = primary.gene_feature(row, 100)
    assert independent.protein_descriptor(row, 100) == feature
    assert feature["status"] == state and feature["ordered_signature"] == signature
    assert len(feature["ordered_instances"]) == len(instances)


def test_invented_random_intervals_and_permutations_agree():
    rng = random.Random(937)
    for _ in range(200):
        tuples = set()
        for _ in range(rng.randrange(8)):
            a, b = sorted([rng.randrange(101), rng.randrange(101)])
            tuples.add((a, b, rng.choice("ABC")))
        rows = list(tuples)
        rng.shuffle(rows)
        row = annotation(rows)
        assert primary.gene_feature(row, 100) == independent.protein_descriptor(row, 100)


@pytest.mark.parametrize("corruption", ["length", "bool_length", "negative", "beyond", "reversed", "bool_coordinate",
                                        "duplicate", "types", "count", "bool_count", "repeat", "namespace"])
def test_invalid_annotation_is_not_silently_dropped(corruption):
    row = annotation([(5, 8, "A")])
    if corruption == "length":
        row["length"] = 99
    elif corruption == "bool_length":
        row["length"] = True
    elif corruption == "negative":
        row["ordered_pfam_instances"][0]["start"] = -1
    elif corruption == "beyond":
        row["ordered_pfam_instances"][0]["end"] = 101
    elif corruption == "reversed":
        row["ordered_pfam_instances"][0]["start"] = 9
    elif corruption == "bool_coordinate":
        row["ordered_pfam_instances"][0]["start"] = True
    elif corruption == "duplicate":
        row["ordered_pfam_instances"] *= 2
    elif corruption == "types":
        row["pfam_types"] = []
    elif corruption == "count":
        row["pfam_instance_count"] = 2
    elif corruption == "bool_count":
        row["pfam_type_count"] = True
    elif corruption == "repeat":
        row["has_repeated_pfam_type"] = True
    elif corruption == "namespace":
        row["ordered_pfam_instances"][0]["domain"] = "other_A"
    for implementation in (primary.gene_feature, independent.protein_descriptor):
        with pytest.raises(ValueError):
            implementation(row, 100)


def test_reordering_is_distinct_from_type_or_repeat_difference():
    inputs = dict(a=[(1, 3, "A"), (5, 8, "B")], b=[(1, 3, "B"), (5, 8, "A")],
                  c=[(1, 3, "A"), (5, 8, "B"), (9, 12, "B")], d=[(1, 3, "A"), (5, 8, "B")])
    genes = {g: primary.gene_feature(annotation(v), 100) for g, v in inputs.items()}
    members = sorted(genes)
    feature = primary.family_feature(members, genes)
    assert feature == independent.family_descriptor(members, genes)
    assert feature["same_multiset_comparable_pairs"] == 3
    assert feature["same_multiset_order_discordant_pairs"] == 2
    assert feature["distinct_usable_signatures"] == 3 and feature["stratum"] == primary.BINS[2]
    genes["c"] = primary.gene_feature(annotation([]), 100)
    assert primary.family_feature(members, genes)["stratum"] == primary.BINS[3]
    assert independent.family_descriptor(["c"], genes)["same_multiset_comparable_pairs"] == 0


@pytest.fixture
def toy_counts():
    memberships = dict(A=["a"], B=["b"], C=["c"])
    bins = dict(zip(primary.BINS, [["A", "B", "C"], ["A", "B"], ["C"], []]))
    rows = []
    for i, cell in enumerate(primary.CELLS):
        for j, family in enumerate(memberships):
            c = dict(TP=10 + 100 * j + i, FP=20 + j - i, FN=3 + 7 * j + i, TN=200 - j)
            rows.append(dict(cell=cell, family=family, counts_without_prior=c,
                             **{k: float(v) for k, v in independent.point(c).items()}))
    counts = dict(memberships=memberships, family_rows=rows,
                  cells=[dict(cell=c, **(dict(timing_eligible=False, timing_admitted=False)
                         if c == primary.CELLS[1] else {})) for c in primary.CELLS])
    features = dict(memberships=memberships, bins=bins)
    projected, diffs = primary.projected_rows(counts, bins)
    report = dict(schema="native_swiss_ordered_pfam_projection_v1", status="projection_constructed_unverified",
                  rows=projected, differences=diffs, bins=bins, **counts, **primary.SCOPE)
    return counts, features, report


def test_exact_rational_scores_all_cells_pp_units_and_empty_bins(toy_counts, tmp_path):
    counts, features, report = toy_counts
    rows, diffs = independent.score_readback(report, counts, features)
    assert len(rows) == 12 and len(diffs) == 8
    assert all(r["F1"] is None for r in rows[-3:])
    assert all(d["F1_pp"] is None for d in diffs[-2:])
    outputs = primary.write_tables(tmp_path, report["rows"], report["differences"])
    independent.tables_readback(outputs, rows, diffs)
    assert "| some_members_unusable | 0 | p0_c0_r0 | NA | NA | NA |" in Path(outputs["table"]["path"]).read_text()
    mean_f1 = sum(r["F1"] for r in counts["family_rows"][:3]) / 3
    pooled = independent.point({k: sum(r["counts_without_prior"][k] for r in counts["family_rows"][:3])
                                for k in ("TP", "FP", "FN", "TN")})
    assert abs(float(rows[0]["F1"]) - mean_f1) > .001
    assert abs(float(rows[0]["F1"]) - float(pooled["F1"])) > .001


@pytest.mark.parametrize("corruption", ["value", "members", "missing", "semantics", "unit", "failed_timing", "count", "bin"])
def test_independent_scores_reject_corruption(toy_counts, corruption):
    counts, features, report = copy.deepcopy(toy_counts)
    if corruption == "value":
        report["rows"][0]["F1"] += .01
    elif corruption == "members":
        report["rows"][0]["family_members"] = []
    elif corruption == "missing":
        report["differences"].pop()
    elif corruption == "semantics":
        report["rows"][1]["prediction_semantics"] = "group_clique"
    elif corruption == "unit":
        report["differences"][0]["F1_pp"] /= 100
    elif corruption == "failed_timing":
        report["cells"] = copy.deepcopy(report["cells"])
        report["cells"][1]["timing_eligible"] = True
    elif corruption == "count":
        report["family_rows"] = copy.deepcopy(report["family_rows"])
        report["family_rows"][0]["counts_without_prior"]["TN"] += 1
    elif corruption == "bin":
        report["bins"] = copy.deepcopy(report["bins"])
        report["bins"][primary.BINS[1]] = []
    with pytest.raises(ValueError):
        independent.score_readback(report, counts, features)


@pytest.mark.parametrize("corruption", ["bool", "negative", "zero", "missing", "duplicate_row", "partition"])
def test_primary_rejects_invalid_counts_or_partition(toy_counts, corruption):
    counts, features, _ = toy_counts
    if corruption in {"bool", "negative", "zero", "missing"}:
        c = counts["family_rows"][0]["counts_without_prior"]
        if corruption == "bool":
            c["TP"] = True
        elif corruption == "negative":
            c["FP"] = -1
        elif corruption == "zero":
            c.update(dict.fromkeys(c, 0))
        else:
            c.pop("TN")
    elif corruption == "duplicate_row":
        counts["family_rows"][1] = counts["family_rows"][0]
    else:
        features["bins"][primary.BINS[2]] = ["A", "C"]
    with pytest.raises(ValueError):
        primary.projected_rows(counts, features["bins"])


@pytest.mark.parametrize("corruption", ["number", "header", "extra_row", "human", "scope"])
def test_tables_reject_corruption(toy_counts, tmp_path, corruption):
    counts, features, report = toy_counts
    rows, diffs = independent.score_readback(report, counts, features)
    outputs = primary.write_tables(tmp_path, report["rows"], report["differences"])
    key = "table" if corruption in {"human", "scope"} else "scores"
    path = Path(outputs[key]["path"])
    text = path.read_text()
    if corruption == "number":
        text = text.replace(str(report["rows"][0]["F1"]), "0.01", 1)
    elif corruption == "header":
        text = text.replace("F1", "bad", 1)
    elif corruption == "extra_row":
        text += text.splitlines()[1] + "\n"
    elif corruption == "human":
        text = text.replace(f"{100 * report['rows'][0]['F1']:.3f}", "0.001", 1)
    else:
        text = text.replace("Initial HMM search on, downstream profiles off", "wrong")
    path.write_text(text)
    with pytest.raises(ValueError):
        independent.tables_readback(outputs, rows, diffs)


def full_invented_fixture(tmp_path, monkeypatch):
    memberships = {f"F{i:02d}": [f"F{i:02d}_g{j:03d}" for j in range(31 if i < 17 else 36)] for i in range(18)}
    genes, runs = {}, []
    for i, (family, members) in enumerate(memberships.items()):
        path = tmp_path / (family + ".faa")
        records = [SeqRecord(Seq("ACDEFGHIKL--"), id=g, description="") for g in members]
        SeqIO.write(records, path, "fasta")
        runs.append(dict(family=family, columns=12, alignment=primary.record(path)))
        for j, g in enumerate(members):
            pattern = [] if i == 2 else [(0, 2, "A"), (4, 6, "B")] if i != 1 or j % 2 == 0 else [(0, 2, "B"), (4, 6, "A")]
            genes[g] = annotation(pattern, length=10)
    docs = dict(annotations=dict(genes=genes, prediction_statistics_evaluated=False),
                alignments=dict(memberships=memberships, runs=runs, failed_families=[], prediction_statistics_evaluated=False),
                domain_report=dict(should_not_parse=True), protocol="invented prospective protocol")
    refs = {}
    for k in ("annotations", "alignments", "domain_report", "protocol"):
        path = tmp_path / (k + (".txt" if k == "protocol" else ".json"))
        path.write_text(docs[k] if k == "protocol" else json.dumps(docs[k]))
        refs[k] = primary.record(path)
    docs["domain_reader"] = dict(proteins_checked=563, report=refs["domain_report"], checked_inputs=[refs["annotations"]])
    path = tmp_path / "domain_reader.json"
    path.write_text(json.dumps(docs["domain_reader"]))
    refs["domain_reader"] = primary.record(path)
    source_path = tmp_path / "invented_source.py"
    source_path.write_text("# invented fixture\n")
    source_ref = primary.record(source_path)
    pins = {k: (Path(refs[k]["path"]).name, refs[k]["sha256"]) for k in primary.FEATURE_PINS}
    monkeypatch.setattr(primary, "FEATURE_PINS", pins)
    monkeypatch.setattr(primary, "source", lambda *args: source_ref)
    monkeypatch.setattr(independent, "PINS", {k: r["sha256"] for k, r in refs.items()})
    monkeypatch.setattr(independent, "committed_source", lambda *args: None)
    return docs


def test_complete_feature_stage_and_independent_readback_after_sorted_json(tmp_path, monkeypatch):
    full_invented_fixture(tmp_path, monkeypatch)
    original = primary.pinned_inputs
    def stage_checked(repo, pins, parse):
        assert set(pins) == {"annotations", "domain_report", "domain_reader", "alignments", "protocol"}
        assert set(parse) == {"annotations", "domain_reader", "alignments"}
        return original(repo, pins, parse)
    monkeypatch.setattr(primary, "pinned_inputs", stage_checked)
    path = tmp_path / "features.json"
    report = primary.construct(tmp_path, path, "invented")
    assert len(report["genes"]) == 563 and len(report["families"]) == 18
    assert [len(report["bins"][k]) for k in primary.BINS] == [18, 16, 1, 1]
    result = independent.verify("features", path, tmp_path, "invented")
    assert result["status"] == "ordered_features_verified" and result["genes_checked"] == 563
    assert result["protein_states"] == dict(usable=532, order_ambiguous=0, zero_pfam=31)


def test_full_invented_four_stage_provenance_and_counts_roundtrip(tmp_path, monkeypatch):
    full_invented_fixture(tmp_path, monkeypatch)
    features_path = tmp_path / "features.json"
    features = primary.construct(tmp_path, features_path, "invented")
    feature_reader_path = tmp_path / "feature_reader.json"
    primary.save(feature_reader_path, independent.verify("features", features_path, tmp_path, "invented"))
    rows = []
    for i, cell in enumerate(primary.CELLS):
        for j, family in enumerate(features["memberships"]):
            c = dict(TP=10 + 30 * j + i, FP=20 + j - i, FN=3 + 7 * j + i, TN=200 - j)
            rows.append(dict(cell=cell, family=family, counts_without_prior=c,
                             **{k: float(v) for k, v in independent.point(c).items()}))
    counts = dict(memberships=features["memberships"], family_rows=rows,
                  cells=[dict(cell=c, **(dict(timing_eligible=False, timing_admitted=False)
                         if c == primary.CELLS[1] else {})) for c in primary.CELLS])
    counts_path = tmp_path / "counts.json"
    primary.save(counts_path, counts)
    count_reader_path = tmp_path / "counts_reader.json"
    primary.save(count_reader_path, dict(report=primary.record(counts_path), family_rows_checked=54))
    refs = dict(counts=primary.record(counts_path), counts_reader=primary.record(count_reader_path))
    monkeypatch.setattr(primary, "COUNT_PINS", {k: (Path(v["path"]).name, v["sha256"]) for k, v in refs.items()})
    independent.PINS.update({k: v["sha256"] for k, v in refs.items()})
    output = tmp_path / "projected"
    report = primary.project(tmp_path, features_path, feature_reader_path, output, "invented")
    result = independent.verify("scores", output / "report.json", tmp_path, "invented")
    assert result["status"] == "ordered_projection_verified"
    assert result["family_rows_checked"] == 54 and result["score_rows_checked"] == 12
    assert result["differences_checked"] == 8 and result["human_numeric_rows_checked"] == 20
    assert report["family_rows"] == counts["family_rows"]
    assert report["scientific_timings_admitted"] is False


@pytest.mark.parametrize("corruption", ["length", "membership", "coordinates", "descriptor", "bin", "scope", "binding", "checked_inputs"])
def test_full_feature_failure_never_produces_partial_success(tmp_path, monkeypatch, corruption):
    docs = full_invented_fixture(tmp_path, monkeypatch)
    if corruption in {"length", "membership", "coordinates"}:
        raw_path = tmp_path / ("alignments.json" if corruption == "membership" else "annotations.json")
        changed = docs["alignments"] if corruption == "membership" else docs["annotations"]
        if corruption == "membership":
            changed["memberships"]["F00"].pop()
        elif corruption == "length":
            changed["genes"]["F00_g000"]["length"] = 9
        else:
            changed["genes"]["F00_g000"]["ordered_pfam_instances"][0]["end"] = 11
        raw_path.write_text(json.dumps(changed))
        key = "alignments" if corruption == "membership" else "annotations"
        primary.FEATURE_PINS[key] = (raw_path.name, primary.record(raw_path)["sha256"])
        if key == "annotations":
            receipt_path = tmp_path / "domain_reader.json"
            receipt = json.loads(receipt_path.read_text())
            receipt["checked_inputs"] = [primary.record(raw_path)]
            receipt_path.write_text(json.dumps(receipt))
            primary.FEATURE_PINS["domain_reader"] = (receipt_path.name, primary.record(receipt_path)["sha256"])
        target = tmp_path / "features.json"
        with pytest.raises(ValueError):
            primary.construct(tmp_path, target, "invented")
        assert not target.exists()
    else:
        path = tmp_path / "features.json"
        report = primary.construct(tmp_path, path, "invented")
        if corruption == "descriptor":
            report["genes"]["F00_g000"]["ordered_signature"] = ["pfam_B", "pfam_A"]
        elif corruption == "bin":
            report["bins"][primary.BINS[1]].pop()
        elif corruption == "scope":
            report["publication_ready"] = True
        elif corruption == "binding":
            report["inputs"]["annotations"]["sha256"] = "0" * 64
        else:
            report["checked_inputs"].pop()
        path.write_text(json.dumps(report, sort_keys=True))
        with pytest.raises(ValueError):
            independent.verify("features", path, tmp_path, "invented")


@pytest.mark.parametrize("module", [primary, independent])
@pytest.mark.parametrize("kind", ["file", "dangling_symlink"])
def test_cli_refuses_occupied_output_before_any_input_access(tmp_path, module, kind):
    target = tmp_path / "occupied"
    if kind == "file":
        target.write_text("untouched")
    else:
        target.symlink_to(tmp_path / "missing")
    args = [sys.executable, module.__file__, "features", "--repo", str(tmp_path), "--source-commit", "unused", "--output", str(target)]
    if module is independent:
        args.extend(["--report", str(tmp_path / "absent")])
    done = subprocess.run(args, capture_output=True, text=True)
    assert done.returncode != 0 and "FileExistsError" in done.stderr
    assert "untouched" == target.read_text() if kind == "file" else target.is_symlink()


def test_features_are_separate_and_reader_does_not_call_primary_implementation():
    text = inspect.getsource(independent)
    assert "import prepare_swiss_ordered_pfam" not in text and "from benchmark_tools" not in text
    assert "itertools.combinations(raw, 2)" in text
    assert "COUNT_PINS" not in inspect.getsource(primary.construct)
    assert all(len(sha) == 64 for _, sha in itertools.chain(primary.FEATURE_PINS.values(), primary.COUNT_PINS.values()))
    for k, (_, sha) in itertools.chain(primary.FEATURE_PINS.items(), primary.COUNT_PINS.items()):
        assert independent.PINS[k] == sha
