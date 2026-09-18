import pytest

from benchmark_tools.diagnose_vgnc_block_influence import incidents, metrics, deleted_scores


def test_cross_block_row_counted_once_per_deletion():
    rows = [("a","b","TP","A","A","s1","s2"), ("a","c","FP","A","C","s1","s2")]
    total, removed = incidents(rows, {"A":"A","C":"C","D":"D"})
    assert total == {"TP":1,"FP":1,"FN":0}
    assert removed["A"] == total
    assert removed["C"] == {"TP":0,"FP":1,"FN":0}
    assert removed["D"] == {"TP":0,"FP":0,"FN":0}
    assert deleted_scores(total, removed["C"])[1] == {"precision":1.,"recall":1.,"f1":1.}


def test_merged_labels_do_not_double_remove():
    _, removed = incidents([("a","b","TP","A","B","s1","s2")], {"A":"A","B":"A"})
    assert removed["A"]["TP"] == 1


def test_undefined_ratios_are_explicit():
    assert metrics(dict(TP=0,FP=0,FN=0)) == dict(precision=None,recall=None,f1=None)


@pytest.mark.parametrize("row", [("a","a","TP","A","A","s","s"), ("a","b","bad","A","A","s","s"), ("a","b","TP","X","A","s","s")])
def test_bad_rows(row):
    with pytest.raises(ValueError):
        incidents([row], {"A":"A"})


def test_conflicting_categories_rejected():
    with pytest.raises(ValueError):
        incidents([("a","b","TP","A","A","s1","s2"), ("b","a","FN","A","A","s2","s1")], {"A":"A"})


def test_negative_remaining_counts_rejected():
    with pytest.raises(ValueError):
        deleted_scores(dict(TP=1,FP=0,FN=0), dict(TP=2,FP=0,FN=0))


def test_retained_result_against_raw_rows_and_complete_table():
    import csv
    import gzip
    import json
    from collections import Counter
    from pathlib import Path

    from benchmark_tools.prepare_ob_candidate_neighborhood import check

    root = Path(__file__).resolve().parents[2]
    report_path = root / "benchmark_tools/results/vgnc_block_influence_20260918.json"
    if not report_path.exists():
        pytest.skip("Retained diagnostic not available")
    report = json.loads(report_path.read_text())
    if not Path(report["table"]["path"]).exists():
        pytest.skip("Local raw diagnostic table not available")
    check(report["table"])
    for item in report["inputs"]:
        check(item)
    assert report["uncertainty_admitted"] is False
    assert report["publication_ready"] is False
    merged = report["reference"]["merged_label_groups"]
    mapping = {label: min(group) for group in merged for label in group}
    with open(report["table"]["path"], newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    assert len(rows) == 4 * report["reference"]["reference_blocks"]
    scores = {}
    for stage in report["stages"]:
        index = stage["stage_index"]
        table = {r["block"]: r for r in rows if int(r["stage"]) == index}
        assert len(table) == report["reference"]["reference_blocks"]
        raw_path = next(item["path"] for item in report["inputs"]
                        if Path(item["path"]).name == f"VGNC_ohmm-checked-v2-{index}_raw.txt.gz")
        with gzip.open(raw_path, "rt", newline="") as stream:
            raw = list(csv.reader(stream, delimiter="\t"))
        # Direct row exclusion is independent of the production incident accumulator.
        for effect in stage["largest_absolute_f1_changes"]:
            block = effect["block"]
            kept = Counter(row[2] for row in raw
                           if mapping.get(row[3], row[3]) != block
                           and mapping.get(row[4], row[4]) != block)
            for category in ("TP", "FP", "FN"):
                assert kept[category] == int(table[block][f"remaining_{category.lower()}"])
                assert stage["full_counts"][category] - kept[category] == effect["removed"][category]
        scores[index] = {}
        for block, row in table.items():
            tp, fp, fn = (int(row[f"remaining_{c}"]) for c in ("tp", "fp", "fn"))
            assert float(row["precision"]) == tp / (tp + fp)
            assert float(row["recall"]) == tp / (tp + fn)
            score = tp / (tp + (fp + fn) / 2)
            assert float(row["f1"]) == score
            scores[index][block] = score
        changes = [s - stage["full_metrics"]["f1"] for s in scores[index].values()]
        assert min(changes) == stage["minimum_f1_change"]
        assert max(changes) == stage["maximum_f1_change"]
    for contrast in report["paired_contrasts"]:
        on, off = scores[contrast["on"]], scores[contrast["off"]]
        assert on.keys() == off.keys()
        values = [on[b] - off[b] for b in on]
        assert min(values) == contrast["minimum_deleted_f1_difference"]
        assert max(values) == contrast["maximum_deleted_f1_difference"]
        assert sum(v > 0 for v in values) == contrast["positive"]
        assert sum(v == 0 for v in values) == contrast["zero"]
        assert sum(v < 0 for v in values) == contrast["negative"]
