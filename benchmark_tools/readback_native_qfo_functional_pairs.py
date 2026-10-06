"""Check native functional-pair overlap with an independent SQLite raw-table join."""

import argparse
import csv
from decimal import Decimal
import gzip
import json
from pathlib import Path
import sqlite3
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.export_native_factorial_progress import load, require
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


def raw_rows(path, metric):
    with gzip.open(path, "rt") as stream:
        if metric == "FAS":
            require(next(stream).rstrip("\n") == "Acc1\tAcc2\tFAS", "FAS raw header differs")
        else:
            headers = [next(stream).rstrip("\n") for _ in range(3)]
            require(headers[0].startswith(f"# {metric} Similarities between orthologs from ")
                    and headers[1].startswith("# Computing timestamp: ")
                    and headers[2] == f"# Protein ID 1<tab>Protein ID 2<tab>{metric} Similarity",
                    "GO/EC raw header differs")
        for row in csv.reader(stream, delimiter="\t"):
            require(len(row) == 3 and row[0] and row[1] and row[0] != row[1], "Malformed raw pair")
            a, b = sorted(row[:2])
            score = Decimal(row[2])
            require(score.is_finite() and 0 <= score <= 1, "Invalid raw score")
            if metric != "FAS":
                scaled = score * 1000000
                require(scaled == scaled.to_integral_value(), "GO/EC score is not an integer millionth")
            yield a, b, float(score) if metric == "FAS" else int(scaled)


def joined(left_path, right_path, metric):
    with sqlite3.connect(":memory:") as db:
        storage = "REAL" if metric == "FAS" else "INTEGER"
        for name, path in (("lhs", left_path), ("rhs", right_path)):
            db.execute(f"CREATE TABLE {name} (a TEXT, b TEXT, score {storage}, PRIMARY KEY (a,b)) WITHOUT ROWID")
            db.executemany(f"INSERT INTO {name} VALUES (?,?,?)", raw_rows(path, metric))
        left = db.execute("SELECT count(*), sum(score) FROM lhs").fetchone()
        right = db.execute("SELECT count(*), sum(score) FROM rhs").fetchone()
        shared = db.execute("SELECT count(*), coalesce(sum(lhs.score != rhs.score),0), "
            "max(abs(lhs.score-rhs.score)), coalesce(sum(lhs.score),0), coalesce(sum(rhs.score),0) "
            "FROM lhs JOIN rhs USING (a,b)").fetchone()
    return dict(left_pairs=left[0], right_pairs=right[0], left_sum=left[1], right_sum=right[1],
        shared_pairs=shared[0], differing_shared_pairs=shared[1], maximum_shared_difference=shared[2],
        shared_left_sum=shared[3], shared_right_sum=shared[4])


def readback(composition_path, composition_sha):
    evidence = []
    composition, ref = load(composition_path, composition_sha, evidence)
    require(composition["schema"] == "native_qfo_functional_pair_composition_v1"
            and composition["source"] == record(Path(__file__).with_name("compare_native_qfo_functional_pairs.py"))
            and composition["cells"] == ["p0_c0_r1", "p0_c0_r0"]
            and all(composition[k] is False for k in ("publication_ready", "uncertainty_admitted",
                "scientific_timings_admitted", "new_scoring_or_admission")), "Changed composition scope/source")
    endpoints = {(row["cell"], row["metric"]): row for row in composition["endpoints"]}
    require(len(endpoints) == len(composition["endpoints"]) == 6 and
            [row["metric"] for row in composition["comparisons"]] == ["GO", "EC", "FAS"],
            "Changed functional endpoint inventory")
    rows = []
    for comparison in composition["comparisons"]:
        metric = comparison["metric"]
        require(comparison["left"] == "p0_c0_r1" and comparison["right"] == "p0_c0_r0", "Changed direction")
        source_rows = [endpoints[cell, metric] for cell in composition["cells"]]
        for row in source_rows:
            check(row["raw"])
            evidence.append(row["raw"])
        actual = joined(*(Path(row["raw"]["path"]) for row in source_rows), metric)
        original = comparison["result"]
        if metric == "FAS":
            for name, value in (("left_sample_pairs", actual["left_pairs"]),
                ("right_sample_pairs", actual["right_pairs"]), ("shared_sample_pairs", actual["shared_pairs"]),
                ("shared_pairs_with_different_serialized_scores", actual["differing_shared_pairs"]),
                ("maximum_shared_absolute_difference", actual["maximum_shared_difference"])):
                require(original[name] == value, "FAS SQL join differs: " + name)
            require(abs(actual["left_sum"]/actual["left_pairs"] - actual["right_sum"]/actual["right_pairs"]
                        - original["original_sample_mean_difference"]) < 1e-12, "FAS SQL mean difference differs")
        else:
            checks = dict(left_pairs="left_pairs", right_pairs="right_pairs", shared_pairs="shared_pairs",
                left_score_sum_millionths="left_sum", right_score_sum_millionths="right_sum",
                shared_left_sum_millionths="shared_left_sum", shared_right_sum_millionths="shared_right_sum",
                shared_pairs_with_different_serialized_scores="differing_shared_pairs",
                maximum_shared_absolute_difference_millionths="maximum_shared_difference")
            require(all(original[name] == actual[key] for name, key in checks.items()), "GO/EC SQL join differs")
        for row, count, total in zip(source_rows, (actual["left_pairs"], actual["right_pairs"]),
                                     (actual["left_sum"], actual["right_sum"])):
            mean = total/count/(1 if metric == "FAS" else 1000000)
            require(row["scored_pairs"] == count and abs(mean-row["raw_mean"]) <= 1e-12,
                    "SQL endpoint count/mean differs")
        rows.append(dict(metric=metric, independent_sql_join=actual))
    for item in evidence:
        check(item)
    return dict(schema="native_qfo_functional_pair_sql_readback_v1", composition=ref,
        source=record(__file__), evidence=evidence, comparisons=rows, sqlite_version=sqlite3.sqlite_version,
        checked_raw_rows=sum(r["independent_sql_join"][name] for r in rows for name in ("left_pairs", "right_pairs")),
        independent_parsing_and_join=True, new_scoring_or_admission=False,
        uncertainty_admitted=False, publication_ready=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("composition", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--composition-sha256", required=True)
    args = parser.parse_args()
    require(not args.output.exists() and not args.output.is_symlink(), "Output already exists")
    result = readback(args.composition, args.composition_sha256)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps(dict(checked_raw_rows=result["checked_raw_rows"], comparisons=len(result["comparisons"]))))
