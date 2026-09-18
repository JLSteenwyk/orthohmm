"""Audit rounded GO/EC pair scores and Darwin's native interval semantics."""

import argparse
from collections import Counter
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_qfo_swiss_comparators import COMPARISON_SHA, METHODS
from benchmark_tools.audit_qfo_swiss_counts import IMAGE_SHA
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check


def read_raw(path, metric):
    n, mean, m2 = 0, 0.0, 0.0
    seen, degrees = set(), Counter()
    with gzip.open(path, "rt") as stream:
        headers = [next(stream, "").rstrip("\n") for _ in range(3)]
        if not headers[0].startswith(f"# {metric} Similarities between orthologs from ") or not headers[1].startswith("# Computing timestamp: ") or headers[2] != f"# Protein ID 1<tab>Protein ID 2<tab>{metric} Similarity":
            raise ValueError("Unexpected raw header")
        for line in stream:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 3:
                raise ValueError("Malformed raw row")
            a, b, value = fields
            if not a or not b or a == b or not re.fullmatch(r"[01]\.\d{6}", value):
                raise ValueError("Invalid pair or score serialization")
            score = float(value)
            if not 0 <= score <= 1:
                raise ValueError("Invalid score")
            pair = tuple(sorted((a, b)))
            if pair in seen:
                raise ValueError("Duplicate undirected pair")
            seen.add(pair)
            degrees.update(pair)
            n += 1
            delta = score - mean
            mean += delta / n
            m2 += delta * (score - mean)
    if n < 2:
        raise ValueError("Insufficient raw observations")
    return {"pairs": n, "mean_from_rounded_raw": mean, "pair_iid_sem_from_rounded_raw": math.sqrt(m2 / (n - 1) / n),
            "proteins": len(degrees), "proteins_in_multiple_pairs": sum(v > 1 for v in degrees.values()),
            "maximum_protein_degree": max(degrees.values())}


def validate(summary, participant, critical):
    n = summary["pairs"]
    if n != participant["metric_x"]:
        raise ValueError("Raw count differs from assessed count")
    if not math.isfinite(critical) or not 1 < critical < 20:
        raise ValueError("Invalid native critical value")
    # Each raw value is rounded to six decimals. The sample-SD perturbation
    # bound is epsilon*sqrt(n/(n-1)); divide by sqrt(n) for the SEM bound.
    mean_tolerance = 5e-7 + 5e-9 + 1e-11
    width_tolerance = critical * 5e-7 / math.sqrt(n - 1) + 5e-11
    width = critical * summary["pair_iid_sem_from_rounded_raw"]
    if not math.isclose(summary["mean_from_rounded_raw"], participant["metric_y"], rel_tol=0, abs_tol=mean_tolerance):
        raise ValueError("Raw mean differs beyond serialization bound")
    if not math.isclose(width, participant["stderr_y"], rel_tol=0, abs_tol=width_tolerance):
        raise ValueError("Native interval differs beyond serialization bound")
    return {"native_t_critical": critical, "native_95_half_width_from_rounded_raw": width,
            "mean_absolute_tolerance": mean_tolerance, "half_width_absolute_tolerance": width_tolerance,
            "reported_mean": participant["metric_y"], "reported_stderr_field": participant["stderr_y"]}


def audit(repo):
    frozen = repo / "benchmark_tools/results/publication_comparison_orthomcl_complete_20260916.json"
    identity = record(frozen)
    if identity["sha256"] != COMPARISON_SHA:
        raise ValueError("Changed comparison")
    comparison = json.loads(frozen.read_text())
    if tuple(r["key"] for r in comparison["methods"]) != METHODS:
        raise ValueError("Changed method inventory")
    image = record(repo / "qfo_benchmark/scoring/container_cache/qfobenchmark-darwin-2022.1.img")
    if image["sha256"] != IMAGE_SHA:
        raise ValueError("Changed Darwin image")
    native_sources = {}
    for filename in ("/darwin/lib/Stat", "/benchmark/GoTest.drw", "/benchmark/EcTest.drw"):
        contents = subprocess.run(["singularity", "exec", image["path"], "cat", filename], capture_output=True, check=True, timeout=30).stdout
        native_sources[filename] = {"sha256": hashlib.sha256(contents).hexdigest(), "bytes": len(contents)}
    rows = []
    for metric in ("GO", "EC"):
        sources = [(r["key"], r["qfo"]["metric_details"][metric]["source"]) for r in comparison["methods"]]
        sources.extend((f"checked_v2_{i}", record(repo / f"qfo_benchmark/scoring/checked_v2_{i}/results/{metric}/{metric}.json")) for i in range(4))
        for method, source in sources:
            check(source)
            data = json.loads(Path(source["path"]).read_text())["datalink"]["inline_data"]
            if data["visualization"]["x_axis"] != "NR_ORTHOLOGS" or data["visualization"]["y_axis"] != "avg Schlicker":
                raise ValueError("Changed endpoint")
            paths = list(Path(source["path"]).parent.glob("*raw.txt.gz"))
            if len(paths) != 1 or len(data["challenge_participants"]) != 1:
                raise ValueError("Ambiguous raw evidence")
            raw = record(paths[0])
            summary = read_raw(paths[0], metric)
            check(raw)
            rows.append({"method": method, "metric": metric, "raw": raw, "aggregate": source,
                         "participant": data["challenge_participants"][0], **summary})
    program = "s := Stat('audit'):\n" + "\n".join(f"printf('CRITICAL\\t{i}\\t%.15g\\n', InverseStudent_t(0.95,{r['pairs'] - 1})):" for i, r in enumerate(rows)) + "\ndone;\n"
    native = subprocess.run(["singularity", "exec", image["path"], "darwin", "-E"], input=program, text=True, capture_output=True, check=True, timeout=60)
    values = re.findall(r"^CRITICAL\t(\d+)\t([0-9.eE+-]+)$", native.stdout, re.MULTILINE)
    if len(values) != len(rows) or [int(i) for i, _ in values] != list(range(len(rows))):
        raise ValueError("Incomplete native critical values")
    for row, (_, value) in zip(rows, values):
        row.update(validate(row, row["participant"], float(value)))
    check(image)
    return {"status": "all_retained_go_ec_arithmetic_verified_with_serialization_bounds", "source": record(__file__),
            "comparison": identity, "container": image, "native_sources": native_sources,
            "critical_value_program": program, "native_stdout": native.stdout, "native_stderr": native.stderr,
            "results": rows, "limitations": ["Native stderr field is a Student-t 95% mean half-width, not one SEM.",
            "Pair independence is not established; these intervals do not provide family-aware method-comparison uncertainty.",
            "Raw scores have six-decimal precision; exact full-precision means cannot be reconstructed.",
            "No underlying ontology/annotation score recomputation or prediction completeness validation is claimed."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.repo)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
