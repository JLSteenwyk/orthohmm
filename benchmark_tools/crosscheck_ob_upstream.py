"""Cross-check comparator readback with the retained upstream evaluator functions."""

import argparse
from contextlib import redirect_stdout
import importlib.util
import io
import json
import math
from pathlib import Path

from benchmark_tools.audit_failed_recovery_refinement import record
from benchmark_tools.normalize_three_kingdoms_orthogroups import iter_fastoma

UPSTREAM_SHA = "81eb1e660c17819549b07eea8a54b4fb42a89180cafeb4569d92195d282f5e6f"
READBACK_SHA = "5dcb4e65f277eb9debf2bbacb1d4e298c678ac13fc981e75c4ba78e4851eeb55"


def compare(values, expected):
    if len(values) != 3 or any(not math.isfinite(x) for x in values):
        raise ValueError("Invalid upstream metrics")
    differences = {key: observed - expected[key] for key, observed in
                   zip(("f_score", "precision", "recall"), values)}
    if any(abs(value) > 1e-10 for value in differences.values()):
        raise ValueError("Upstream/reimplementation metric disagreement")
    return differences


def audit(root, evaluator):
    source = root / "benchmark_tools/results/retained_ob_comparator_readback_20260926.json"
    inputs = [record(source), record(evaluator)]
    if [r["sha256"] for r in inputs] != [READBACK_SHA, UPSTREAM_SHA]:
        raise ValueError("Changed retained readback or upstream evaluator")
    prior = json.loads(source.read_text())
    if prior["all_scores_agree"] is not True or len(prior["rows"]) != 5:
        raise ValueError("Require complete comparator readback")
    inputs.extend(prior["checked_records"])
    refbase = evaluator.parent / "RefOGs"
    retained_refs = {item["path"] for item in prior["checked_records"]
                     if Path(item["path"]).name.startswith("RefOG")}
    actual_refs = {str(path) for path in refbase.glob("RefOG*.txt")} | {
        str(path) for path in (refbase / "low_certainty_assignments").glob("RefOG*.txt")}
    if actual_refs != retained_refs:
        raise ValueError("Upstream reference directory differs from pinned readback")
    fasta = sorted((evaluator.parent / "Input").glob("*fa"))
    if len(fasta) != 12:
        raise ValueError("Changed upstream input inventory")
    inputs.extend(record(path) for path in fasta)
    inputs.extend(record(Path(__file__).with_name(name).resolve()) for name in
                  ("crosscheck_ob_upstream.py", "normalize_three_kingdoms_orthogroups.py"))
    for item in inputs:
        if record(item["path"]) != item:
            raise ValueError("Changed scoring input")
    spec = importlib.util.spec_from_file_location("retained_upstream_orthobench", evaluator)
    upstream = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(upstream)
    refdir = str(evaluator.parent / "RefOGs") + "/"
    references = upstream.read_refogs(refdir)
    uncertain = upstream.read_uncertain_refogs(refdir + "low_certainty_assignments/")
    genes = upstream.get_expected_genes()
    rows = []
    for row in prior["rows"]:
        path = Path(row["prediction"]["path"])
        groups = ([set(members) for _, members in iter_fastoma(path)] if row["key"] == "fastoma_0_3_5"
                  else upstream.read_orthogroups_smart(str(path)))
        log = io.StringIO()
        with redirect_stdout(log):
            upstream.check_orthogroups(groups, genes)
            values = upstream.calculate_benchmarks_pairwise(references, uncertain, groups)
        text = log.getvalue()
        if "ERROR:" in text:
            raise ValueError("Upstream validation reported an error: " + text)
        rows.append(dict(key=row["key"], prediction=row["prediction"], upstream_scores=values,
                         differences=compare(values, row["score"]), stdout=text,
                         reader="iter_fastoma" if row["key"] == "fastoma_0_3_5" else "upstream.read_orthogroups_smart"))
    if any(record(item["path"]) != item for item in inputs):
        raise ValueError("Evidence changed during upstream cross-check")
    return dict(status="upstream_ob_comparator_functions_agree", rows=rows, checked_records=inputs,
                expected_input_genes=len(genes), publication_ready=False,
                limitations=["Upstream evaluator functions executed, not its CLI wrapper or native inference.",
                    "FastOMA retains the shared two-column adapter; parser independence is not claimed there.",
                    "Current input membership checks do not establish historical consumption or timing.",
                    "Missing native genes are retained; upstream warning text is not suppressed."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--evaluator", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    report = audit(args.root.resolve(), args.evaluator.resolve())
    with args.output.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
