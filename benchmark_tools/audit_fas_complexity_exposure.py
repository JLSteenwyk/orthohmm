"""Intersect saved FAS strata with native complexity flags, without causal attribution."""

import argparse
from collections import Counter
import gzip
import json
import math
from pathlib import Path
import sys

import ijson

from benchmark_tools.audit_qfo_fas_samples import read_sample
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


def summarize(pairs, rejected, known):
    counts = Counter()
    sums = Counter()
    for (a, b), value in pairs:
        if a not in known or b not in known:
            raise ValueError("Saved score endpoint absent from annotation panel")
        if not math.isfinite(value) or not 0 <= value <= 1:
            raise ValueError("Invalid saved score")
        category = "touches_flagged_protein" if a in rejected or b in rejected else "neither_flagged"
        counts[category] += 1
        sums[category] += value
    return {key: dict(pairs=counts[key], mean=sums[key] / counts[key] if counts[key] else None)
            for key in ("touches_flagged_protein", "neither_flagged")}


def verify_lookup(wanted, entries):
    seen_keys, observed = set(), {}
    scanned = collisions = 0
    for key, directions in entries:
        scanned += 1
        a, b = key.split("_")
        pair = tuple(sorted((a, b)))
        if pair not in wanted:
            continue
        if key in seen_keys:
            raise ValueError("Duplicate relevant JSON key")
        seen_keys.add(key)
        try:
            scores = [float(v) for v in directions]
        except (ValueError, TypeError):
            continue
        if len(scores) != 2 or any(not math.isfinite(v) or not 0 <= v <= 1 for v in scores):
            raise ValueError("Unexpected precomputed directional scores")
        collisions += pair in observed
        observed[pair] = sum(scores) / 2
        if wanted[pair] is None:
            raise ValueError("Logged new score is present in precomputed lookup")
    expected = {pair for pair, value in wanted.items() if value is not None}
    if set(observed) != expected:
        raise ValueError("Saved precomputed pairs missing from lookup")
    if any(not math.isclose(value, wanted[pair], rel_tol=0, abs_tol=1e-12)
           for pair, value in observed.items()):
        raise ValueError("Saved precomputed score differs from lookup")
    return dict(entries_scanned=scanned, unique_saved_precomputed_pairs=len(observed),
                unique_saved_new_pairs=len(wanted) - len(observed),
                relevant_canonical_overwrites=collisions)


def audit(repo):
    results = repo / "benchmark_tools/results"
    panel_ref = record(results / "fas_native_path_panel_20260928.json")
    review_ref = record(results / "fas_native_path_panel_review_20260928.json")
    attrition_ref = record(results / "qfo_fas_sample_attrition_20260928.json")
    panel = json.loads(Path(panel_ref["path"]).read_text())
    review = json.loads(Path(review_ref["path"]).read_text())
    attrition = json.loads(Path(attrition_ref["path"]).read_text())
    if review["panel"] != panel_ref or review["status"] != "native_path_panel_record_completeness_verified":
        raise ValueError("Unreviewed annotation panel")
    if attrition["status"] != "fas_requested_sample_attrition_bounded" or len(attrition["methods"]) != 8:
        raise ValueError("Require eight-method corrected attrition audit")
    refs = [panel_ref, review_ref, attrition_ref, record(__file__),
            record(Path(__file__).with_name("audit_qfo_fas_samples.py")),
            *panel["checked_inputs"], *attrition["checked_records"]]
    for ref in refs:
        check(ref)
    known, rejected = set(), set()
    for item in panel["files"]:
        check(item["result"])
        refs.append(item["result"])
        rows = json.loads(Path(item["result"]["path"]).read_text())["proteins"]
        ids = {r["protein"] for r in rows}
        if len(ids) != len(rows) or known & ids:
            raise ValueError("Duplicate annotation identifiers")
        known.update(ids)
        rejected.update(r["protein"] for r in rows if r["rejected_by_path_limit"])
    if len(known) != review["unique_protein_identifiers"] or len(rejected) != review["unique_rejected_identifiers"]:
        raise ValueError("Panel counts differ from review")
    methods, wanted = [], {}
    for method in attrition["methods"]:
        check(method["raw"])
        pairs = list(read_sample(Path(method["raw"]["path"])).items())
        pre = method["strata"]["precomputed"]["saved"]
        new = method["strata"]["missing"]["saved"]
        if len(pairs) != pre + new:
            raise ValueError("Saved sample count differs from attrition audit")
        for index, (pair, score) in enumerate(pairs):
            expected = score if index < pre else None
            if pair in wanted and wanted[pair] != expected:
                raise ValueError("Shared pair has inconsistent stratum or precomputed score")
            wanted[pair] = expected
        methods.append(dict(method=method["method"], raw=method["raw"],
            precomputed=summarize(pairs[:pre], rejected, known),
            newly_calculated=summarize(pairs[pre:], rejected, known)))
    environment = json.loads(Path(panel["checked_inputs"][0]["path"]).read_text())
    candidates = [r for r in environment["reference_files"] if r["path"].endswith("/fas_precomputed.json.gz")]
    if len(candidates) != 1:
        raise ValueError("Ambiguous frozen precomputed lookup")
    lookup = candidates[0]
    check(lookup)
    refs.append(lookup)
    with gzip.open(lookup["path"], "rb") as stream:
        lookup_check = verify_lookup(wanted, ijson.kvitems(stream, "", use_float=True))
    parser_sources = sorted({m.__file__ for name, m in sys.modules.items()
                             if name.startswith("ijson") and getattr(m, "__file__", None)})
    refs.extend(record(p) for p in parser_sources)
    for ref in refs:
        check(ref)
    return dict(status="saved_fas_complexity_exposure_checked", checked_inputs=refs,
        annotation_identifiers=len(known), flagged_identifiers=len(rejected), methods=methods,
        precomputed_lookup_check=lookup_check, ijson_version=ijson.__version__,
        historical_omissions_attributed=False, uncertainty_admitted=False,
        limitations=["Stratum assignments and saved precomputed means are independently checked against the frozen lookup; underlying FAS score correctness is not established.",
            "Saved-score exposure does not identify omitted sampled pairs or estimate eligible-population exposure.",
            "Precomputed scores are not invalidated by the current new-scoring complexity cutoff.",
            "Conditional subgroup means are descriptive, not replacement endpoints or causal effects."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    value = audit(args.repo.resolve())
    with args.output.open("x") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
