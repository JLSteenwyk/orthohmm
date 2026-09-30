"""Independently check retained FAS completion arithmetic and reuse provenance."""

import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import shlex
import subprocess


COMMIT = "bd9863bb98f463f7aba1a34a96e0982e83c4a71b"
PRIOR_SHA = "7734520254562f63dbbf4fed81a2f750fd43530e5d03cd22e85888b467989461"
MANIFEST_SHA = "042aa221d01554114a1b0b413ca3e2ff56dad5f8c06362c69a1bf49f2de642fc"
QUERY = ("SELECT DISTINCT p1.uniprot_id, p2.uniprot_id FROM orthologs "
         "JOIN proteomes as p1 ON orthologs.prot_nr1 = p1.prot_nr "
         "JOIN proteomes as p2 ON orthologs.prot_nr2 = p2.prot_nr "
         "WHERE p1.uniprot_id < p2.uniprot_id")


def pin(path):
    path = Path(path).resolve(strict=True)
    digest, size = hashlib.sha256(), 0
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            size += len(chunk)
            digest.update(chunk)
    return dict(path=str(path), bytes=size, sha256=digest.hexdigest())


def check(ref):
    if pin(ref["path"]) != ref:
        raise ValueError("Evidence identity changed")


def read(ref):
    check(ref)
    value = json.loads(Path(ref["path"]).read_text())
    check(ref)
    return value


def terminal_fields(raw):
    if len(raw.splitlines()) != 1:
        raise ValueError("Require one retained controller record")
    pairs = [token.split("=", 1) for token in shlex.split(raw) if "=" in token]
    fields = dict(pairs)
    expected = dict(JobId="22383", JobState="COMPLETED", ExitCode="0:0",
                    Restarts="0", Requeue="0", NodeList="bizon", NumCPUs="1",
                    NumTasks="1", MinMemoryNode="16G", TimeLimit="03:00:00")
    if (len(fields) != len(pairs) or any(fields.get(k) != v for k, v in expected.items())
            or "ArrayJobId" in fields or "ArrayTaskId" in fields):
        raise ValueError("Require successful original completion allocation")
    return fields


def readback(report, manifest, attrition, prior, partials):
    methods, rows = manifest["methods"], report["methods"]
    if (len(methods) != 8 or len(rows) != 8 or len(partials) != 6
            or len(attrition["methods"]) != 8
            or [r["method"] for r in rows] != [m["key"] for m in methods]
            or len({m["key"] for m in methods}) != 8
            or [r["method"] for r in attrition["methods"]] != [m["key"] for m in methods]
            or report["status"] != "retained_fas_eligible_populations_completed_with_reuse"
            or report["query"] != QUERY or report["lookup"] != prior["lookup"]
            or prior["status"] != "timed_out_incomplete_population_recount"
            or prior["job_id"] != 22382):
        raise ValueError("Incomplete or changed fixed comparison")
    for key in ("uncertainty_admitted", "benchmark_scores_changed", "publication_ready"):
        if report[key] is not False or prior[key] is not False:
            raise ValueError("Unsupported scientific admission")
    reuse = report["reuse"]
    expected_reuse = dict(prior_receipt=reuse["prior_receipt"], original_job_id=22382,
        original_state="TIMEOUT", original_final_stability_pass_completed=False,
        historical_parser_hash_identity_established=False,
        completion_input_stability_pass_completed=True,
        fresh_database_recounts=2, reused_database_recounts=6)
    if (reuse != expected_reuse or reuse["prior_receipt"]["sha256"] != PRIOR_SHA
            or any(type(reuse[k]) is not type(v) for k, v in expected_reuse.items())):
        raise ValueError("Changed reuse or historical provenance claim")
    result = []
    for index, (row, method, native) in enumerate(zip(rows, methods, attrition["methods"])):
        if (row["native_counts_match"] is not True
                or row["saved_lookup_strata_and_values_match"] is not True
                or row["database_historically_hash_bound"] is not False
                or row["reused_prior_recount"] is not (index < 6)):
            raise ValueError("Invalid row scope or reported checks")
        if index < 6:
            source = prior["completed_partial_rows"][index]["partial_row"]
            original = {k: v for k, v in row.items()
                        if k not in {"reused_prior_recount", "prior_row_source"}}
            if row["prior_row_source"] != source or original != partials[index]:
                raise ValueError("Reused row differs from retained original")
        elif "prior_row_source" in row:
            raise ValueError("Fresh row claims partial provenance")
        fields = ("distinct_query_pairs", "skipped_alias_pairs", "precomputed",
                  "missing", "unannotated", "eligible_pairs")
        if any(type(row[k]) is not int or row[k] < 0 for k in fields):
            raise ValueError("Noninteger population count")
        p, m, n = row["precomputed"], row["missing"], row["eligible_pairs"]
        s = row["precomputed_score_sum"]
        logged = dict(precomputed=native["strata"]["precomputed"]["population"],
                      missing=native["strata"]["missing"]["population"],
                      unannotated=native["unannotated_pairs_logged"])
        if (n <= 0 or p <= 0 or p + m != n
                or n != method["details"]["FAS"]["assessed_relations"]
                or row["distinct_query_pairs"] != n + row["unannotated"] + row["skipped_alias_pairs"]
                or row["native_logged_counts"] != logged
                or any(row[k] != v for k, v in logged.items())
                or type(s) not in (int, float) or not math.isfinite(s) or not 0 <= s <= p):
            raise ValueError("Population conservation or native count mismatch")
        # Rational arithmetic checks the stored aggregate, not a repeated lookup scan.
        total = Fraction(s)
        bounds = [float(total / n), float((total + m) / n)]
        mean = float(total / p)
        observed = row["hypothetical_full_mean_bounds"]
        if (len(observed) != 2 or any(type(v) not in (int, float) or
                not math.isclose(v, b, rel_tol=0, abs_tol=1e-15) for v, b in zip(observed, bounds))
                or not math.isclose(row["precomputed_mean"], mean, rel_tol=0, abs_tol=1e-15)):
            raise ValueError("Stored mean or completion-bound arithmetic differs")
        result.append(dict(method=row["method"], eligible_pairs=n, precomputed=p,
                           missing=m, missing_fraction=m / n, precomputed_mean=mean,
                           hypothetical_full_mean_bounds=bounds,
                           reused_prior_recount=index < 6))
    return result


def run(repo, report_path, report_sha, controller, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    report_ref, terminal_ref, source_ref = pin(report_path), pin(controller), pin(__file__)
    if report_ref["sha256"] != report_sha:
        raise ValueError("Explicit report digest differs")
    terminal_fields(Path(controller).read_text())
    report = read(report_ref)
    refs = {r["path"]: r for r in report["checked_records"]}
    if len(refs) != len(report["checked_records"]):
        raise ValueError("Duplicate report evidence identity")
    for ref in refs.values():
        check(ref)
    manifest_ref = refs[str(repo / "benchmark_tools/results/qfo_corrected_comparison_20260926_v7/manifest.json")]
    if manifest_ref["sha256"] != MANIFEST_SHA:
        raise ValueError("Frozen manifest differs")
    prior_ref = report["reuse"]["prior_receipt"]
    prior = read(prior_ref)
    attrition = read(refs[str(repo / "benchmark_tools/results/qfo_fas_sample_attrition_20260928.json")])
    manifest = read(manifest_ref)
    partials = [read(item["partial_row"]) for item in prior["completed_partial_rows"]]
    submission_ref = pin(repo / "benchmark_tools/results/qfo_fas_population_completion_submission_20260930.json")
    submission = read(submission_ref)
    if submission["job_id"] != 22383 or submission["source_commit"] != COMMIT or submission["returncode"] != 0:
        raise ValueError("Submission binding differs")
    for ref in submission["sources"]:
        check(ref)
        path = Path(ref["path"]).relative_to(repo)
        blob = subprocess.check_output(["git", "show", f"{COMMIT}:{path}"], cwd=repo)
        if len(blob) != ref["bytes"] or hashlib.sha256(blob).hexdigest() != ref["sha256"]:
            raise ValueError("Submitted source Git binding differs")
    rows = readback(report, manifest, attrition, prior, partials)
    for row, native in zip(report["methods"], attrition["methods"]):
        database = row["database"]
        if refs[database["path"]] != database:
            raise ValueError("Row database missing from stable inputs")
        for suffix in ("-wal", "-shm", "-journal"):
            sidecar = Path(database["path"] + suffix)
            if sidecar.exists() or sidecar.is_symlink():
                raise ValueError("Database sidecar present")
        command = native["command"]
        check(command)
        tokens = shlex.split(Path(command["path"]).read_text(), comments=True)
        if tokens.count("--db") != 1 or "--limited-species" in tokens:
            raise ValueError("Changed native query scope")
        db = (Path(command["path"]).parent / tokens[tokens.index("--db") + 1]).resolve()
        if str(db) != database["path"]:
            raise ValueError("Native command database differs")
    for ref in (report_ref, terminal_ref, submission_ref, source_ref):
        check(ref)
    result = dict(schema="qfo_fas_completion_independent_readback_v1",
        status="eight_method_aggregate_and_reuse_readback_passed", job_id=22383,
        report=report_ref, controller=terminal_ref, submission=submission_ref,
        readback_source=source_ref,
        checked_report_records=len(refs), submitted_git_sources_checked=len(submission["sources"]),
        methods=rows, joins_repeated=False, lookup_rescored=False,
        historical_parser_identity_established=False, uncertainty_admitted=False,
        benchmark_scores_changed=False, publication_ready=False)
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path("."))
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--report-sha256", required=True)
    parser.add_argument("--controller", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.repo.resolve(), args.report, args.report_sha256, args.controller, args.output)
