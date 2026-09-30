"""Complete the fixed FAS panel by reusing six rows and scanning two databases."""

import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
import shlex
import sqlite3
import subprocess
import sys

import ijson
import numpy as np

from benchmark_tools import audit_fas_population as original
from benchmark_tools.audit_qfo_fas_samples import read_sample
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.render_fas_population import validate, validate_row

PRIOR = "benchmark_tools/results/qfo_fas_population_terminal_22382.json"
PRIOR_SHA = "7734520254562f63dbbf4fed81a2f750fd43530e5d03cd22e85888b467989461"
ORIGINAL_COMMIT = "0591bcbe993130ad214ca21af34400b7ba6a5694"


def unique_records(records):
    result = {}
    for pin in records:
        if pin["path"] in result and result[pin["path"]] != pin:
            raise ValueError("Conflicting evidence identities")
        result[pin["path"]] = pin
    return list(result.values())


def read_pin(pin):
    check(pin)
    result = json.loads(Path(pin["path"]).read_text())
    check(pin)
    return result


def verify_git_pin(pin, repo):
    file_pin = {key: pin[key] for key in ("path", "bytes", "sha256")}
    check(file_pin)
    path = Path(pin["path"]).relative_to(repo)
    blob = subprocess.check_output(["git", "show", f"{ORIGINAL_COMMIT}:{path}"], cwd=repo)
    if hashlib.sha256(blob).hexdigest() != pin["sha256"] or len(blob) != pin["bytes"]:
        raise ValueError("Original producer/input Git binding changed")


def validate_prior(prior, manifest, rows):
    methods = manifest["methods"]
    if (prior.get("schema") != "qfo_fas_population_terminal_partial_v1"
            or prior.get("status") != "timed_out_incomplete_population_recount"
            or prior.get("job_id") != 22382 or prior.get("source_commit") != ORIGINAL_COMMIT
            or prior.get("complete_report_produced") is not False
            or prior.get("uncertainty_admitted") is not False
            or prior.get("benchmark_scores_changed") is not False
            or prior.get("publication_ready") is not False
            or len(methods) != 8 or len(rows) != 6
            or [r["method"] for r in rows] != [m["key"] for m in methods[:6]]
            or [r["method"] for r in prior["completed_partial_rows"]] != [m["key"] for m in methods[:6]]
            or prior["missing_method_rows"] != [m["key"] for m in methods[6:]]):
        raise ValueError("Require the frozen six-row timeout and two missing identities")
    for row, method in zip(rows, methods):
        validate_row(row, method)


def check_sample(saved, precomputed_count, identifiers, codes, scores):
    if type(precomputed_count) is not int or not 0 <= precomputed_count <= len(saved):
        raise ValueError("Invalid saved precomputed stratum")
    for a, b in saved:
        for accession in (a, b):
            if accession not in identifiers:
                identifiers[accession] = len(identifiers)
    wanted = np.array([original.encode(identifiers[a], identifiers[b]) for a, b in saved],
                      dtype=np.uint64)
    present, values = original.match(codes, scores, wanted)
    return (bool(np.all(present[:precomputed_count]))
            and not bool(np.any(present[precomputed_count:]))
            and all(math.isclose(a, b, rel_tol=0, abs_tol=1e-12)
                    for a, b in zip(values, list(saved.values())[:precomputed_count])))


def module_records():
    return [record(module.__file__) for name, module in list(sys.modules.items())
            if (name == "numpy" or name.startswith("numpy.")
                or name == "ijson" or name.startswith("ijson."))
            and getattr(module, "__file__", None)]


def check_database_state(pin):
    for suffix in ("-wal", "-shm", "-journal"):
        sidecar = Path(pin["path"] + suffix)
        if sidecar.exists() or sidecar.is_symlink():
            raise ValueError("Database sidecar prevents recount reuse")
    check(pin)


def run(repo, output):
    if np.__version__ != "2.2.6" or ijson.__version__ != "3.5.0":
        raise ValueError("Declared completion parser/library versions changed")
    output.mkdir(parents=True, exist_ok=False)
    prior_ref = record(repo / PRIOR)
    if prior_ref["sha256"] != PRIOR_SHA:
        raise ValueError("Frozen timeout receipt changed")
    prior = read_pin(prior_ref)
    manifest_ref = record(repo / original.MANIFEST)
    if manifest_ref["sha256"] != original.MANIFEST_SHA:
        raise ValueError("Frozen eight-method manifest changed")
    manifest = read_pin(manifest_ref)
    prior_rows = [read_pin(item["partial_row"]) for item in prior["completed_partial_rows"]]
    validate_prior(prior, manifest, prior_rows)
    for pin in prior["source_git_proof"]:
        verify_git_pin(pin, repo)
    attrition_ref = record(repo / "benchmark_tools/results/qfo_fas_sample_attrition_20260928.json")
    review_ref = record(repo / "benchmark_tools/results/fas_native_path_panel_review_20260928.json")
    attrition, review = read_pin(attrition_ref), read_pin(review_ref)
    for pin in (manifest_ref, attrition_ref, review_ref):
        verify_git_pin(pin, repo)
    panel = read_pin(review["panel"])
    if (len(attrition["methods"]) != 8
            or [m["method"] for m in attrition["methods"]] != [m["key"] for m in manifest["methods"]]
            or attrition["status"] != "fas_requested_sample_attrition_bounded"
            or review["status"] != "native_path_panel_record_completeness_verified"):
        raise ValueError("Incomplete comparison/sample/annotation bindings")
    refs = [prior_ref, *prior["checked_records"], manifest_ref, attrition_ref, review_ref,
            review["panel"], *panel["checked_inputs"], record(__file__),
            record(original.__file__), record(Path(__file__).with_name("render_fas_population.py")),
            record(Path(__file__).with_name("audit_qfo_fas_samples.py")),
            record(Path(__file__).with_name("prepare_ob_candidate_neighborhood.py")),
            *module_records()]
    annotated = set()
    for item in panel["files"]:
        value = read_pin(item["result"])
        refs.append(item["result"])
        annotated.update(p["protein"] for p in value["proteins"])
    if len(annotated) != review["unique_protein_identifiers"]:
        raise ValueError("Annotation universe changed")
    environment = read_pin(panel["checked_inputs"][0])
    lookup_refs = [r for r in environment["reference_files"]
                   if r["path"].endswith("/fas_precomputed.json.gz")]
    if len(lookup_refs) != 1:
        raise ValueError("Ambiguous lookup")
    refs.extend(lookup_refs)
    tasks = []
    for method, frozen in zip(attrition["methods"], manifest["methods"]):
        refs.extend([frozen["admission"], method["command"], method["log"], method["raw"]])
        tokens = shlex.split(Path(method["command"]["path"]).read_text(), comments=True)
        if tokens.count("--db") != 1 or "--limited-species" in tokens:
            raise ValueError("Native query scope changed")
        db = (Path(method["command"]["path"]).parent / tokens[tokens.index("--db") + 1]).resolve()
        db_ref = record(db)
        check_database_state(db_ref)
        tasks.append((method, frozen, db_ref))
        refs.append(db_ref)
    refs = unique_records(refs)
    for pin in refs:
        check(pin)
    for row, (_, _, db_ref) in zip(prior_rows, tasks):
        if row["database"] != db_ref:
            raise ValueError("Reused database differs from native command/current bytes")
    identifiers = {}
    with gzip.open(lookup_refs[0]["path"], "rb") as stream:
        codes, scores, lookup = original.lookup_index(
            ijson.kvitems(stream, "", use_float=True), identifiers)
    if lookup != prior["lookup"]:
        raise ValueError("Complete lookup differs from original pass")
    (output / "lookup.json").write_text(json.dumps(lookup, indent=2) + "\n")
    rows = []
    for index, (method, frozen, db_ref) in enumerate(tasks):
        expected = {k: method["strata"][k]["population"] for k in ("precomputed", "missing")}
        expected["unannotated"] = method["unannotated_pairs_logged"]
        if index < 6:
            row = dict(prior_rows[index], reused_prior_recount=True,
                       prior_row_source=prior["completed_partial_rows"][index]["partial_row"])
        else:
            counts = original.summarize_population(original.database_pairs(Path(db_ref["path"])),
                                                   identifiers, annotated, codes, scores)
            row = dict(method=method["method"], database=db_ref,
                       database_historically_hash_bound=False, reused_prior_recount=False,
                       native_logged_counts=expected, native_counts_match=all(
                           counts[k] == v for k, v in expected.items()), **counts)
        saved = read_sample(Path(method["raw"]["path"]))
        row["saved_lookup_strata_and_values_match"] = check_sample(
            saved, method["strata"]["precomputed"]["saved"], identifiers, codes, scores)
        if row["native_logged_counts"] != expected:
            raise ValueError("Reused counts differ from native logs")
        validate_row(row, frozen)
        check_database_state(db_ref)
        rows.append(row)
        (output / (method["method"] + ".json")).write_text(json.dumps(row, indent=2) + "\n")
        print(method["method"], "reused" if index < 6 else "recounted",
              "eligible", row["eligible_pairs"], flush=True)
    refs = unique_records([*refs, *module_records()])
    for pin in refs:
        check(pin)
    for _, _, db_ref in tasks:
        check_database_state(db_ref)
    report = dict(status="retained_fas_eligible_populations_completed_with_reuse", methods=rows,
        lookup=lookup, query=original.QUERY, checked_records=refs,
        reuse=dict(prior_receipt=prior_ref, original_job_id=22382, original_state="TIMEOUT",
                   original_final_stability_pass_completed=False,
                   historical_parser_hash_identity_established=False,
                   completion_input_stability_pass_completed=True,
                   fresh_database_recounts=2, reused_database_recounts=6),
        versions=dict(numpy=np.__version__, ijson=ijson.__version__, sqlite=sqlite3.sqlite_version),
        benchmark_scores_changed=False, uncertainty_admitted=False, publication_ready=False,
        limitations=["Six fixed rows are reused; this is not a fresh eight-database pass or repair of job 22382.",
                     "Current before/after stability does not establish the historical parser file identity.",
                     "Database hashes are retained identities, not historical prediction/conversion proof.",
                     "Bounds assume hypothetical uncomputed scores in [0,1]; not native FAS or confidence intervals.",
                     "Lookup/sampling validity and biological generalization remain separate."])
    validate(report, manifest)
    (output / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.repo.resolve(), args.output.resolve())
