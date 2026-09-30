"""Recount retained FAS eligibility and bound hypothetical full-population scores."""

import argparse
from array import array
import gzip
import json
import math
from pathlib import Path
import shlex
import sqlite3
import sys

import ijson
import numpy as np

from benchmark_tools.audit_qfo_fas_samples import read_sample
from benchmark_tools.map_corrected_vgnc_blocks import MANIFEST, MANIFEST_SHA
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

QUERY = ("SELECT DISTINCT p1.uniprot_id, p2.uniprot_id FROM orthologs "
         "JOIN proteomes as p1 ON orthologs.prot_nr1 = p1.prot_nr "
         "JOIN proteomes as p2 ON orthologs.prot_nr2 = p2.prot_nr "
         "WHERE p1.uniprot_id < p2.uniprot_id")


def encode(a, b):
    if not 0 <= a < 2**32 or not 0 <= b < 2**32:
        raise ValueError("Accession index does not fit uint32")
    return (min(a, b) << 32) | max(a, b)


def lookup_index(entries, identifiers):
    raw, codes, scores = array("Q"), array("Q"), array("d")
    scanned = invalid = 0
    for key, directions in entries:
        scanned += 1
        a, b = key.split("_")
        for accession in (a, b):
            if accession not in identifiers:
                identifiers[accession] = len(identifiers)
        ai, bi = identifiers[a], identifiers[b]
        code = encode(ai, bi)
        raw.append((ai << 32) | bi)
        if not isinstance(directions, (list, tuple)) or len(directions) != 2:
            raise ValueError("Unsupported directional array shape")
        try:
            values = [float(v) for v in directions]
        except ValueError:
            invalid += 1
            continue
        if len(values) != 2 or any(not math.isfinite(v) or not 0 <= v <= 1 for v in values):
            raise ValueError("Unsupported directional score; no omission/absence inference")
        codes.append(code)
        scores.append(sum(values) / 2)
        if scanned % 10_000_000 == 0:
            print("lookup entries", scanned, flush=True)
    raw_values = np.frombuffer(raw, dtype=np.uint64)
    raw_values.sort()
    if np.any(raw_values[1:] == raw_values[:-1]):
        raise ValueError("Duplicate JSON keys require native dict-order reconstruction")
    del raw_values, raw
    code_values = np.frombuffer(codes, dtype=np.uint64)
    score_values = np.frombuffer(scores, dtype=np.float64)
    order = np.argsort(code_values, kind="stable")
    ordered_codes = code_values[order]
    last = np.r_[ordered_codes[1:] != ordered_codes[:-1], True] if len(order) else np.array([], dtype=bool)
    selected = order[last]
    return code_values[selected], score_values[selected], dict(
        entries_scanned=scanned, invalid_value_entries=invalid,
        valid_canonical_pairs=len(selected), canonical_overwrites=len(codes) - len(selected))


def match(codes, scores, wanted):
    positions = np.searchsorted(codes, wanted)
    present = positions < len(codes)
    present[present] = codes[positions[present]] == wanted[present]
    return present, scores[positions[present]]


def summarize_population(pairs, identifiers, annotated, codes, scores, batch_size=50000):
    counts = dict(distinct_query_pairs=0, skipped_alias_pairs=0, precomputed=0,
                  missing=0, unannotated=0)
    sums, batch = [], []

    def consume(rows):
        counts["distinct_query_pairs"] += len(rows)
        accepted = [(a, b) for a, b in rows if "_" not in a and "_" not in b]
        counts["skipped_alias_pairs"] += len(rows) - len(accepted)
        for a, b in accepted:
            if not isinstance(a, str) or not isinstance(b, str) or not a < b:
                raise ValueError("Native distinct-query ordering differs")
            for accession in (a, b):
                if accession not in identifiers:
                    identifiers[accession] = len(identifiers)
        wanted = np.array([encode(identifiers[a], identifiers[b]) for a, b in accepted], dtype=np.uint64)
        present, values = match(codes, scores, wanted)
        counts["precomputed"] += int(present.sum())
        sums.append(math.fsum(values))
        for (a, b), exists in zip(accepted, present):
            if not exists:
                counts["missing" if a in annotated and b in annotated else "unannotated"] += 1
        if counts["distinct_query_pairs"] % 5_000_000 == 0:
            print("database pairs", counts["distinct_query_pairs"], flush=True)

    for pair in pairs:
        batch.append(pair)
        if len(batch) == batch_size:
            consume(batch)
            batch = []
    if batch:
        consume(batch)
    pre_sum = math.fsum(sums)
    eligible = counts["precomputed"] + counts["missing"]
    if not eligible:
        raise ValueError("No eligible population")
    return dict(**counts, eligible_pairs=eligible, precomputed_score_sum=pre_sum,
                precomputed_mean=pre_sum / counts["precomputed"] if counts["precomputed"] else None,
                hypothetical_full_mean_bounds=[pre_sum / eligible, (pre_sum + counts["missing"]) / eligible])


def database_pairs(path):
    for suffix in ("-wal", "-shm", "-journal"):
        if Path(str(path) + suffix).exists():
            raise ValueError("Database sidecar needs separate provenance review")
    with sqlite3.connect(path.as_uri() + "?mode=ro", uri=True) as connection:
        connection.execute("PRAGMA query_only=ON")
        connection.execute("PRAGMA temp_store=FILE")
        connection.execute("PRAGMA cache_size=-65536")
        cursor = connection.execute(QUERY)
        while rows := cursor.fetchmany(50000):
            yield from rows


def run(repo, output):
    output.mkdir(parents=True, exist_ok=False)
    manifest_ref = record(repo / MANIFEST)
    if manifest_ref["sha256"] != MANIFEST_SHA:
        raise ValueError("Frozen eight-method comparison changed")
    manifest = json.loads(Path(manifest_ref["path"]).read_text())
    attrition_ref = record(repo / "benchmark_tools/results/qfo_fas_sample_attrition_20260928.json")
    attrition = json.loads(Path(attrition_ref["path"]).read_text())
    review_ref = record(repo / "benchmark_tools/results/fas_native_path_panel_review_20260928.json")
    review = json.loads(Path(review_ref["path"]).read_text())
    check(review["panel"])
    panel = json.loads(Path(review["panel"]["path"]).read_text())
    if (len(manifest["methods"]) != 8 or len(attrition["methods"]) != 8
            or [m["key"] for m in manifest["methods"]] != [m["method"] for m in attrition["methods"]]
            or attrition["status"] != "fas_requested_sample_attrition_bounded"
            or review["status"] != "native_path_panel_record_completeness_verified"):
        raise ValueError("Require complete admitted comparison/sample/annotation bindings")
    refs = [manifest_ref, attrition_ref, review_ref, review["panel"], record(__file__),
            record(Path(__file__).with_name("audit_qfo_fas_samples.py")),
            record(Path(__file__).with_name("prepare_ob_candidate_neighborhood.py")),
            *panel["checked_inputs"]]
    annotated = set()
    for item in panel["files"]:
        check(item["result"])
        refs.append(item["result"])
        annotated.update(r["protein"] for r in json.loads(Path(item["result"]["path"]).read_text())["proteins"])
    if len(annotated) != review["unique_protein_identifiers"]:
        raise ValueError("Annotation universe differs from retained review")
    environment = json.loads(Path(panel["checked_inputs"][0]["path"]).read_text())
    lookup = [r for r in environment["reference_files"] if r["path"].endswith("/fas_precomputed.json.gz")]
    if len(lookup) != 1:
        raise ValueError("Ambiguous frozen lookup")
    refs.append(lookup[0])
    for ref in refs:
        check(ref)
    identifiers = {}
    with gzip.open(lookup[0]["path"], "rb") as stream:
        codes, scores, lookup_summary = lookup_index(ijson.kvitems(stream, "", use_float=True), identifiers)
    (output / "lookup.json").write_text(json.dumps(lookup_summary, indent=2) + "\n")
    rows = []
    for method in attrition["methods"]:
        frozen = next(m for m in manifest["methods"] if m["key"] == method["method"])
        check(frozen["admission"])
        refs.append(frozen["admission"])
        for ref in (method["command"], method["log"], method["raw"]):
            check(ref)
            refs.append(ref)
        tokens = shlex.split(Path(method["command"]["path"]).read_text(), comments=True)
        if tokens.count("--db") != 1 or "--limited-species" in tokens:
            raise ValueError("Unexpected native scoring query scope")
        task = Path(method["command"]["path"]).parent
        db = (task / tokens[tokens.index("--db") + 1]).resolve()
        db_ref = record(db)
        refs.append(db_ref)
        counts = summarize_population(database_pairs(db), identifiers, annotated, codes, scores)
        expected = {key: method["strata"][key]["population"] for key in ("precomputed", "missing")}
        expected["unannotated"] = method["unannotated_pairs_logged"]
        matches = all(counts[key] == value for key, value in expected.items())
        matches = matches and counts["eligible_pairs"] == frozen["details"]["FAS"]["assessed_relations"]
        saved = read_sample(Path(method["raw"]["path"]))
        sample_codes = np.array([encode(identifiers[a], identifiers[b]) for a, b in saved], dtype=np.uint64)
        presence, values = match(codes, scores, sample_codes)
        pre = method["strata"]["precomputed"]["saved"]
        sample_match = (bool(np.all(presence[:pre])) and not bool(np.any(presence[pre:]))
                        and all(math.isclose(a, b, rel_tol=0, abs_tol=1e-12)
                                for a, b in zip(values, list(saved.values())[:pre])))
        row = dict(method=method["method"], database=db_ref, database_historically_hash_bound=False,
                   native_logged_counts=expected, native_counts_match=matches,
                   saved_lookup_strata_and_values_match=sample_match, **counts)
        rows.append(row)
        (output / (method["method"] + ".json")).write_text(json.dumps(row, indent=2) + "\n")
        print(method["method"], "eligible", counts["eligible_pairs"], "log match", matches, flush=True)
        check(db_ref)
        if not matches or not sample_match:
            raise ValueError("Recount differs from native logs/saved sample; retain discrepancy")
    for name in ("ijson", "numpy"):
        refs.extend(record(m.__file__) for key, m in list(sys.modules.items())
                    if (key == name or key.startswith(name + ".")) and getattr(m, "__file__", None))
    for ref in refs:
        check(ref)
    report = dict(status="retained_fas_eligible_populations_recounted", methods=rows,
                  lookup=lookup_summary, checked_records=refs, query=QUERY,
                  versions=dict(numpy=np.__version__, ijson=ijson.__version__, sqlite=sqlite3.sqlite_version),
                  uncertainty_admitted=False, benchmark_scores_changed=False, publication_ready=False,
                  limitations=[
                      "New hashes identify retained database bytes, not independently established historical database pins.",
                      "Counts use the native DISTINCT accession query; correctness of upstream predictions/conversion is not proved.",
                      "Bounds condition on all uncomputed eligible scores having hypothetical values in [0,1]; they are not confidence intervals.",
                      "Precomputed values are reproduced, not validated against underlying annotations or FAS calculations.",
                      "Not a replacement for native sampled FAS, evidence of family generalization, or sampling-law validation."])
    (output / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.repo.resolve(), args.output.resolve())
