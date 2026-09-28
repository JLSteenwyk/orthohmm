"""Bind corrected FAS samples to execution outputs and reproduce their statistics."""

import argparse
import json
import math
from pathlib import Path

from benchmark_tools.audit_qfo_fas_samples import read_sample, summarize, render_table
from benchmark_tools.map_corrected_vgnc_blocks import MANIFEST, MANIFEST_SHA
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_qfo_scored_pair_panel import execution_binding


def audit(repo):
    identity = record(repo / MANIFEST)
    if identity["sha256"] != MANIFEST_SHA:
        raise ValueError("Changed corrected comparison manifest")
    manifest = json.loads(Path(identity["path"]).read_text())
    methods = manifest["methods"]
    if (len(methods) != 8 or len({m["key"] for m in methods}) != 8
            or any(m["status"] != "admitted" for m in methods)):
        raise ValueError("Require eight distinct admitted methods")
    checked = [identity, record(__file__),
               record(Path(__file__).with_name("audit_qfo_fas_samples.py")),
               record(Path(__file__).with_name("run_qfo_scored_pair_panel.py"))]
    rows = []
    for method in methods:
        check(method["admission"])
        checked.append(method["admission"])
        admission = json.loads(Path(method["admission"]["path"]).read_text())
        execution_ref, execution = execution_binding(admission)
        checked.append(execution_ref)
        candidates = [r for r in admission["metric_files"]
                      if r["path"].endswith("/results/FAS/FAS.json")]
        if len(candidates) != 1:
            raise ValueError("Ambiguous admitted FAS endpoint")
        endpoint = candidates[0]
        check(endpoint)
        checked.append(endpoint)
        aggregate = json.loads(Path(endpoint["path"]).read_text())["datalink"]["inline_data"]
        if (aggregate["visualization"]["x_axis"] != "NR_ORTHOLOGS"
                or aggregate["visualization"]["y_axis"] != "FAS"):
            raise ValueError("Unexpected FAS axes")
        participants = aggregate["challenge_participants"]
        paths = list(Path(endpoint["path"]).parent.glob("*raw.txt.gz"))
        if len(participants) != 1 or len(paths) != 1:
            raise ValueError("Ambiguous participant or raw FAS file")
        participant = participants[0]
        expected_participant = method.get("participant")
        if expected_participant is None:
            expected_participant = admission["assessment"]["participant"]
        if not expected_participant or participant.get("participant_id") != expected_participant:
            raise ValueError("FAS participant differs from corrected comparison")
        if (not math.isclose(participant["metric_y"], method["scores"]["FAS"],
                             rel_tol=0, abs_tol=1e-12)
                or participant["metric_x"] != method["details"]["FAS"]["assessed_relations"]):
            raise ValueError("FAS endpoint differs from corrected comparison")
        raw = record(paths[0])
        historical = [r for r in execution["outputs"] if r["path"] == raw["path"]]
        if len(historical) != 1 or historical[0] != raw:
            raise ValueError("Raw FAS table differs from historical execution output pin")
        checked.append(raw)
        rows.append(dict(method=method["key"], aggregate=endpoint, raw=raw,
            historical_execution_output_bound=True, **summarize(read_sample(paths[0]), participant)))
    for ref in checked:
        check(ref)
    return dict(status="corrected_fas_sample_arithmetic_verified", methods=rows,
        checked_records=checked, uncertainty_admitted=False, publication_ready=False,
        limitations=[
            "Eligible-pair counts are admitted endpoint values, not independently recounted predictions.",
            "Saved-sample arithmetic does not establish annotation or FAS implementation correctness.",
            "Method-specific sampling fractions and missing-score handling prevent assuming representative common samples.",
            "Native pair-IID SEM is reproduced, not endorsed as family-aware uncertainty or a paired comparison interval.",
            "Historical random states and stratum-specific inclusion probabilities are not reconstructed."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--table", type=Path, required=True)
    args = parser.parse_args()
    if any(p.exists() or p.is_symlink() for p in (args.output, args.table)):
        raise FileExistsError("Outputs must be fresh")
    result = audit(args.repo.resolve())
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    with args.table.open("x") as stream:
        stream.write(render_table(result).replace("# Retained FAS Samples", "# Corrected FAS Samples")
                     .replace("audit_qfo_fas_samples.py", "audit_corrected_fas_samples.py"))
