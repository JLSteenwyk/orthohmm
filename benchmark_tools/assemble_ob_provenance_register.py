"""Join retained OrthoBench scores with explicit provenance/resource limits."""

import argparse
import json
from pathlib import Path
import re

from benchmark_tools.prepare_ob_candidate_neighborhood import record,check

SOURCES = {
    "scores":("current_benchmark_scores_20260926_v2/manifest.json","8d56dae74721c4530ab2b86c5008b2d0ebac07878e574fd85b1ad84e4b908f07"),
    "readback":("retained_ob_comparator_readback_20260926.json","5dcb4e65f277eb9debf2bbacb1d4e298c678ac13fc981e75c4ba78e4851eeb55"),
    "of":("ob_orthofinder_provenance_20260926.json","509195fa01994a148ca0af5403be2a64ba702028b753b5172492505ed3126682"),
    "matrix":("ob_matrix_provenance_20260926.json","f68706495773a015518ef0bea2be73090da79c5baf1a9290a81f7aad95521fb3"),
    "fastoma":("ob_fastoma_provenance_20260926.json","516a7d7dc125309e2dbfd2dab89eec7af346b12bf961ac8c3d5e965d98a768b3"),
    "orthomcl":("ob_orthomcl_provenance_20260926.json","cbbadd37dc35dcf41d07a5b7ce29269dc9c4b0a49ebb0180f16d879bf33182d5"),
}


def workflow_duration(lines):
    values = set()
    for line in lines:
        if line.strip().startswith("Duration"):
            match = re.fullmatch(r"Duration\s*:\s*(\d+)m (\d+)s",line.strip())
            if not match or int(match[2]) >= 60:
                raise ValueError("Unsupported workflow duration format")
            values.add(int(match[1])*60+int(match[2]))
    if len(values) != 1:
        raise ValueError("Missing or conflicting workflow durations")
    return values.pop()


def details(key, reports):
    row = dict(input_sequences_equal=None,input_note="Not consolidated here; retained replay evidence is separate",
        output_semantics="refined orthogroups",inference_wall_seconds=None,conversion_wall_seconds=None,
        timing_basis="No consolidated end-to-end time for this historical score row",
        controlled_timing=False,peak_process_rss_kib=None,command=None,limitations=[])
    if key.startswith("orthohmm_"):
        if key == "orthohmm_phylogeny_satellite_v2":
            row["output_semantics"] = "post-membership-filter root HOGs"
        row["limitations"] = ["Historical score, not the fresh installed reproduction or cached incremental duration"]
    elif key.startswith("orthofinder_"):
        source = reports["of"]
        row.update(input_sequences_equal=source["checks"]["processed_sequences_match"],
                   input_note="Staged bytes and processed protein sequences match",command=source["command"])
        if key == "orthofinder_3_1_5_full":
            resources = source["full_inference_resources"]
            row.update(output_semantics="final root HOGs",inference_wall_seconds=resources["elapsed_seconds"],
                peak_process_rss_kib=resources["max_process_rss_kib"],resources=resources,
                timing_basis="GNU-time full native run; shared host")
        else:
            row.update(output_semantics="pre-phylogeny clustering checkpoint",command=source["checkpoint_conversion_command"],
                conversion_wall_seconds=source["checkpoint_conversion_resources"]["elapsed_seconds"],
                timing_basis="Conversion only; sequence-only inference time unknown")
    elif key in {"sonicparanoid_2_0_9","proteinortho_6_3_6"}:
        source = next(r for r in reports["matrix"]["rows"] if r["key"] == key)
        native = source["native_metadata"]
        row.update(input_sequences_equal=source["input_hashes_match"],output_semantics="native groups plus explicit singleton padding",
            native_assigned_genes=source["native_assigned_genes"],singleton_padding=source["singleton_padding"],
            inference_wall_seconds=native["tool_reported_elapsed_seconds"],native_metadata=native,
            input_note="Staged bytes match" if source["input_hashes_match"] else "177 proteins differ; 869 asterisks deleted",
            timing_basis="Tool-reported run duration" if native["tool_reported_elapsed_seconds"] is not None else "No verified per-invocation duration")
        row["limitations"] = (["Exact argv, CPU time and RSS unavailable"] if key.startswith("sonic") else
                              ["Two retained invocations; accuracy effect of input sanitization not measured"])
    elif key == "fastoma_0_3_5":
        source = reports["fastoma"]
        row.update(input_sequences_equal=source["all_protein_sequences_match"],input_note="IDs and residues match; renamed/reordered files",
            output_semantics="final orthologous groups, not root-HOG diagnostic",
            inference_wall_seconds=workflow_duration(source["retained_run_summary"]),
            timing_basis="Nextflow workflow duration including failed attempts",native_metadata=source["retained_run_summary"],
            task_exit_counts=source["task_exit_counts"],native_assigned_genes=source["final_groups"]["assigned_genes"],
            limitations=["Supplied species tree", "Ten retained exit-137 tasks; final scored collection succeeded",
                         "Historical image/database identities and peak RSS unresolved"])
    elif key == "orthomcl_1_4":
        source = next(r for r in reports["orthomcl"]["runs"] if r["run"] == "Jul_25")
        row.update(input_sequences_equal=source["sequences_exact"],input_note="177 proteins differ; 869 asterisks deleted",
            output_semantics="July 2025 native OrthoMCL groups",native_assigned_genes=source["assigned_genes"],
            inference_wall_seconds=source["logged_duration"]["seconds"],command=source["commands"],
            timing_basis="July 2025 native-log interval; 8 BLAST threads",run="Jul_25",
            limitations=["Do not substitute April 2026 input identity or 32-thread timing", "CPU time and RSS unknown"])
    else:
        raise ValueError(f"Unknown method {key}")
    return row


def supplement_orthohmm(rows, audit):
    if (audit["status"] != "retained_orthohmm_ob_provenance_audited"
            or audit["controlled_comparative_resources"] is not False
            or audit["complete_transitive_provenance"] is not False):
        raise ValueError("Unexpected historical audit scope")
    selected = {r["key"]: r for r in rows if r["key"].startswith("orthohmm_")}
    names = {"orthohmm_high_sensitivity": "high_sensitivity",
             "orthohmm_phylogeny_satellite_v2": "phylogeny"}
    if len(selected) != 2 or set(selected) != set(names):
        raise ValueError("Expected both historical OrthoHMM rows")
    for key, name in names.items():
        if selected[key]["prediction"] != audit[name]["prediction"]:
            raise ValueError("Historical prediction binding differs")
    high, phylo = audit["high_sensitivity"], audit["phylogeny"]
    selected["orthohmm_high_sensitivity"].update(command=high["command"],
        input_note="Cached-hit identity verified; historical FASTA bytes not proven",
        inference_wall_seconds=None, replay_wall_seconds=high["wall_s"],
        timing_basis="Cached-hit downstream replay only; initial search excluded",
        resources={"peak_process_rss_gib": high["peak_process_rss_gib"],
                   "scope": high["scope"], "timings": high["timings"]},
        historical_audit=high,
        limitations=["Replay time is not full inference time", "Historical FASTA checksums unavailable"])
    selected["orthohmm_phylogeny_satellite_v2"].update(command=phylo["command"],
        input_sequences_equal=True, input_note="All 12 historical input hashes match retained FASTAs",
        inference_wall_seconds=phylo["resources"]["wall_s"], resources=phylo["resources"],
        timing_basis="Historical full inference, shared host; scoring excluded",
        historical_audit=phylo,
        limitations=["Retained metrics are not immutable execution attestation",
                     "Sampled summed RSS is not unique physical memory; external runtime audit incomplete"])


def assemble(repo, output, include_orthohmm_audit=False):
    if output.exists():
        raise FileExistsError(output)
    reports, records = {}, []
    for key,(name,checksum) in SOURCES.items():
        path = repo / "benchmark_tools/results" / name
        item = record(path)
        if item["sha256"] != checksum:
            raise ValueError(f"Changed source {key}")
        reports[key] = json.loads(path.read_text())
        records.append(item)
    readback = {r["key"]:r["prediction"] for r in reports["readback"]["rows"]}
    rows = []
    for source in reports["scores"]["rows"]:
        prediction = readback.get(source["key"],source["orthobench_retained_evidence"].get("prediction_provenance"))
        if prediction is None:
            raise ValueError("Missing score prediction identity")
        check(prediction)
        rows.append(dict(key=source["key"],label=source["label"],weighted_refog_f1=source["scores"]["OrthoBench"],
            prediction=prediction,**details(source["key"],reports)))
    if len(rows) != 8 or len({r["key"] for r in rows}) != 8:
        raise ValueError("Expected exactly eight score rows")
    selected = next(r for r in rows if r["key"] == "orthomcl_1_4")
    july = next(r for r in reports["orthomcl"]["runs"] if r["run"] == "Jul_25")
    if selected["prediction"] != july["input_records"][0]:
        raise ValueError("OrthoMCL row would mix run identities")
    if include_orthohmm_audit:
        audit_path = repo / "benchmark_tools/results/ob_orthohmm_retained_provenance_20260927.json"
        item = record(audit_path)
        if item["sha256"] != "39f492ddfb05166c805e3c420dd98fcbef5dcffb31a8a367b1f3ee17774c759b":
            raise ValueError("Changed historical OrthoHMM audit")
        audit = json.loads(audit_path.read_text())
        records.extend([item, *audit["checked_records"]])
        for checked in records:
            check(checked)
        supplement_orthohmm(rows, audit)
    for item in records:
        check(item)
    output.mkdir(parents=True)
    result = dict(status="orthobench_provenance_register_partial",rows=rows,inputs=records,source=record(__file__),
        complete_transitive_provenance=False,controlled_comparative_resources=False,scores_recomputed=False,
        limitations=["Consolidates existing retained audits; does not independently attest historical execution.",
            "Unknown is not zero; process RSS is not aggregate concurrent-process memory.",
            "Mixed run-duration/conversion/log-interval scopes must not be ranked as speed benchmarks.",
            "OrthoHMM historical full-run provenance/resources remain separate from fresh reproduction and are not filled by proxy."])
    (output / "register.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    lines = ["# OrthoBench Score and Provenance Register", "", "Descriptive retained evidence; not a controlled runtime ranking.", "",
        "| Method | F1 (%) | Input sequence identity | Inference wall (s) | Timing scope |",
        "|---|---:|---|---:|---|"]
    for row in rows:
        seconds = "Unknown" if row["inference_wall_seconds"] is None else f'{row["inference_wall_seconds"]:,.3f}'
        lines.append(f'| {row["label"]} | {100*row["weighted_refog_f1"]:.4f} | {row["input_note"]} | {seconds} | {row["timing_basis"]} |')
    lines.extend(["", "OrthoFinder sequence-checkpoint conversion alone took 0.42 s; its inference time is unknown.",
        "FastOMA used a supplied tree and retains 90 successful/10 exit-137 tasks. The successful final collection is scored.",
        "OrthoMCL refers to July 2025 (8 BLAST threads), not April 2026 (32 threads).", "",
        "Prediction hashes, native coverage, explicit singleton padding, available commands/settings and resource semantics are in register.json.",
        "Proteinortho and July OrthoMCL inputs fail exact identity. Equal scores do not bound the effect of their preprocessing.",
        "No score replacement, missing-time imputation or complete-publication-provenance claim is made."])
    if include_orthohmm_audit:
        high = next(r for r in rows if r["key"] == "orthohmm_high_sensitivity")
        lines.extend(["", f"OrthoHMM high-sensitivity cached-hit replay: {high['replay_wall_seconds']:.6f} s; "
                      "excluded from the full-inference column because initial search was not timed here.",
                      "OrthoHMM phylogeny timing is the historical run, not the fresh installed reproduction."])
    (output / "register.md").write_text("\n".join(lines)+"\n")
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--include-orthohmm-audit", action="store_true")
    a = p.parse_args()
    assemble(a.repo.resolve(),a.output.resolve(),a.include_orthohmm_audit)
