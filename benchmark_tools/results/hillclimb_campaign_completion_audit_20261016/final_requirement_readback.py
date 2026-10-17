"""Final original-scope audit; a bounded plateau is not an accuracy victory."""

import collections
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys


def record(path):
    path = Path(path).resolve()
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1048576), b""):
            digest.update(block)
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": digest.hexdigest()}


def same(pin):
    actual = record(pin["path"])
    assert (actual["bytes"], actual["sha256"]) == (pin["bytes"], pin["sha256"]), pin["path"]
    return actual


def load(pin):
    same(pin)
    return json.loads(Path(pin["path"]).read_text())


def jobs_in(value, found):
    if isinstance(value, dict):
        for key, child in value.items():
            if (key in {"job_id", "producer_job_id", "original_job", "parent_job", "failed_job",
                        "successful_graph_recovery_job", "JobIDRaw"}
                    and type(child) in (str, int) and str(child).isdigit() and 24180 <= int(child) <= 24388):
                found.add(str(child))
            if key != "retained_details":
                jobs_in(child, found)
    elif isinstance(value, list):
        for child in value:
            jobs_in(child, found)
    elif isinstance(value, str) and "|" in value:
        states = {"COMPLETED", "FAILED", "CANCELLED", "TIMEOUT", "RUNNING", "PENDING", "OUT_OF_MEMORY"}
        for row in csv.reader(value.splitlines(), delimiter="|"):
            if len(row) >= 3 and row[0].isdigit() and 24180 <= int(row[0]) <= 24388 and row[1] in states:
                found.add(row[0])


def main(root, output):
    if output.exists():
        raise FileExistsError("Preserve previous final-audit attempts")
    ledger_path = Path("benchmark_tools/results/hillclimb_campaign_status_20261016.json")
    ledger_bytes = ledger_path.read_bytes()
    ledger = json.loads(ledger_bytes)
    objective_path = Path(ledger["active_objective"])
    objective_bytes = objective_path.read_bytes()
    assert all((str(i) + ". ").encode() in objective_bytes for i in range(1, 9))
    evidence = {key: record(value["path"]) for key, value in ledger.items()
                if isinstance(value, dict) and {"path", "bytes", "sha256"} <= value.keys()}
    for key, pin in evidence.items():
        same(ledger[key])
    data = {key: load(pin) for key, pin in evidence.items() if pin["path"].endswith(".json")}
    comparison = data["authoritative_comparison"]
    core = data["latest_core_protocol_readback"]
    supporting = data["latest_supporting_protocol_readback"]
    diagnostic = data["latest_endpoint_diagnostic_readback"]
    resources = data["latest_resource_interval_readback"]
    retained = data["latest_retained_comparator_identity_audit"]
    regression = data["latest_campaign_regression"]
    assert supporting["status"] == "supporting_protocol_evidence_review_verified" and supporting["unresolved"] == []
    assert supporting["protocol_count"] == 67 and supporting["group_count"] == 20
    assert core["status"] == "core_round_protocol_git_source_and_screening_readback_verified"
    assert len(core["rounds"]) == 10 and sum(len(r["plans"]) for r in core["rounds"]) == 16
    assert all(p["exact_plan_present_in_execution_revision"] for r in core["rounds"] for p in r["plans"])
    assert core["protected_union_strict_gate_passed"] is False
    methods = {m["key"]: m for m in comparison["methods"]}
    assert len(methods) == 28 and sum(m["complete_eight_endpoints"] for m in methods.values()) == 25
    assert {key: len(m["scores"]) for key, m in methods.items() if not m["complete_eight_endpoints"]} == {
        "deferred_fixed:control": 6, "deferred_fixed:deferred500": 6, "protected_union:protected_union": 1}
    cells = [cell for m in methods.values() for cell in m["cells"]]
    assert len(cells) == 91 and dict(collections.Counter(c["status"] for c in cells)) == ledger["counts"]["cell_status_counts"]
    assert diagnostic["current_rows"] == 224 and diagnostic["scored_current_rows"] == 213
    assert diagnostic["missing_current_endpoints"] == 11 and diagnostic["scored_current_F1_rows"] == 132
    assert resources["rows"] == 222 and resources["historical_cells"] == 91
    assert retained["method_configurations"] == 8 and retained["method_dataset_rows"] == 24
    assert retained["all_eight_method_records_and_24_cells_unchanged"] is True
    assert retained["all_retained_scores_versions_semantics_and_resources_reproduced"] is True
    assert retained["historical_input_consumption_proven"] is False
    assert regression["passed"] == 2349 and all(regression[k] == 0 for k in ("failures", "errors", "skipped"))
    assert supporting["primary_regression"]["executed_modules"] == 99
    assert supporting["reader_tests"]["counts"] == dict(tests=78, failures=0, errors=0, skipped=0)
    assert not ledger["defaults_promoted"] and not comparison["defaults_promoted"]
    selected = load(record("benchmarks/results/hillclimb_campaign_regression_v2_20261016/selection.json"))
    production = []
    for pin in selected["production_files"]:
        relative = Path(pin["path"]).relative_to(Path(selected["snapshot"])).as_posix() if not Path(pin["path"]).is_absolute() else None
        if relative is None:
            relative = Path(pin["path"]).relative_to(Path(selected["snapshot"]).resolve()).as_posix()
        current = record(relative)
        assert (current["bytes"], current["sha256"]) == (pin["bytes"], pin["sha256"])
        production.append(current)
    assert len(production) == 39
    libraries = []
    for pin in selected["libraries"]:
        current = record(Path("orthohmm/search/csrc") / Path(pin["path"]).name)
        assert (current["bytes"], current["sha256"]) == (pin["bytes"], pin["sha256"])
        libraries.append(current)
    for args in (("diff", "de89fd84..HEAD", "--", "orthohmm"), ("diff", "--", "orthohmm"),
                 ("diff", "--cached", "--", "orthohmm")):
        assert subprocess.check_output(["git", *args]) == b""
    samples = selected["primary_samples_before"]
    for field, command in (("status_sha256", ["git", "status", "--short", "--", "tests/samples"]),
            ("tracked_diff_sha256", ["git", "diff", "--binary", "--", "tests/samples"])):
        assert hashlib.sha256(subprocess.check_output(command)).hexdigest() == samples[field]
    assert subprocess.check_output(["git", "rev-list", "--left-right", "--count", "HEAD...@{upstream}"]).split() == [b"0", b"0"]
    jobs = set()
    for group in supporting["groups"]:
        for receipt in group["receipts"]:
            jobs_in(load(receipt["source"]), jobs)
    jobs.update(("24192", "24379", "24384", "24385", "24386", "24387", "24388"))
    accounting = subprocess.check_output(["sacct", "-j", ",".join(sorted(jobs)),
        "--format=JobIDRaw,State,ExitCode,AllocCPUS,Elapsed,NodeList", "--noheader", "--parsable2"], text=True)
    fields = ("JobIDRaw", "State", "ExitCode", "AllocCPUS", "Elapsed", "NodeList")
    accounting_rows = [dict(zip(fields, line.split("|"))) for line in accounting.splitlines()]
    parents = {row["JobIDRaw"]: row for row in accounting_rows if row["JobIDRaw"].isdigit()}
    assert jobs <= parents.keys(), sorted(jobs - parents.keys())
    terminal = {"COMPLETED", "FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY", "NODE_FAIL", "PREEMPTED", "BOOT_FAIL", "DEADLINE", "REVOKED"}
    assert all(parents[job]["State"].split()[0] in terminal for job in jobs)
    assert all(parents[job]["State"] == "COMPLETED" and parents[job]["ExitCode"] == "0:0"
               for job in ("24180", "24181", "24182", "24188", "24388"))
    assert parents["24387"]["State"] == "FAILED"
    scientific_winners = []
    reference = methods["orthofinder_3_1_5_full"]["scores"]
    for key, method in methods.items():
        if (key.startswith("orthohmm") or ":" in key) and method["complete_eight_endpoints"]:
            if all(method["scores"][endpoint] >= value for endpoint, value in reference.items()):
                scientific_winners.append(key)
    assert scientific_winners == []
    assert ledger["default_recommendation"] and len(ledger["next_evidence_driven_hypotheses"]) == 2
    requirements = [
        {"id": 1, "state": "fulfilled", "requirement": "Continue and close original evaluations without duplicate inference",
         "proof": "Original array24180/24181/24182 and supplement24188 are terminal successful.28 configurations/91 histories preserve failures, unavailable cells and separately admitted recoveries. Fixed-tree diagnostic gaps and union's stopped extension are explicit, not unrun required successful experiments.",
         "evidence": ["authoritative_comparison", "latest_core_protocol_readback", "latest_supporting_protocol_readback"]},
        {"id": 2, "state": "fulfilled_with_bounded_negative_findings", "requirement": "Diagnose fresh HMM recovery and broader duplication-aware family/retention behavior",
         "proof": "Six fresh search arms, both significance thresholds, cap400 and shared-hit complete vectors, protected-union fresh equivalence/pilot, native stage localization, held-tree substitution, copy-regime/component and satellite/tie mechanisms are retained. No broadly promising single mechanism warrants additional interaction sweeps; negative gates constrain further experiments rather than force expensive profile/tree changes.",
         "evidence": ["latest_stage_localization", "latest_core_protocol_readback", "latest_supporting_protocol_readback"]},
        {"id": 3, "state": "fulfilled", "requirement": "Bounded precommitted rounds, immutable sources and opt-in policies",
         "proof": "Sixteen core plan versions are present at executed revisions; all ten original screens and Pareto relations were independently reconstructed.67 protocols have explicit evidence and conditional-stop dispositions.16381 frozen source records and16325 Git code identities checked. No accuracy/default or private-kernel promotion. Adaptive denominator correctness repair is disclosed separately from unchanged fixed-cap controls.",
         "evidence": ["latest_core_protocol_readback", "latest_supporting_protocol_readback", "latest_campaign_regression"]},
        {"id": 4, "state": "fulfilled_with_nonblocking_historical_comparability_limits", "requirement": "Separate endpoints, full retained comparator panel, native validation and diagnostic/resource reporting",
         "proof": "Eight endpoint slots per configuration,213 scores/11NA;132 scored F1 entries have direct P/R.276 historical observations and222 scoped resource rows retain coverage/pairs/failures/CPU/wall/memory where available.24 retained tool/dataset rows and OrthoFinder3.1.5 full-versus-sequence semantics are checked. Contemporary trials bind identical inputs and scorers to matched controls; historical tool input consumption/transitive execution cannot universally be proven, remains disclosed, and limits strict all-tool matched-history claims.",
         "evidence": ["authoritative_score_table", "latest_endpoint_diagnostic_readback", "latest_resource_interval_readback", "latest_retained_comparator_identity_audit", "latest_orthofinder_output_semantics_audit"]},
        {"id": 5, "state": "fulfilled_with_nonblocking_confirmation_limits", "requirement": "Scientific validity, uncertainty, supporting TreeFam and feasible confirmation search",
         "proof": "All panels remain development-exposed. CGOB/ParaBench overlap and bounded public archive search are documented; no favorable remainder or unavailable confirmation is invented. Paired OrthoBench uncertainty preserves all family sufficient statistics and seeds. GO/EC/FAS remain conditional similarities, not F1 or proof of orthology. Unseeded FAS gains and unavailable original TreeFam family archives limit claims without blocking useful work.",
         "evidence": ["latest_uncertainty", "latest_supporting_protocol_readback", "latest_endpoint_diagnostic_readback"]},
        {"id": 6, "state": "fulfilled", "requirement": "Local bounded Threadripper execution and honest resource scopes",
         "proof": "Retained allocations, actual terminal accounting and competing-activity observations use local bizon without a DGX or quiet-window dependency. Native, cached, failed, conversion, scoring and recovery intervals remain separate. Process/parent/sampled-tree memory and missing values are not summed or imputed. No isolated speed ranking is claimed.",
         "evidence": ["latest_resource_interval_readback", "latest_supporting_protocol_readback"]},
        {"id": 7, "state": "fulfilled", "requirement": "Risk-scaled regression tests and scoped commits/pushes preserving user work",
         "proof": "2349 cases across the exact99 selected modules pass with zero errors/failures/skips; both scientific fixtures ran.78 current reporting tests pass, including the class/module correction.39 scientific production files and three original CPU libraries match the regression baseline.1174 user sample files and exact dirty-status/binary-diff fingerprints remain unchanged. Prior failed tests/readers are preserved; audit source is committed/pushed. Final result push must be verified before goal accounting is closed.",
         "evidence": ["latest_campaign_regression", "latest_supporting_protocol_readback"]},
        {"id": 8, "state": "fulfilled_as_documented_bounded_plateau", "requirement": "Reproducible ledger, separate-metric recommendation and honest conclusion",
         "proof": "Complete tables, ten-round negative/mixed decisions, explicit no-default-promotion recommendation and two evidence-driven future hypotheses are retained. No complete tested OrthoHMM configuration matches or exceeds OrthoFinder full on every endpoint. This is a bounded development plateau, not a global optimum, universal superiority, independent confirmation or a per-dataset winner composite.",
         "evidence": ["authoritative_comparison", "authoritative_score_table", "latest_core_protocol_readback", "latest_supporting_protocol_readback"]}
    ]
    assert [r["id"] for r in requirements] == list(range(1, 9))
    for requirement in requirements:
        assert all(key in evidence for key in requirement["evidence"])
    for pin in evidence.values():
        same(pin)
    assert ledger_path.read_bytes() == ledger_bytes and objective_path.read_bytes() == objective_bytes
    output.mkdir(parents=True)
    (output / "ledger_before_final_audit.json").write_bytes(ledger_bytes)
    (output / "objective.txt").write_bytes(objective_bytes)
    (output / "selected_job_accounting.tsv").write_text("\t".join(fields) + "\n" + "\n".join(
        "\t".join(row[key] for key in fields) for row in accounting_rows) + "\n")
    contrast_path = output / "canonical_vs_orthofinder_full.tsv"
    with contrast_path.open("x", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["dataset", "endpoint", "statistic", "orthohmm_canonical", "orthofinder_full", "score_difference", "precision_difference", "recall_difference"])
        for row in diagnostic["canonical_vs_orthofinder_full"]:
            writer.writerow([row["dataset"], row["endpoint"], row["statistic"], row["orthohmm"], row["orthofinder_full"],
                *["NA" if row["differences"][key] is None else row["differences"][key] for key in ("score", "precision", "recall")]])
    result = dict(status="original_eight_requirement_audit_verified_bounded_plateau_not_accuracy_victory", requirements=requirements,
        original_objective=record(output / "objective.txt"), ledger_snapshot=record(output / "ledger_before_final_audit.json"),
        evidence=evidence, recipe=record(__file__), production_files=production, original_cpu_libraries=libraries,
        selected_terminal_jobs=len(jobs), job_accounting=record(output / "selected_job_accounting.tsv"),
        per_metric_contrast=record(contrast_path), no_complete_orthohmm_full_vector_winner=True,
        default_recommendation=ledger["default_recommendation"], future_hypotheses=ledger["next_evidence_driven_hypotheses"],
        historical_comparability_limits=retained["limitations"], independent_confirmation_established=False,
        isolated_speed_ranking_established=False, unresolved_required_work=[],
        qualification="The original objective permits an evidence-supported bounded plateau and nonblocking scientific limitations. This audit does not claim the desired universal OrthoFinder victory or strict immutable historical input consumption for every retained tool.",
        supporting_reader_scope_correction="The supporting reader skipped retained_details recursion, but other explicit retained records still led it to hash the OrthoMCL BPO/index; the6122-file/65.345GiB totals and per-file dispositions are authoritative. Its limitation about inherited comparator records must not be read as claiming those files were not hashed at all.",
        scientific_inference_rerun=False, defaults_promoted=False, goal_completed=False,
        final_accounting_gate="Commit/push these audit results and final ledger, verify upstream synchronization and stable evidence, then close goal accounting.")
    with (output / "requirements.json").open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"requirements": len(requirements), "selected_terminal_jobs": len(jobs),
        "production_files": len(production), "result": str(output / "requirements.json")}))


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
