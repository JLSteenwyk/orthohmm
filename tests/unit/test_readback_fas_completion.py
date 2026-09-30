import copy

import pytest

from benchmark_tools.readback_fas_completion import PRIOR_SHA, QUERY, check, pin, readback, terminal_fields


@pytest.fixture
def evidence():
    methods = [dict(key=f"method_{i}", details={"FAS": {"assessed_relations": 10}}) for i in range(8)]
    natives = [dict(method=m["key"], strata={"precomputed": {"population": 6},
               "missing": {"population": 4}}, unannotated_pairs_logged=0) for m in methods]
    rows = [dict(method=m["key"], native_counts_match=True,
        saved_lookup_strata_and_values_match=True, database_historically_hash_bound=False,
        reused_prior_recount=i < 6, distinct_query_pairs=12, skipped_alias_pairs=2,
        precomputed=6, missing=4, unannotated=0, eligible_pairs=10,
        native_logged_counts={"precomputed": 6, "missing": 4, "unannotated": 0},
        precomputed_score_sum=3.0, precomputed_mean=0.5,
        hypothetical_full_mean_bounds=[0.3, 0.7]) for i, m in enumerate(methods)]
    sources = [dict(partial_row={"path": str(i)}) for i in range(6)]
    partials = [{k: v for k, v in r.items() if k != "reused_prior_recount"} for r in rows[:6]]
    for row, source in zip(rows, sources):
        row["prior_row_source"] = source["partial_row"]
    flags = dict(uncertainty_admitted=False, benchmark_scores_changed=False, publication_ready=False)
    prior = dict(status="timed_out_incomplete_population_recount", job_id=22382,
                 lookup={"entries": 20}, completed_partial_rows=sources, **flags)
    report = dict(status="retained_fas_eligible_populations_completed_with_reuse", methods=rows,
        query=QUERY, lookup=prior["lookup"], **flags,
        reuse=dict(prior_receipt={"sha256": PRIOR_SHA}, original_job_id=22382,
                   original_state="TIMEOUT", original_final_stability_pass_completed=False,
                   historical_parser_hash_identity_established=False,
                   completion_input_stability_pass_completed=True,
                   fresh_database_recounts=2, reused_database_recounts=6))
    return report, {"methods": methods}, {"methods": natives}, prior, partials


def test_complete_readback(evidence):
    rows = readback(*evidence)
    assert len(rows) == 8
    assert [r["reused_prior_recount"] for r in rows] == [True] * 6 + [False] * 2
    assert rows[-1]["missing_fraction"] == 0.4


@pytest.mark.parametrize("mutation", ["panel", "bounds", "mean", "native", "reuse",
                                     "historical", "query", "uncertainty", "count", "prior"])
def test_changed_evidence_rejected(evidence, mutation):
    report, manifest, attrition, prior, partials = copy.deepcopy(evidence)
    row = report["methods"][-1]
    if mutation == "panel":
        report["methods"].pop()
    elif mutation == "bounds":
        row["hypothetical_full_mean_bounds"][1] = 0.8
    elif mutation == "mean":
        row["precomputed_mean"] = 0.6
    elif mutation == "native":
        attrition["methods"][-1]["strata"]["missing"]["population"] = 3
    elif mutation == "reuse":
        row["reused_prior_recount"] = True
    elif mutation == "historical":
        report["reuse"]["historical_parser_hash_identity_established"] = True
    elif mutation == "query":
        report["query"] += " LIMIT 10"
    elif mutation == "uncertainty":
        report["uncertainty_admitted"] = True
    elif mutation == "count":
        row["unannotated"] = False
    else:
        partials[0]["precomputed_score_sum"] = 2.5
    with pytest.raises(ValueError):
        readback(report, manifest, attrition, prior, partials)


def test_pin_detects_same_size_change(tmp_path):
    path = tmp_path / "input.json"
    path.write_text("{}")
    ref = pin(path)
    check(ref)
    path.write_text("[]")
    with pytest.raises(ValueError):
        check(ref)


TERMINAL = ("JobId=22383 JobState=COMPLETED ExitCode=0:0 Restarts=0 Requeue=0 "
            "NodeList=bizon NumCPUs=1 NumTasks=1 MinMemoryNode=16G TimeLimit=03:00:00 \n")


def test_successful_original_terminal():
    assert terminal_fields(TERMINAL)["JobState"] == "COMPLETED"


@pytest.mark.parametrize("raw", [TERMINAL.replace("COMPLETED", "RUNNING"),
    TERMINAL.replace("0:0", "0:15"), TERMINAL.replace("Restarts=0", "Restarts=1"),
    TERMINAL.replace("22383", "22382"), TERMINAL.replace("NumCPUs=1", "NumCPUs=2"),
    TERMINAL.replace("03:00:00", "04:00:00"), TERMINAL + "JobState=COMPLETED\n",
    TERMINAL.strip() + " JobState=COMPLETED\n", TERMINAL.strip() + " ArrayJobId=22383\n"])
def test_terminal_scope_rejected(raw):
    with pytest.raises(ValueError):
        terminal_fields(raw)
