"""Verify adaptive-cap arithmetic with explicit stubbed Phase A outcomes."""

import argparse
import json
from pathlib import Path
import sys
from types import SimpleNamespace
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_membership_hillclimb import write_json
from orthohmm.search import engine, prefilter
from orthohmm.search.sequences import SpeciesSequences

CASES = ((500,400), (512,512), (750,600), (1024,410))


def probe_case(query_count, successes):
    if type(query_count) is not int or query_count <= 0 or type(successes) is not int or successes < 0:
        raise ValueError("Require positive query count and nonnegative explicit successes")
    nominal = min(query_count,500)
    step = max(1,query_count//nominal)
    expected_probes = len(range(0,query_count,step))
    if successes > expected_probes:
        raise ValueError("Stub successes exceed actual probe count")
    query = SpeciesSequences("synthetic-query",[f"q{i}" for i in range(query_count)],
        np.zeros(4*query_count,dtype=np.uint8),4*np.arange(query_count,dtype=np.int64),
        np.full(query_count,4,dtype=np.int32))
    target = SpeciesSequences("synthetic-target",["target"],np.zeros(4,dtype=np.uint8),
        np.zeros(1,dtype=np.int64),np.full(1,4,dtype=np.int32))
    index = SimpleNamespace(alpha_size=20, kmer_offsets=np.zeros(1,dtype=np.int32),
        kmer_entries=np.zeros(1,dtype=np.int32),kmer_freqs=np.zeros(1,dtype=np.int32),freq_threshold=1)
    calls, observed_hits = 0,0
    def explicit_phase_a(*args,**kwargs):
        nonlocal calls,observed_hits
        calls += 1
        success = calls <= successes
        observed_hits += int(success)
        return (np.zeros(1 if success else 0,dtype=np.int32),np.zeros(1 if success else 0,dtype=np.float64))
    with patch.object(prefilter,"_score_query_two_stage",explicit_phase_a):
        cap, coverage = engine._adaptive_max_candidates(query,target,4,4,False,prepared_index=index)
    if calls != expected_probes or observed_hits != successes:
        raise ValueError("Actual adaptive probe loop differs from prespecified source")
    observed_coverage = successes/calls
    observed_cap = 20 if observed_coverage > .8 else 60 if observed_coverage > .5 else 100 if observed_coverage > .2 else 200
    return {"queries":query_count,"explicit_stub_successes":successes,"nominal_sample_size":nominal,
        "actual_probe_calls":calls,"stride":step,"legacy_reported_coverage":coverage,
        "observed_fraction":observed_coverage,"legacy_cap":cap,"observed_fraction_cap":observed_cap,
        "reported_coverage_exceeds_one":coverage>1,"cap_would_change":cap!=observed_cap}


def run(output):
    if output.exists():
        raise ValueError("Preserve existing adaptive denominator audit")
    source = record(engine.__file__)
    plan_path = Path(__file__).with_name("hillclimb_adaptive_coverage_denominator_audit_20261014.json")
    plan = json.loads(plan_path.read_text())
    if plan["status"] != "prespecified_before_reproducer_execution":
        raise ValueError("Unexpected adaptive arithmetic audit plan")
    rows = [probe_case(*case) for case in CASES]
    if rows[2]["legacy_cap"] != 20 or rows[2]["observed_fraction_cap"] != 60 or not rows[1]["reported_coverage_exceeds_one"]:
        raise ValueError("Prespecified unchanged-source reproducer no longer demonstrates the issue")
    check(source)
    write_json(output,{"status":"adaptive_probe_denominator_inflation_reproduced_without_source_changes",
        "source":source,"reader":record(__file__),"plan":record(plan_path),"rows":rows,
        "scope":"Explicit synthetic/stubbed branch arithmetic only; no biological match outcomes, HMM scoring, benchmark rescoring or source/default modification.",
        "fixed_cap_panels_affected":False,"defaults_promoted":False,"goal_completed":False})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    args = parser.parse_args()
    run(args.output.resolve())
