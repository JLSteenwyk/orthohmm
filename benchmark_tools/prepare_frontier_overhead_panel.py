"""Derive a prospective paired native overhead panel from frozen commands."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.launch_dgx_native_run import read_pinned
from benchmark_tools.prepare_ob_candidate_neighborhood import record

PLAN_SHA = "7096348236f8372ef7f6ad12a3e829eae3dfee812b17c3b0bdd582480538ff1b"
SPEC_SHA = "fe96893489935124b675e3fb06317870cc7c155d9765e08067ad5b38b83df362"
ROOT = "/home/jlsteenwyk/projects/orthohmm-publication"
METHODS = ("orthohmm_high_sensitivity", "orthohmm_satellite_v2", "orthofinder_full")


def relocate(value, old, new):
    if isinstance(value, str):
        return new + value[len(old):] if value == old or value.startswith(old + "/") else value
    if isinstance(value, dict):
        return {key: relocate(item, old, new) for key, item in value.items()}
    if isinstance(value, list):
        return [relocate(item, old, new) for item in value]
    return value


def build(plan_path, spec_path):
    plan = read_pinned(plan_path, PLAN_SHA)
    spec = read_pinned(spec_path, SPEC_SHA)
    if spec["runs"] != plan["runs"] or spec["original_plan"]["sha256"] != PLAN_SHA:
        raise ValueError("Execution commands differ from pinned plan")
    templates = plan["runs"][:3]
    if tuple(r["native_method"] for r in templates) != METHODS:
        raise ValueError("Unexpected native method order")
    if any(r["dataset"]["proteins"] != 73266 or r["proteomes"] != 4 or r["repeat"] != 0 for r in templates):
        raise ValueError("Require original four-proteome templates")
    rows = []
    # Rotate methods between paired blocks; reverse arm order in the middle block.
    for pair in range(3):
        for offset in range(3):
            method_index = (pair + offset) % 3
            template = templates[method_index]
            modes = ("boundary", "periodic") if pair != 1 else ("periodic", "boundary")
            for mode in modes:
                index = len(rows)
                old = f"{ROOT}/scaling_native_v1/run_{method_index:02d}"
                new = f"{ROOT}/frontier_overhead_v1/run_{index:02d}"
                run = relocate(template, old, new)
                run.update(index=index, repeat=pair)
                rows.append(dict(index=index, pair=pair, mode=mode, method=template["native_method"],
                                 original_index=method_index, run=run))
    result = {key: spec[key] for key in ("runtime_manifests", "enumerator", "environment_paths",
                                        "environment_overrides", "unset_environment")}
    result.update(schema=1, status="prospective_native_frontier_overhead_plan",
        source=record(__file__), original_plan=record(plan_path), original_execution=record(spec_path),
        core_commit=plan["core_commit"], runs=rows,
        order=next(order for order in spec["orders"] if order["proteomes"] == 4),
        resource_plan=dict(host="spark-7ff0", partition="spark", cpus=20, memory_gib=96,
                           exclusive=True, concurrency=1, gpu=False, scheduler_limit_s=3600),
        native_timeout_s=900, interval_s=1., eligibility_delay_s=60,
        overhead_statistic="periodic native command wall / paired boundary native command wall - 1",
        engineering_budget=dict(per_method_median_max=.05, every_pair_max=.10),
        execution_authorized=False, scientific_timings_admitted=False, publication_ready=False,
        limitations=["Three paired replicates per method; engineering diagnostic, not a significance test.",
            "Representative smallest scaling input only; no claim about eight/twelve-proteome overhead.",
            "Boundary-only is not uninstrumented and cannot establish interval-level quietness.",
            "Report all failures, missing pairs and screen flags; no selective repeats or overhead subtraction.",
            "Budget results and environmental/screen validity must be reported separately.",
            "Requires a verified launcher/recipe and frozen protocol before submission."])
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--spec", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = build(args.plan.resolve(), args.spec.resolve())
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
