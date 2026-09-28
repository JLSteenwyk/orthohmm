"""Stress retained report structures; synthetic slots are never native evidence."""

import argparse
from copy import deepcopy
import json
from pathlib import Path
import resource
import time

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_dgx_step_separation import save


def repeat(values, count):
    if not isinstance(values, list) or not values:
        raise ValueError("Require nonempty retained template list")
    return [deepcopy(values[i % len(values)]) for i in range(count)]


def expand(lineage, context, count):
    if type(count) is not int or not 2 <= count <= 85802:
        raise ValueError("Require 2-85802 synthetic observation slots")
    if lineage.get("schema") != "threadripper_scaling_v3":
        raise ValueError("Require retained Threadripper v3 report")
    a, b = deepcopy(lineage), deepcopy(context)
    a["point_records"] = repeat(a["point_records"], count)
    s = a["screening"]
    original = s["original_screening"]
    threshold = original["original_threshold_screen"]
    for owner, key in ((s, "intervals"), (original, "hierarchy_intervals"),
                       (threshold, "intervals"), (b["context"], "intervals")):
        owner[key] = repeat(owner[key], count - 1)
    # Include a diagnostic flag for every slot to avoid assuming a quiet run.
    s["narrow_flagged_intervals"] = list(range(count - 1))
    threshold["flagged_intervals"] = list(range(count - 1))
    b["context"]["observations"] = count
    for value in (a, b):
        value.update(status="synthetic_report_size_probe", synthetic=True,
                     scientific_timings_admitted=False, publication_ready=False)
        value["limitations"] = ["Repeated report entries, not observations or valid replay evidence."]
    return a, b


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--measurement", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--count", type=int, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    paths = [args.measurement / name for name in ("lineage_report.json", "root_context_report.json")]
    inputs = [record(p) for p in paths]
    started = time.monotonic()
    a, b = expand(*(json.loads(p.read_text()) for p in paths), args.count)
    construction = time.monotonic() - started
    constructed_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    started = time.monotonic()
    outputs = []
    for name, value in (("synthetic_lineage.json", a), ("synthetic_context.json", b)):
        path = args.output / name
        save(path, value)
        outputs.append(record(path))
    serialization = time.monotonic() - started
    for ref in inputs:
        check(ref)
    save(args.output / "result.json", dict(status="synthetic_report_storage_measured",
        observation_slots=args.count, inputs=inputs, outputs=outputs,
        construction_s=construction, serialization_and_hashing_s=serialization,
        constructed_max_rss_kib=constructed_rss,
        final_max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        source=record(__file__), writer_source=record(Path(__file__).with_name("probe_dgx_step_separation.py")),
        scientific_timings_admitted=False, publication_ready=False,
        limitations=["Linux process high-water RSS, not cgroup memory or page cache.",
            "Repeated fixture report rows; no chronological observations or valid timing replay.",
            "Measures retained report construction and serialization, not evaluator temporaries, live sampling or native slowdown.",
            "Full-day slot count is not a bound on reports with larger process/step inventories."]))


if __name__ == "__main__":
    main()
