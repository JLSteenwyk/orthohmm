"""Measure raw JSON retention only, not native workload or collector overhead."""

import argparse
import json
from pathlib import Path
import sys
import time
import tracemalloc

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.disk_observation_sequence import DiskObservations
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.probe_dgx_step_separation import save


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--point", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--count", type=int, default=1000)
    parser.add_argument("--mode", choices=("list", "disk"), required=True)
    args = parser.parse_args()
    if not 1 <= args.count <= 10000:
        raise ValueError("Require bounded diagnostic count")
    args.output.mkdir(parents=True, exist_ok=False)
    raw = args.point.read_text()
    points = [] if args.mode == "list" else DiskObservations(args.output)
    tracemalloc.start()
    started = time.monotonic()
    for _ in range(args.count):
        points.append(json.loads(raw))
    wall = time.monotonic() - started
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    if points[0] != json.loads(raw) or points[-1] != json.loads(raw):
        raise ValueError("Retained endpoints differ")
    save(args.output / "result.json", dict(mode=args.mode, count=len(points),
        wall_s=wall, traced_current_bytes=current, traced_peak_bytes=peak,
        raw=record(args.point), source=record(__file__),
        sequence_source=record(Path(__file__).with_name("disk_observation_sequence.py")),
        scientific_timings_admitted=False,
        limitations=["Repeated decoding of one retained snapshot, not live collection.",
                     "tracemalloc measures Python allocations, not process RSS, page cache or total job memory.",
                     "Disk mode writes every raw point; list mode measures the old retained history, not its additional writes.",
                     "No claim about native slowdown, full-day peak memory, thread-count scaling or post-run interval reports."]))


if __name__ == "__main__":
    main()
