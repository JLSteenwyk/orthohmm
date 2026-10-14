"""Invoke the unchanged production CLI with observable, fail-closed CPU routes."""

import argparse
import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import numpy as np
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_membership_hillclimb import write_json


def install_guard(engine, profiles, viterbi, directory):
    library, available = viterbi._load_c_hmm()
    if (not available or library is None
            or not hasattr(library, "batch_hmm_viterbi_multipair_avx2_c")
            or engine.is_cuda_available() or profiles.is_cuda_available()):
        raise ValueError("Require available CPU multipair route and hidden CUDA")
    identity = record(Path(library._name))
    directory.mkdir(parents=True, exist_ok=False)
    counters = {}

    def forbidden(*args, **kwargs):
        raise AssertionError("Full native validation forbids silent non-multipair fallback")

    def observed(original, route):
        def call(*args, **kwargs):
            scores = original(*args, **kwargs)
            if scores.dtype != np.int32 or scores.shape != (len(args[8]),):
                raise ValueError("Native kernel returned an invalid score array")
            pid = os.getpid()
            state = counters.setdefault(pid, {"pid": pid, "library": identity, "routes": {}})
            stats = state["routes"].setdefault(route, {"calls": 0, "scored_pairs": 0})
            stats["calls"] += 1
            stats["scored_pairs"] += len(scores)
            write_json(directory / f"routes_{pid}.json", state)
            return scores
        return call

    for module, route in ((engine, "initial_search"), (profiles, "profile_scoring")):
        module.batch_viterbi_multipair_c = observed(module.batch_viterbi_multipair_c, route)
        for name in ("batch_viterbi_c", "batch_viterbi_score", "batch_viterbi_cuda"):
            setattr(module, name, forbidden)
    return identity


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, add_help=False)
    parser.add_argument("--trace-directory", required=True, type=Path)
    args, production_args = parser.parse_known_args(argv)
    if os.environ.get("CUDA_VISIBLE_DEVICES") != "":
        raise ValueError("Hide CUDA before importing the native engine")
    from orthohmm.search import engine, profile_expansion, viterbi
    from orthohmm.orthohmm import main as production_main
    source = Path(__file__).resolve().parent.parent
    if Path(viterbi.__file__).resolve() != source / "orthohmm/search/viterbi.py":
        raise ValueError("Wrong frozen native implementation")
    library = install_guard(engine, profile_expansion, viterbi, args.trace_directory.resolve())
    # Production main currently parses sys.argv, even when its argv argument is supplied.
    previous_argv = sys.argv
    try:
        sys.argv = ["orthohmm", *production_args]
        production_main(production_args)
    finally:
        sys.argv = previous_argv
    check(library)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
