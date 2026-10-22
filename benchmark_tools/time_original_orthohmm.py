"""Observe stage boundaries in an unchanged original OrthoHMM source tree."""

import argparse
from functools import wraps
import importlib
import json
from pathlib import Path
import sys
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--metrics", type=Path, required=True)
    parser.add_argument("arguments", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    source = args.source.resolve()
    sys.path.insert(0, str(source))
    module = importlib.import_module("orthohmm.orthohmm")
    if not Path(module.__file__).resolve().is_relative_to(source):
        raise RuntimeError("Original implementation imported from the wrong source")
    stages = {}
    names = dict(execute_phmmer_search="search", determine_edge_thresholds="edge_thresholds",
                 determine_network_edges="network_edges", execute_mcl="clustering",
                 generate_orthogroup_clusters_file="orthogroup_materialization",
                 generate_orthogroup_files="output")

    def observed(function, name):
        @wraps(function)
        def call(*values, **kwargs):
            started = time.perf_counter()
            try:
                return function(*values, **kwargs)
            finally:
                stages[name] = dict(wall_s=time.perf_counter() - started)
        return call

    for function, name in names.items():
        setattr(module, function, observed(getattr(module, function), name))
    arguments = args.arguments[1:] if args.arguments[:1] == ["--"] else args.arguments
    sys.argv = ["original-orthohmm", *arguments]
    status = "failed"
    try:
        module.main()
        status = "complete"
    except SystemExit as error:
        if error.code in (None, 0):
            status = "complete"
        raise
    finally:
        args.metrics.write_text(json.dumps(dict(status=status, stages=stages,
            source=str(source), module=str(Path(module.__file__).resolve()),
            scope="Observer wraps top-level calls only; original source and search commands are unchanged."),
            indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
