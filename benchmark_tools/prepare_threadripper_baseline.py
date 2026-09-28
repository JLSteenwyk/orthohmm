"""Amend cwd-dependent package metadata without changing native runtime files."""

import argparse
from copy import deepcopy
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.run_simulation_methods import execution_environment, read_frozen
from benchmark_tools.snapshot_orthohmm_input_order import record
from benchmark_tools.probe_dgx_step_separation import save

QUERY = "import importlib.metadata as m,json,sys; names=sorted({d.metadata['Name'] for d in m.distributions() if d.metadata['Name']}); print(json.dumps({'python':sys.version,'packages':{n:m.version(n) for n in names}}))"


def amend(baseline, observed):
    expected = deepcopy(baseline["environments"])
    if expected["orthofinder"]["packages"].pop("orthohmm") != "0.5.0":
        raise ValueError("Unexpected historical repository metadata version")
    if observed != expected:
        raise ValueError("Environment differs beyond known cwd-only metadata")
    result = deepcopy(baseline)
    result["environments"] = observed
    result["threadripper_amendment"] = dict(
        reason="Historical OrthoFinder inventory included repository-local orthohmm.egg-info via cwd.",
        execution_cwd=baseline["core_root"],
        difference={"orthofinder": {"orthohmm": {"historical": "0.5.0", "execution": None}}},
        runtime_files_changed=False, scientific_execution_authorized=False)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    baseline = read_frozen(args.baseline, args.sha256)
    env, _ = execution_environment(baseline)
    for key in ("LD_LIBRARY_PATH", "LD_PRELOAD", "LD_AUDIT"):
        env.pop(key, None)
    interpreters = dict(orthohmm=baseline["tool_entrypoints"]["orthohmm_python"]["absolute_path"],
        orthofinder=str(Path(baseline["tool_entrypoints"]["orthofinder"]["absolute_path"]).parent / "python"))
    observed = {name: json.loads(subprocess.check_output([exe, "-B", "-c", QUERY],
        env=env, cwd=baseline["core_root"], text=True)) for name, exe in interpreters.items()}
    result = amend(baseline, observed)
    result["threadripper_amendment"].update(source=record(__file__), baseline=record(args.baseline))
    save(args.output, result)


if __name__ == "__main__":
    main()
