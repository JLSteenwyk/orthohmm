"""Freeze scientific DGX execution only after the pinned engineering gates pass."""

import argparse
from copy import deepcopy
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.launch_dgx_native_run import read_pinned

PROJECT = Path("/home/jlsteenwyk/projects/orthohmm-publication")


def authorize(results):
    plan = read_pinned(results / "dgx_scaling_commands_20260917.json",
                       "7096348236f8372ef7f6ad12a3e829eae3dfee812b17c3b0bdd582480538ff1b")
    smoke = read_pinned(results / "dgx_launcher_smoke_spec_20260917.json",
                        "8af8a049c69b7e479e8443e1e0e14b74c78d67b345c03d7f557768488d057fa2")
    admission = read_pinned(results / "dgx_native_launcher_smokes_admitted_20260917.json",
                            "85f05ed52652076d904fc0870fbc087c60d22be58f1a6316226d1fe72c1a393e")
    overhead = read_pinned(results / "dgx_verified_overhead_20260917.json",
                           "0184d54be2f55720aff35f27b6db0123de28ca7dfcdabd50cbf1d8d79d60b437")
    order = read_pinned(results / "dgx_native_input_order_20260917.json",
                        "7972eb5e1224c0f14dc09a36b26d38f777fddfad3b6075cb5d20999a55a5a451")
    if admission["status"] != "all_native_launcher_smokes_admitted" or not overhead["all_protocol_gates_met"]:
        raise ValueError("Engineering gates have not passed")
    for field in ("environment_overrides", "environment_paths", "unset_environment", "resource_plan"):
        if smoke[field] != plan[field]:
            raise ValueError("Scientific environment differs from tested smoke")
    for run in plan["runs"]:
        matching = [d for d in order["datasets"] if d["input_directory"] == run["dataset"]["input_directory"]]
        if len(matching) != 1 or sorted(matching[0]["inputs_in_native_order"], key=lambda r: r["path"]) != sorted(run["dataset"]["inputs"], key=lambda r: r["path"]):
            raise ValueError("Scientific order/input identities differ")
    spec = deepcopy(smoke)
    recipe = PROJECT / "native_launcher_recipe_v2/benchmark_tools"
    spec.update(purpose="scientific_scaling", execution_authorized=True,
                scientific_timing_runs_authorized=True, runs=plan["runs"], orders=order["datasets"],
                enumerator=order["source"], native_timeout_s=85800,
                validated_launcher_smoke={"path": str(recipe / "dgx_native_launcher_smokes_admitted_20260917.json"),
                                          "sha256": "85f05ed52652076d904fc0870fbc087c60d22be58f1a6316226d1fe72c1a393e"},
                original_plan={"path": str(recipe / "dgx_scaling_commands_20260917.json"),
                               "sha256": "7096348236f8372ef7f6ad12a3e829eae3dfee812b17c3b0bdd582480538ff1b"},
                overhead={"path": str(recipe / "dgx_verified_overhead_20260917.json"),
                          "sha256": "0184d54be2f55720aff35f27b6db0123de28ca7dfcdabd50cbf1d8d79d60b437"},
                limitations=["27 fresh sequential matched-resource runs; authorization is not completed evidence.",
                             "Native timeout23h50m leaves10minutes for preparation/checks within each24h allocation.",
                             "Hashing/copying warm caches; these are not cold-cache timings.",
                             "Preserve every failure; independently admit native outputs and resource evidence after execution.",
                             "No pooling with x86 timings, overhead correction or outcome-based default changes."])
    return spec


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = authorize(args.results)
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
