"""Derive the local 27-run storage plan without launching or preparing inference."""

import argparse
from copy import deepcopy
from pathlib import Path

from benchmark_tools.prepare_scaling_inputs import planned_runs
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.materialize_threadripper_inputs import ORDER_SHA
from benchmark_tools.snapshot_orthohmm_input_order import record
from benchmark_tools.probe_dgx_step_separation import save

PLAN_SHA = "25345e7a5d49e7474b09188498dc4760b298512266cd21510fba56ad6401e53d"


def derive(plan, order, output, inputs):
    output, inputs = Path(output), Path(inputs)
    for path in (output, inputs):
        if not path.is_absolute() or path.resolve() != path or path.exists() or path.is_symlink():
            raise ValueError("Require fresh absolute roots without symlink traversal")
    if output.is_relative_to("/dev/shm") or not inputs.is_relative_to("/dev/shm") or inputs == Path("/dev/shm"):
        raise ValueError("Require persistent output and private tmpfs input roots")
    keys = ("index", "method", "proteomes", "repeat")
    actual = [{key: row[key] for key in keys} for row in plan["runs"]]
    if actual != planned_runs():
        raise ValueError("Changed 27-run identities or order")
    if [r["proteomes"] for r in order["datasets"]] != [4, 8, 12]:
        raise ValueError("Changed enumeration dataset inventory")
    orders = {row["proteomes"]: row for row in order["datasets"]}
    rows = []
    for original in plan["runs"]:
        row = deepcopy(original)
        source = orders[row["proteomes"]]
        expected = {Path(r["path"]).name: (r["bytes"], r["sha256"]) for r in row["dataset"]["inputs"]}
        named = source["native_order"]
        if (len(expected) != row["proteomes"] or len(named) != len(expected) or set(named) != set(expected)
                or [Path(r["path"]).name for r in source["inputs_in_native_order"]] != named
                or any(expected[Path(r["path"]).name] != (r["bytes"], r["sha256"])
                       for r in source["inputs_in_native_order"])):
            raise ValueError("Retained input identities differ")
        old_root = Path(row["measurement_directory"]).parent
        new_root = output / f"run_{row['index']:02d}"
        fresh_input = inputs / f"run_{row['index']:02d}" / "input"
        config = row["configuration"]
        old_input = config.get("copy_inputs_to", row["dataset"]["input_directory"])

        def relocate(value):
            if value == old_input or value == row["dataset"]["input_directory"]:
                return str(fresh_input)
            path = Path(value)
            if path.is_absolute() and path.is_relative_to(old_root):
                return str(new_root / path.relative_to(old_root))
            return value

        row["native_argv"] = [relocate(v) for v in row["native_argv"]]
        config["argv"] = [relocate(v) for v in config["argv"]]
        for key in ("output", "metrics"):
            if key in config:
                config[key] = relocate(config[key])
        config.update(copy_inputs_from=row["dataset"]["input_directory"], copy_inputs_to=str(fresh_input))
        row.update(measurement_directory=str(new_root / "measurement"),
                   prepared_input_directory=str(fresh_input), input_creation_order=named,
                   original_native_argv=original["native_argv"])
        if row["native_method"] == "orthofinder_full":
            if "-o" in row["native_argv"] or "-op" in row["native_argv"]:
                raise ValueError("Unexpected output override or preparation-only flag")
            row["native_argv"] += ["-o", config["output"]]
            config["argv"] = list(row["native_argv"])
            row["expected_native_order"] = sorted(named)
        else:
            row["expected_native_order"] = named
        rows.append(row)
    return dict(status="threadripper_commands_prepared_unrun", runs=rows,
        source_plan_sha256=PLAN_SHA, source_order_sha256=ORDER_SHA,
        environment_overrides=plan["environment_overrides"], prepend_path=plan["prepend_path"],
        resolved_executables=plan["resolved_executables"], scientific_execution_authorized=False,
        resources=dict(host="bizon", native_workers=32, native_affinity=list(range(32)),
                       scheduler_cpu_slots=64, memory_bytes=128*1024**3, native_timeout_s=85800),
        limitations=["Plan only: inputs and output directories are not created; nothing submitted.",
                     "Requires new per-run preparation for all methods; historical preparer is not compatible.",
                     "Runtime/input rechecks, full-pipeline/overhead validation and quiet-host gates remain mandatory."])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--input-root", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    args = parser.parse_args()
    base = Path(__file__).parent / "results"
    plan = base / "publication_scaling_commands_20260917.json"
    order = base / "dgx_native_input_order_20260917.json"
    result = derive(read_frozen(plan, PLAN_SHA), read_frozen(order, ORDER_SHA), args.output_root, args.input_root)
    result["sources"] = [record(p) for p in (plan, order, Path(__file__))]
    save(args.manifest, result)


if __name__ == "__main__":
    main()
