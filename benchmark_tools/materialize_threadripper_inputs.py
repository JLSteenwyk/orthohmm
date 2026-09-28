"""Copy frozen inputs in retained order, then test the actual native enumerator."""

import argparse
import json
from pathlib import Path
import shutil

from benchmark_tools.snapshot_orthohmm_input_order import snapshot, record
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.run_simulation_tree_mode_control import METHOD_SHA
from benchmark_tools.prepare_scaling_commands import INPUT_SHA
from benchmark_tools.probe_dgx_step_separation import save

ORDER_SHA = "7972eb5e1224c0f14dc09a36b26d38f777fddfad3b6075cb5d20999a55a5a451"


def materialize(inputs, retained, baseline, output, reverse_creation=False):
    output = Path(output).absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    for manifest in (inputs, retained):
        if [d["proteomes"] for d in manifest["datasets"]] != [4, 8, 12]:
            raise ValueError("Require complete ordered 4/8/12 inventory")
    plans = []
    for dataset, reference in zip(inputs["datasets"], retained["datasets"]):
        files = {Path(row["path"]).name: row for row in dataset["inputs"]}
        names = reference["native_order"]
        if (len(files) != dataset["proteomes"] or len(files) != len(dataset["inputs"])
                or len(names) != len(files) or set(names) != set(files)):
            raise ValueError("Membership or unique basename inventory differs")
        expected = reference["inputs_in_native_order"]
        if [Path(r["path"]).name for r in expected] != names:
            raise ValueError("Retained ordered hashes differ from names")
        for name, row in zip(names, expected):
            if Path(name).name != name:
                raise ValueError("Require basenames")
            actual = record(Path(files[name]["path"]))
            if actual != files[name] or any(actual[k] != row[k] for k in ("bytes", "sha256")):
                raise ValueError("Input bytes differ")
        plans.append((dataset, names, files))
    output.mkdir(parents=True)
    copied = []
    for dataset, names, files in plans:
        directory = output / str(dataset["proteomes"])
        directory.mkdir()
        rows = []
        for name in (list(reversed(names)) if reverse_creation else names):
            target = directory / name
            shutil.copy2(files[name]["path"], target)
            row = record(target)
            if any(row[k] != files[name][k] for k in ("bytes", "sha256")):
                raise ValueError("Copy differs from input")
            rows.append(row)
        copied.append(dict(proteomes=dataset["proteomes"], input_directory=str(directory), inputs=rows))
    runtime = {"records": [dict(path=r["absolute_path"], bytes=r["bytes"], sha256=r["sha256"])
                           for r in baseline["core_sources"]]}
    native = snapshot(baseline["core_root"], {"datasets": copied}, runtime)
    matches = [a["native_order"] == b["native_order"]
               for a, b in zip(native["datasets"], retained["datasets"])]
    result = dict(status="native_order_matched" if all(matches) else "native_order_mismatched",
                  datasets=copied, native_snapshot=native, matches_retained=matches,
                  source=record(__file__), reverse_creation=reverse_creation,
                  scientific_execution_authorized=False,
                  limitations=["Creation order is not assumed to guarantee enumeration; actual enumerator is tested.",
                               "Recheck immediately before and after each inference; fresh OrthoFinder copies need validation."])
    save(output / "manifest.json", result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--reverse-creation", action="store_true")
    args = parser.parse_args()
    results = Path(__file__).resolve().parent / "results"
    paths = [(results / "publication_scaling_inputs_20260916.json", INPUT_SHA),
             (results / "dgx_native_input_order_20260917.json", ORDER_SHA),
             (results / "publication_variable_native_methods_20260916.json", METHOD_SHA)]
    manifests = [read_frozen(path, digest) for path, digest in paths]
    result = materialize(*manifests, args.output, reverse_creation=args.reverse_creation)
    save(args.output / "provenance.json", {"sources": [record(path) for path, _ in paths]})
    print(result["status"])
    if result["status"] != "native_order_matched":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
