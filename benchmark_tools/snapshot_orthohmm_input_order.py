"""Call the frozen OrthoHMM enumerator and verify exact scaling input bytes."""

import argparse
import hashlib
import json
from pathlib import Path
import sys
import types


def record(path):
    path = Path(path).resolve()
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": digest.hexdigest()}


def snapshot(core, inputs, runtime):
    source = Path(core).resolve() / "orthohmm/files.py"
    source_record = record(source)
    matches = [row for row in runtime["records"] if row["path"] == str(source)]
    if len(matches) != 1 or matches[0].get("sha256") != source_record["sha256"]:
        raise ValueError("Enumerator source is not the frozen runtime source")
    # Execute the verified source directly, bypassing its excluded bytecode cache.
    module = types.ModuleType("frozen_orthohmm_files")
    module.__file__ = str(source)
    exec(compile(source.read_bytes(), str(source), "exec"), module.__dict__)
    datasets = []
    for dataset in inputs["datasets"]:
        directory = Path(dataset["input_directory"]).resolve()
        expected = {Path(row["path"]).name: row for row in dataset["inputs"]}
        if len(expected) != dataset["proteomes"] or len(expected) != len(dataset["inputs"]):
            raise ValueError("Duplicate or incorrect input inventory")
        names = module.fetch_fasta_files(str(directory))
        if len(names) != len(expected) or set(names) != set(expected):
            raise ValueError("Native enumeration differs from frozen input membership")
        observed = []
        for name in names:
            if Path(name).name != name:
                raise ValueError("Native enumeration returned a non-basename")
            row = record(directory / name)
            if row != expected[name]:
                raise ValueError("Input bytes or path differ from frozen manifest")
            observed.append(row)
        if names != module.fetch_fasta_files(str(directory)):
            raise ValueError("Native enumeration changed during snapshot")
        datasets.append({"proteomes": dataset["proteomes"], "input_directory": str(directory),
                         "native_order": names, "inputs_in_native_order": observed})
    if source_record != record(source):
        raise ValueError("Enumerator source changed during snapshot")
    return {"status": "orthohmm_native_enumeration_captured", "source": source_record,
            "interpreter": record(sys.executable), "python_version": sys.version,
            "datasets": datasets, "scientific_execution_authorized": False,
            "limitations": ["Actual frozen OrthoHMM enumerator, not a sorted metadata listing.",
                            "Must recheck immediately before every run and after inference.",
                            "OrthoFinder separately sorts accepted basenames; its fresh copies need independent validation."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("core", "inputs", "runtime", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--verify", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = snapshot(args.core, json.loads(args.inputs.read_text()), json.loads(args.runtime.read_text()))
    result["input_manifest"] = record(args.inputs)
    result["runtime_manifest"] = record(args.runtime)
    if args.verify and result != json.loads(args.verify.read_text()):
        raise ValueError("Frozen native enumeration or provenance changed")
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
