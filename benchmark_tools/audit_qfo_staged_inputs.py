"""Audit retained comparator staging files against the original QfO manifest."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.inventory_swiss_sequences import PREPARED_SHA
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check


STAGES = {
    "sonicparanoid": ("qfo_benchmark/results/sonicparanoid/run/input", ".fasta"),
    "proteinortho": ("qfo_benchmark/results/proteinortho/run/input", ".fasta"),
    "fastoma": ("qfo_benchmark/results/fastoma/run/input/proteome", ".fa"),
}


def compare_stage(inputs, directory, extension):
    expected = {}
    for item in inputs:
        name = Path(item["path"]).stem + extension
        if name in expected:
            raise ValueError("Duplicate expected staged filename")
        expected[name] = item
    if not expected:
        raise ValueError("Empty original inventory")
    # Check all common FASTA suffixes so an extra proteome cannot hide under another suffix.
    inventory = {p.name for p in directory.iterdir()
                 if p.suffix.lower() in {".fasta", ".fa", ".faa", ".fas", ".fna"}}
    if inventory != set(expected):
        raise ValueError("Staged FASTA inventory differs from frozen originals")
    rows = []
    for name, original in sorted(expected.items()):
        check(original)
        path = directory / name
        if not path.is_file():
            raise ValueError("Staged FASTA is not a file")
        native = record(path)
        rows.append({"original": original, "staged": native, "staged_name": name,
                     "is_symlink": path.is_symlink(),
                     "same_resolved_path": native["path"] == original["path"],
                     "identical_file_bytes": all(native[k] == original[k] for k in ("bytes", "sha256"))})
    for row in rows:
        check(row["original"])
        if record(directory / row["staged_name"]) != row["staged"]:
            raise ValueError("Staged file changed during audit")
    final_inventory = {p.name for p in directory.iterdir()
                       if p.suffix.lower() in {".fasta", ".fa", ".faa", ".fas", ".fna"}}
    if final_inventory != inventory:
        raise ValueError("Staged inventory changed during audit")
    return {"directory": str(directory.resolve()), "files": rows, "proteomes": len(rows),
            "identical_files": sum(row["identical_file_bytes"] for row in rows),
            "all_files_identical": all(row["identical_file_bytes"] for row in rows)}


def audit(root, prepared):
    identity = record(prepared)
    if identity["sha256"] != PREPARED_SHA:
        raise ValueError("Changed frozen original input manifest")
    inputs = json.loads(prepared.read_text())["input_fastas"]
    methods = {name: compare_stage(inputs, root / directory, extension)
               for name, (directory, extension) in STAGES.items()}
    check(identity)
    return {"status": "retained_comparator_staged_fastas_compared",
            "source": record(__file__), "prepared": identity, "methods": methods,
            "limitations": ["Current retained staging copies, not authenticated historical execution traces.",
                            "No audit of downstream native sequence transformations, binary search databases or predictions.",
                            "Original-release parity does not establish corrected-release compatibility.",
                            "No rerun, input normalization, scoring change or claim of end-to-end publication readiness."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "prepared", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.root, args.prepared)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
