"""Finite-fixture FastME probe using the distance-tree arguments used by STAG."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import subprocess

from Bio import Phylo


def record(path):
    return {"path": str(path.resolve()), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def canonical_tree(path, expected):
    tree = Phylo.read(path, "newick")
    names = [tip.name for tip in tree.get_terminals()]
    if len(names) != len(set(names)) or set(names) != set(expected):
        raise ValueError("Tree leaves differ from input")
    splits = {}
    for node in tree.find_clades():
        if node is tree.root:
            continue
        if node.branch_length is None or not math.isfinite(node.branch_length):
            raise ValueError("Missing or non-finite branch length")
        a = tuple(sorted(tip.name for tip in node.get_terminals()))
        b = tuple(sorted(set(names) - set(a)))
        key = min((a, b), key=lambda x: (len(x), x))
        splits[key] = splits.get(key, 0.0) + node.branch_length
    return [{"split": list(key), "length": value} for key, value in sorted(splits.items())]


def example_matrices(text):
    rows = [line.split() for line in text.splitlines() if line.strip()]
    examples = []
    while rows:
        if len(rows[0]) != 1:
            raise ValueError("Invalid matrix header")
        n = int(rows[0][0])
        if n < 4 or len(rows) < n + 1 or any(len(row) != n + 1 for row in rows[1:n + 1]):
            raise ValueError("Require square matrices")
        block, rows = rows[:n + 1], rows[n + 1:]
        matrix = "\n".join(" ".join(row) for row in block) + "\n"
        examples.append((f"archive_{len(examples)}", matrix, [row[0] for row in block[1:]]))
    if not examples:
        raise ValueError("No matrices")
    return examples


def probe(binary, example, output):
    if output.exists():
        raise FileExistsError(output)
    before, example_record = record(binary), record(example)
    version = subprocess.run([str(binary), "-V"], capture_output=True, text=True, check=True, timeout=30)
    if "FastME 2.1.4" not in version.stdout + version.stderr:
        raise ValueError("Unexpected FastME version")
    examples = example_matrices(example.read_text())
    names = [f"s{i}" for i in range(8)]
    for name in ("additive", "tied"):
        matrix = [[0 if i == j else (2 * (i ^ j).bit_length() if name == "additive" else 2)
                   for j in range(8)] for i in range(8)]
        text = "8\n" + "\n".join(names[i] + " " + " ".join(map(str, row)) for i, row in enumerate(matrix)) + "\n"
        examples.append((name, text, names))
    output.mkdir(parents=True)
    runs = []
    for name, matrix, labels in examples:
        for repeat in range(2):
            directory = output / f"{name}_{repeat}"
            directory.mkdir()
            data, tree = directory / "input.mat", directory / "tree.nwk"
            data.write_text(matrix)
            argv = [str(binary), "-i", str(data), "-o", str(tree), "-I", str(directory / "info.txt"), "-w", "O", "-s", "-n"]
            completed = subprocess.run(argv, capture_output=True, text=True, timeout=60, cwd=directory)
            row = {"fixture": name, "repeat": repeat, "argv": argv, "input": record(data),
                   "exit_code": completed.returncode, "stdout": completed.stdout, "stderr": completed.stderr}
            if completed.returncode == 0:
                row.update(tree=record(tree), canonical_tree=canonical_tree(tree, labels))
            runs.append(row)
    if record(binary) != before or record(example) != example_record:
        raise ValueError("Probe input changed")
    result = {"status": "finite_distance_probe", "machine": platform.machine(), "binary": before,
              "example": example_record, "source": record(Path(__file__)), "runs": runs,
              "limitations": ["Finite example/synthetic matrices and two repeats, not general cross-platform equivalence.",
                              "This tests the STAG command shape, not end-to-end STAG species-tree inference.",
                              "Not a scientific timing result; no inference accuracy is evaluated."]}
    (output / "report.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    if any(row["exit_code"] for row in runs):
        raise RuntimeError("Native FastME failure retained in report")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("binary", "example", "output"):
        parser.add_argument("--" + flag, required=True, type=Path)
    args = parser.parse_args()
    probe(args.binary.resolve(), args.example.resolve(), args.output.resolve())
