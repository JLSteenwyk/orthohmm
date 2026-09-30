"""Read back tiny numerical fixtures from the frozen Python scoring source."""

from __future__ import annotations

import argparse
import ast
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import subprocess

import numpy as np


REVISION = "7f3a9e40dd7e79f842cc2c11fb8b548f9a802806"
PATHS = (
    "orthohmm/accuracy.py", "orthohmm/orthohmm.py",
    "orthohmm/search/profile.py", "orthohmm/search/viterbi.py",
    "orthohmm/search/engine.py", "orthohmm/search/evalue.py",
    "orthohmm/search/matrices.py", "orthohmm/search/msa_profile.py",
    "orthohmm/search/profile_expansion.py",
)


def _definitions(source, names, namespace):
    """Run selected trusted definitions without imports or JIT decorators."""
    tree = ast.parse(source)
    nodes = [node for node in tree.body
             if isinstance(node, (ast.FunctionDef, ast.ClassDef))
             and node.name in names]
    if {node.name for node in nodes} != set(names):
        raise ValueError("frozen definition absent")
    for node in nodes:
        if isinstance(node, ast.FunctionDef):
            node.decorator_list = []
    future = ast.ImportFrom(module="__future__",
                            names=[ast.alias(name="annotations")], level=0)
    module = ast.fix_missing_locations(ast.Module(body=[future, *nodes],
                                                 type_ignores=[]))
    exec(compile(module, "<frozen-selected-definitions>", "exec"), namespace)


def verify(repo: Path):
    blobs = {path: subprocess.check_output(
        ["git", "show", f"{REVISION}:{path}"], cwd=repo) for path in PATHS}
    env = {"np": np, "dataclass": dataclass, "ALPHABET_SIZE": 20,
           "int32": np.int32, "NEG_INF": np.int32(-1000000)}
    exec(compile(blobs["orthohmm/search/matrices.py"], "<frozen-matrices>",
                 "exec"), env)
    _definitions(blobs["orthohmm/search/profile.py"],
                 {"ProfileHMM", "build_profile"}, env)
    _definitions(blobs["orthohmm/search/viterbi.py"],
                 {"_viterbi_score_one"}, env)
    _definitions(blobs["orthohmm/search/evalue.py"], {"estimate_evalue"}, env)
    _definitions(blobs["orthohmm/search/msa_profile.py"],
                 {"encode_alignment", "compute_pssm"}, env)
    matrix = env["get_matrix"]("BLOSUM62")
    background = env["get_background_freqs"]("BLOSUM62")
    seq = np.arange(4, dtype=np.uint8)  # ACDE
    profile = env["build_profile"](seq, 4, matrix, background)
    np.testing.assert_array_equal(profile.transitions, [0, -12, -12, -1, -3, -1, -3])
    np.testing.assert_array_equal(profile.insert_emissions, np.full(20, -1))
    np.testing.assert_array_equal(profile.match_emissions, matrix[:4])
    unknown = env["build_profile"](np.array([20], dtype=np.uint8), 1,
                                    matrix, background)
    np.testing.assert_array_equal(unknown.match_emissions, np.zeros((1, 20)))
    score = int(env["_viterbi_score_one"](profile.match_emissions,
                profile.insert_emissions, profile.transitions, seq, 4, 4, 64))
    assert score == 24
    lam, k = env["get_ka_params"]("BLOSUM62")
    assert (lam, k) == (0.3176, 0.134)
    evalue = float(env["estimate_evalue"](score, 4, 1000, lam, k))
    assert math.isclose(evalue, 0.134 * 4000 * math.exp(-0.3176 * 24),
                        rel_tol=1e-14)
    assert env["estimate_evalue"](0, 4, 1000, lam, k) == 1e10
    encoded = env["encode_alignment"](["AC", "A-"])
    pssm, mask, consensus = env["compute_pssm"](encoded, matrix, background)
    np.testing.assert_array_equal(mask, [True, False])
    pseudo = np.power(2.0, matrix[0].astype(float) / 2) * background
    pseudo /= pseudo.sum()
    counts = np.zeros(20)
    counts[0] = 2
    probability = (counts + 5 * pseudo) / 7
    probability /= probability.sum()
    expected = np.clip(np.rint(2 * np.log2(probability / background)), -128, 127)
    np.testing.assert_array_equal(pssm[0], expected.astype(np.int8))
    assert int(consensus[0]) == int(np.argmax(probability))
    return {
        "revision": REVISION,
        "source": {path: {"bytes": len(blob),
                          "sha256": hashlib.sha256(blob).hexdigest()}
                   for path, blob in blobs.items()},
        "numpy_version": np.__version__,
        "fixtures": {"transitions": profile.transitions.tolist(),
                     "ACDE_self_raw_score": score, "ACDE_normalized_score": 6.0,
                     "ACDE_E_with_1000_target_residues": evalue,
                     "nonpositive_E": 1e10, "half_gap_column_mask": mask.tolist(),
                     "MSA_pseudocount_formula_matches": True,
                     "unknown_query_emissions_zero": True},
        "selected_python_definitions_without_jit": True,
        "native_backend_validation": False, "calibration_established": False,
        "biological_validation": False, "frozen_method_changed": False,
        "publication_ready": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = verify(Path(__file__).resolve().parents[1])
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")


if __name__ == "__main__":
    main()
