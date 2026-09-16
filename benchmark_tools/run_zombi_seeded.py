"""Run pinned Zombi with explicit reproducible Pyvolve family seeds."""

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import runpy
import subprocess
import sys


ZOMBI_COMMIT = "8db13ee4ba007f46c17f38586d31e5aa617c1647"


def family_seed(seed, seqfile):
    if not isinstance(seed, int) or seed <= 0 or not seqfile:
        raise ValueError("Positive master seed and sequence filename required")
    # Basename makes independent output directories reproduce the same family.
    key = f"{seed}:{Path(seqfile).name}".encode("ascii")
    return int.from_bytes(hashlib.sha256(key).digest()[:16], "big")


def seeded_evolver(original, seed):
    class SeededEvolver(original):
        def __call__(self, **kwargs):
            if kwargs.get("seed") is None:
                kwargs["seed"] = family_seed(seed, kwargs.get("seqfile"))
            return super().__call__(**kwargs)
    return SeededEvolver


def family_lengths(seed, families):
    if type(seed) is not int or seed <= 0:
        raise ValueError("Positive integer length seed required")
    names = list(families)
    if not names or len(set(names)) != len(names) or any(
            not isinstance(n, str) or not n.isascii() or not n.isdigit() or str(int(n)) != n or int(n) <= 0
            for n in names):
        raise ValueError("Unique canonical positive-integer family IDs required")
    return {name: 100 + int.from_bytes(hashlib.sha256(
        f"orthohmm-sim-length-v2:{seed}:{name}".encode("ascii")).digest()[:8], "big") % 401
        for name in sorted(names, key=int)}


def length_aware_simulator(original, lengths):
    class FamilyLengthSimulator(original):
        def run(self, tree_file, sequences_folder):
            name = Path(tree_file).name
            suffix = "_completetree.nwk"
            if not name.endswith(suffix) or name[:-len(suffix)] not in lengths:
                raise ValueError("Tree has no frozen family length")
            previous = self.size
            self.size = lengths[name[:-len(suffix)]]
            try:
                return super().run(tree_file, sequences_folder)
            finally:
                self.size = previous
    return FamilyLengthSimulator


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--mode", choices=["T", "G", "S"], required=True)
    parser.add_argument("--parameters", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--family-lengths", type=Path)
    args = parser.parse_args()
    source, parameters, output = [p.resolve() for p in (args.source, args.parameters, args.output)]
    if any(any(c.isspace() for c in str(p)) for p in (source, parameters, output)):
        raise ValueError("Upstream shell operations require whitespace-free paths")
    if (output / args.mode).exists():
        raise FileExistsError("Refusing to replace an existing Zombi stage")
    commit = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
    if commit != ZOMBI_COMMIT:
        raise ValueError("Unexpected Zombi revision")
    subprocess.run(["git", "-C", str(source), "diff", "--exit-code", "HEAD", "--", "*.py"], check=True)
    if importlib.metadata.version("Pyvolve") != "1.1.0":
        raise ValueError("Pyvolve 1.1.0 required")
    sys.path.insert(0, str(source))
    import AuxiliarFunctions as af
    import pyvolve
    parsed = af.read_parameters(str(parameters))
    if args.seed <= 0 or int(parsed["SEED"]) != args.seed:
        raise ValueError("Parameter seed must equal positive master seed")
    if args.family_lengths is not None:
        if args.mode != "S" or parsed["SEQUENCE"] != "amino-acid":
            raise ValueError("Family lengths require ordinary amino-acid S mode")
        mapping = json.loads(args.family_lengths.read_text())
        trees = sorted((output / "G/Gene_trees").glob("*_completetree.nwk"))
        families = [p.name.removesuffix("_completetree.nwk") for p in trees]
        expected = family_lengths(args.seed, families)
        if mapping != {"schema_version": 1, "seed": args.seed, "lengths": expected}:
            raise ValueError("Family length mapping differs from frozen rule or tree set")
        import SequenceSimulator as sequence_module
        sequence_module.SequenceSimulator = length_aware_simulator(sequence_module.SequenceSimulator, expected)
    pyvolve.Evolver = seeded_evolver(pyvolve.Evolver, args.seed)
    sys.argv = [str(source / "Zombi.py"), args.mode, str(parameters), str(output)]
    runpy.run_path(str(source / "Zombi.py"), run_name="__main__")


if __name__ == "__main__":
    main()
