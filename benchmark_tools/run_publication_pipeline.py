"""Explicit experimental frozen pipeline with candidate-only canonical ordering.

Run by absolute script path with the recovery interpreter and -I. This does
not modify the installed package or promote the policy to a production default.
"""

import argparse
from contextlib import contextmanager
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import sys


def record(path):
    path = Path(path).absolute()
    return dict(path=str(path), bytes=path.stat().st_size,
                sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def save(path, value):
    with Path(path).open("x") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")


@contextmanager
def candidate_policy(pipeline, receipts):
    from benchmark_tools.candidate_hit_order_policy import canonical_hit_order_v1, POLICY
    from benchmark_tools.probe_ob_canonical_candidates import array_identity
    original = pipeline._expand_phylogeny_candidates

    def expand(output_directory, gene_names, gene_to_species, accuracy_hits, profile="satellite_v1"):
        if receipts or profile != "satellite_v2":
            raise ValueError("Require exactly one satellite_v2 expansion")
        ordered = canonical_hit_order_v1(gene_names, *accuracy_hits)
        receipts.append(dict(policy=POLICY, genes=len(gene_names), hits=len(ordered[0]),
                             original=[array_identity(a) for a in accuracy_hits],
                             canonical=[array_identity(a) for a in ordered], scores_modified=False))
        return original(output_directory, gene_names, gene_to_species, ordered, profile=profile)

    pipeline._expand_phylogeny_candidates = expand
    try:
        yield
    finally:
        pipeline._expand_phylogeny_candidates = original


def arguments(args):
    if args.cpu < 1:
        raise ValueError("CPU count must be positive")
    return [str(args.input), "-o", str(args.output / "inference"), "-c", str(args.cpu),
            "--threads_per_worker", "1", "--search_mode", "builtin", "--clustering", "leiden",
            "--accuracy_profile", "high_sensitivity", "--cpm_resolution", "0.1",
            "--phylogeny", "reconcile", "--species_tree_mode", "infer",
            "--species_tree_rooting", "min_variance", "--phylogeny_candidates", "satellite_v2",
            "--phylogeny_root_rule", "species_overlap", "--phylogeny_pair_rule", "positive_paralogy",
            "--aligner", str(args.aligner), "--tree_builder", str(args.tree_builder)]


def run(args):
    # Resolve the scientific package before making the harness importable.
    import orthohmm
    from orthohmm import orthohmm as pipeline, refinement
    package = Path(orthohmm.__file__).resolve().parent
    if not package.is_relative_to(Path(sys.prefix).resolve()):
        raise ValueError("Require isolated installed OrthoHMM, not the checkout")
    expected = {"orthohmm.py": "2afb89b9dc683e64e58208188f720e07701d4760ff09d1ac3a53c7a8075b84bb",
                "refinement.py": "991f1eb6a5f73d0442529ed19095a34b7c6ba8bff8dfe43a1127e24ec73fb26d"}
    for module in (pipeline, refinement):
        if record(module.__file__)["sha256"] != expected[Path(module.__file__).name]:
            raise ValueError("Changed frozen scientific implementation")
    versions = {n: importlib.metadata.version(n) for n in ("numpy", "igraph", "leidenalg", "DendroPy")}
    if versions != dict(numpy="2.2.6", igraph="1.0.0", leidenalg="0.11.0", DendroPy="5.1.0"):
        raise ValueError("Require validated recovery dependency versions")
    for key, value in dict(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                          MKL_NUM_THREADS="1", PYTHONHASHSEED="0").items():
        if os.environ.get(key) != value:
            raise ValueError("Set reproducible environment before interpreter startup: " + key)
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from benchmark_tools import candidate_hit_order_policy as policy
    if record(policy.__file__)["sha256"] != "cc70b884e40fe73a3c25ef9ae60a2133508127323a11382a8c610e0d192ec811":
        raise ValueError("Changed canonical ordering policy")
    argv = arguments(args)
    inputs = [record(p) for p in sorted(args.input.iterdir()) if p.is_file()]
    if not inputs:
        raise ValueError("Empty input directory")
    sources = [record(p) for p in sorted(package.rglob("*.py"))]
    tools = [record(args.aligner), record(args.tree_builder)]
    args.output.mkdir(parents=True, exist_ok=False)
    save(args.output / "started.json", dict(command=sys.argv, executable=sys.executable,
        package=str(package), versions=versions, scientific_sources=sources, inputs=inputs,
        tools=tools, arguments=argv, policy=record(policy.__file__), source=record(__file__),
        environment={k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
            "MKL_NUM_THREADS", "PYTHONHASHSEED", "MAFFT_BINARIES", "PATH")},
        production_default_changed=False, attempts=1, checkpoint_reuse=False))
    receipts = []
    try:
        (args.output / "inference").mkdir()
        with candidate_policy(pipeline, receipts):
            parsed = pipeline.create_parser().parse_args(argv)
            pipeline.execute(**pipeline.process_args(parsed))
        if len(receipts) != 1:
            raise ValueError("Candidate ordering policy was not applied")
        for item in inputs + sources + tools:
            if record(item["path"]) != item:
                raise ValueError("Input, source or executable changed during inference")
        save(args.output / "complete.json", dict(status="native_complete_pending_scientific_readback",
            policy_application=receipts, started=record(args.output / "started.json"),
            accuracy_evaluated=False, publication_ready=False))
    except BaseException as error:
        save(args.output / "failure.json", dict(type=type(error).__name__, error=str(error),
             policy_application=receipts, retry=False))
        if isinstance(error, SystemExit):
            raise RuntimeError("Frozen pipeline exited before verified completion") from error
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("input", "output", "aligner", "tree-builder"):
        parser.add_argument("--" + name, type=lambda p: Path(p).absolute(), required=True)
    parser.add_argument("--cpu", type=int, default=1)
    run(parser.parse_args())
