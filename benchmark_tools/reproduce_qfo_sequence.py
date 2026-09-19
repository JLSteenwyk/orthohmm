"""Recompute committed corrected sequence-control statistics outside the checkout."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

MODULES = ("reproduce_qfo_sequence.py", "bootstrap_qfo_sequence.py", "bootstrap_qfo_swiss_stages.py",
    "audit_qfo_swiss_counts.py", "prepare_ob_candidate_neighborhood.py", "report_qfo_recovered_stages.py",
    "validate_qfo_native_assessment.py", "qfo_summarize_scores.py")
COUNTS = "qfo_sequence_swiss_counts_20260918.json"
RESULT = "qfo_sequence_swiss_bootstrap_20260918.json"
PROTOCOL = "QFO_SEQUENCE_UNCERTAINTY_PROTOCOL_20260918.md"
HASHES = {
    COUNTS: "3f04b05435f5967630e6e109e91b9e51d8618d940eb594e07681f5de84d0ae8c",
    RESULT: "e4e4696e6f090d7d83b99e5190bd21606a980fac294c62cdf260582d822d89f3",
    PROTOCOL: "e853c2f9fbfe469e34fc6ebd0d00692c9fd9c91327750931e01131c0adb5e277"}
NUMERICAL = {"alpha", "comparisons", "families", "multiplicity_endpoints", "numpy_version",
    "point_estimates", "quantile_method", "replicates", "rng", "seed", "units"}
COMMON = {"status", "publication_ready", "uncertainty_admitted", "limitations"}
PROVENANCE = {"admission_inventory", "baseline_audit", "checked_inputs", "counts", "helpers", "protocol", "source"}


def identity(path):
    content = path.read_bytes()
    return dict(path=str(path.resolve()), bytes=len(content), sha256=hashlib.sha256(content).hexdigest())


def compare(expected, observed):
    if (set(expected) != NUMERICAL | COMMON | PROVENANCE or set(observed) != NUMERICAL | COMMON
            or expected["status"] != "paired_corrected_sequence_swiss_intervals"
            or expected["uncertainty_admitted"] is not True or expected["publication_ready"] is not False
            or observed["status"] != "corrected_sequence_swiss_intervals_pending_source_admission"
            or observed["uncertainty_admitted"] is not False or observed["publication_ready"] is not False):
        raise ValueError("Unexpected statistics/admission schema")
    if {k: expected[k] for k in NUMERICAL} != {k: observed[k] for k in NUMERICAL}:
        raise ValueError("Reproduced scientific statistics differ")


def worker(export, output):
    sys.path.insert(0, str(export))
    from benchmark_tools.bootstrap_qfo_sequence import bootstrap
    results = export / "benchmark_tools/results"
    for name, digest in HASHES.items():
        if identity(results / name)["sha256"] != digest:
            raise ValueError("Changed committed analysis input")
    expected = json.loads((results / RESULT).read_text())
    counts = json.loads((results / COUNTS).read_text())
    observed = bootstrap(counts, replicates=100000, seed=20260923)
    compare(expected, observed)
    output.mkdir(exist_ok=False)
    (output / "statistics.json").write_text(json.dumps(observed, indent=2, sort_keys=True) + "\n")


def reproduce(repo, revision, python, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    commit = subprocess.check_output(["git", "rev-parse", "--verify", revision + "^{commit}"], cwd=repo, text=True).strip()
    export = output / "source"
    export.mkdir(parents=True)
    paths = ["LICENSE.md", "benchmark_tools/swiss_analysis_requirements.txt",
        "benchmark_tools/swiss_analysis_requirements.lock", *["benchmark_tools/" + n for n in MODULES],
        *["benchmark_tools/results/" + n for n in HASHES]]
    exported = []
    for name in paths:
        content = subprocess.check_output(["git", "show", commit + ":" + name], cwd=repo)
        path = export / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)
        exported.append(identity(path))
    env = {k: v for k, v in os.environ.items() if k not in
        ("PYTHONPATH", "PYTHONHOME", "LD_PRELOAD", "LD_LIBRARY_PATH", "LD_AUDIT")}
    env.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", PYTHONDONTWRITEBYTECODE="1")
    environment = json.loads(subprocess.check_output([str(python), "-I", "-c",
        "import sys,json,importlib.metadata as m; print(json.dumps(dict(python=sys.version,prefix=sys.prefix,base_prefix=sys.base_prefix,packages={d.metadata['Name']:d.version for d in m.distributions()})))"],
        env=env, text=True))
    if environment["prefix"] == environment["base_prefix"]:
        raise ValueError("Require an isolated analysis virtual environment")
    command = [str(python), "-I", "-B", str(export / "benchmark_tools/reproduce_qfo_sequence.py"),
        "--worker", str(export), "--output", str(output / "reproduced")]
    completed = subprocess.run(command, cwd=export, env=env, capture_output=True, text=True)
    (output / "worker.stdout").write_text(completed.stdout)
    (output / "worker.stderr").write_text(completed.stderr)
    completed.check_returncode()
    observed = json.loads((output / "reproduced/statistics.json").read_text())
    compare(json.loads((export / "benchmark_tools/results" / RESULT).read_text()), observed)
    for item in exported:
        if identity(Path(item["path"])) != item:
            raise ValueError("Export changed during reproduction")
    return dict(status="relocated_corrected_sequence_statistics_reproduced", source_commit=commit,
        runner=identity(Path(__file__)), exported=exported, environment=environment, command=command,
        exit_code=completed.returncode, scientific_results_exact_match=True, publication_ready=False,
        generated=[identity(p) for p in sorted((output / "reproduced").iterdir())],
        limitations=["Statistics-only reproduction from pinned audited family counts, not new source admission.",
            "Native inference, pair conversion, reference construction and official scoring were not rerun.",
            "Historical absolute paths are retained provenance, not accessed analysis inputs.",
            "Isolated execution on the same host is not cross-platform or hermetic OS validation.",
            "No full publication release, raw-data redistribution clearance or method superiority is established."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", type=Path)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument("--revision")
    parser.add_argument("--python", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    if args.worker:
        worker(args.worker.resolve(), args.output.absolute())
    else:
        if not args.revision or not args.python or not args.report:
            parser.error("Require --revision, --python and --report")
        if args.report.exists() or args.report.is_symlink():
            raise FileExistsError(args.report)
        result = reproduce(args.repo.resolve(), args.revision, args.python.absolute(), args.output.absolute())
        with args.report.open("x") as stream:
            json.dump(result, stream, indent=2, sort_keys=True)
            stream.write("\n")
