"""Reproduce corrected factorial statistics from committed counts, not raw scoring."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

MODULES = ("reproduce_corrected_factorial.py", "bootstrap_qfo_factorial.py", "bootstrap_qfo_swiss_stages.py",
    "audit_qfo_swiss_counts.py", "prepare_ob_candidate_neighborhood.py", "report_qfo_recovered_stages.py",
    "validate_qfo_native_assessment.py", "qfo_summarize_scores.py")
COUNTS = "qfo_corrected_factorial_complete_20260919/swiss_counts.json"
RESULT = "qfo_corrected_factorial_complete_20260919/swiss_bootstrap.json"
HASHES = {
    COUNTS: "c9d8bf02ef6f287c56d073fa61f165c982caf834a94fc6e44e6589f03e46eba6",
    RESULT: "f777dead1294b3810c0479877aade7fb4411c0544696f671226ff198432e9211",
    "QFO_FACTORIAL_PROTOCOL_20260917.md": "f8946e12cefcf84abbee0fb9492f240c05508e045efe00a3006304d34c1fd115",
    "QFO_CORRECTED_RELEASE_PROTOCOL_20260917.md": "b3603bc9b51ce1708f49a0ee816eb02d3cb52c66fb00f07eeb0d66f840a03e1e"}
NUMERICAL = {"alpha", "comparisons", "families", "multiplicity_endpoints", "numpy_version",
    "point_estimates", "quantile_method", "replicates", "rng", "seed", "units"}
COMMON = {"status", "publication_ready", "limitations"}
PROVENANCE = {"checked_inputs", "corrected_protocol", "counts", "helpers", "input_release", "protocol", "source"}
ENGINE_SHA = "be09876ae4de31b923818385bd3d40d8d5e7df177215fc45a2101c7d30d1c1fa"


def identity(path):
    content = path.read_bytes()
    return dict(path=str(path.resolve()), bytes=len(content), sha256=hashlib.sha256(content).hexdigest())


def compare(expected, observed):
    if (set(expected) != NUMERICAL | COMMON | PROVENANCE or set(observed) != NUMERICAL | COMMON
            or expected["status"] != "paired_corrected_qfo_factorial_swiss_intervals"
            or expected["input_release"] != "QfO 2020_04 corrected UP000008143"
            or expected["publication_ready"] is not False or observed["publication_ready"] is not False
            or observed["status"] != "paired_qfo_factorial_swiss_intervals"):
        raise ValueError("Unexpected corrected statistics or admission schema")
    if {k: expected[k] for k in NUMERICAL} != {k: observed[k] for k in NUMERICAL}:
        raise ValueError("Reproduced scientific statistics differ")


def worker(export, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    sys.path.insert(0, str(export))
    from benchmark_tools.bootstrap_qfo_factorial import bootstrap
    results = export / "benchmark_tools/results"
    for name, digest in HASHES.items():
        if identity(results / name)["sha256"] != digest:
            raise ValueError("Changed committed analysis input")
    if identity(export / "benchmark_tools/bootstrap_qfo_factorial.py")["sha256"] != ENGINE_SHA:
        raise ValueError("Changed numerical engine")
    expected = json.loads((results / RESULT).read_text())
    counts = json.loads((results / COUNTS).read_text())
    if (counts["status"] != "corrected_qfo_factorial_swiss_counts_verified"
            or counts["uncertainty_admitted"] is not False or counts["publication_ready"] is not False):
        raise ValueError("Require unresampled corrected counts")
    # Identical schema adaptation to the admitted corrected wrapper; no count or arithmetic changes.
    observed = bootstrap({**counts, "status": "qfo_factorial_swiss_counts_verified"}, replicates=100000, seed=20260922)
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
    command = [str(python), "-I", "-B", str(export / "benchmark_tools/reproduce_corrected_factorial.py"),
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
    return dict(status="relocated_corrected_factorial_statistics_reproduced", source_commit=commit,
        runner=identity(Path(__file__)), exported=exported, environment=environment, command=command,
        exit_code=completed.returncode, scientific_results_exact_match=True,
        raw_source_admission_repeated=False, publication_ready=False,
        generated=[identity(p) for p in sorted((output / "reproduced").iterdir())],
        limitations=["Statistics-only reproduction from pinned family counts, not new raw-source admission.",
            "The shared numerical engine retains its generic status; these inputs are corrected-release only.",
            "Native inference, conversion, reference construction and official scoring were not rerun.",
            "Historical absolute paths are provenance, not accessed inputs.",
            "Same-host isolated execution is not cross-platform validation, full release or rights clearance."])


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
