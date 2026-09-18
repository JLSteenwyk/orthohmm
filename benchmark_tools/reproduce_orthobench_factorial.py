"""Reproduce committed OrthoBench factorial statistics outside the checkout."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

MODULES = ("reproduce_orthobench_factorial.py", "bootstrap_orthobench_factorial.py",
           "bootstrap_orthobench.py", "orthobench_stage_diagnostics.py", "score_orthobench_partition.py")
DATA = "benchmark_tools/results/orthobench_factorial_results_20260916.json"
METADATA = {"scores", "official_scores", "native_validation", "scheduler", "accounting_raw", "predictions",
            "references", "coverage_resources", "official_scorer", "assembler", "publication_ready", "timing_scope"}


def identity(path):
    content = path.read_bytes()
    return dict(path=str(path.resolve()), bytes=len(content), sha256=hashlib.sha256(content).hexdigest())


def compare(expected, observed):
    if set(expected) != set(observed) | METADATA or expected["publication_ready"] is not False:
        raise ValueError("Unexpected retained factorial schema")
    if {k: v for k, v in expected.items() if k not in METADATA} != observed:
        raise ValueError("Reproduced scientific statistics differ")


def worker(export, output):
    sys.path.insert(0, str(export))
    from benchmark_tools.bootstrap_orthobench_factorial import CELLS, factorial_bootstrap, render_report
    expected = json.loads((export / DATA).read_text())
    if set(expected["scores"]) != set(CELLS) or expected["failed_cells"]:
        raise ValueError("Require the complete retained factorial")
    cells = {cell: {"status": "complete", "refog_records": expected["scores"][cell]["refog_records"]} for cell in CELLS}
    result = factorial_bootstrap(cells)
    compare(expected, result)
    output.mkdir(exist_ok=False)
    (output / "statistics.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    (output / "statistics.md").write_text(render_report(result))


def reproduce(repo, revision, python, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    commit = subprocess.check_output(["git", "rev-parse", "--verify", revision + "^{commit}"], cwd=repo, text=True).strip()
    export = output / "source"
    export.mkdir(parents=True)
    paths = ["LICENSE.md", DATA, "benchmark_tools/swiss_analysis_requirements.txt",
             "benchmark_tools/swiss_analysis_requirements.lock", *["benchmark_tools/" + n for n in MODULES]]
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
        raise ValueError("Require an analysis virtual environment")
    command = [str(python), "-I", "-B", str(export / "benchmark_tools/reproduce_orthobench_factorial.py"),
               "--worker", str(export), "--output", str(output / "reproduced")]
    completed = subprocess.run(command, cwd=export, env=env, capture_output=True, text=True)
    (output / "worker.stdout").write_text(completed.stdout)
    (output / "worker.stderr").write_text(completed.stderr)
    completed.check_returncode()
    observed = json.loads((output / "reproduced/statistics.json").read_text())
    compare(json.loads((export / DATA).read_text()), observed)
    for item in exported:
        if identity(Path(item["path"])) != item:
            raise ValueError("Export changed during reproduction")
    return dict(status="relocated_orthobench_factorial_statistics_reproduced", source_commit=commit,
        runner=identity(Path(__file__)), exported=exported, environment=environment, command=command,
        exit_code=completed.returncode, scientific_results_exact_match=True, publication_ready=False,
        generated=[identity(p) for p in sorted((output / "reproduced").iterdir())],
        limitations=["Statistics-only reproduction from retained family sufficient statistics.",
            "Native inference, conversion, reference construction and official scoring were not rerun.",
            "Historical absolute paths are provenance only and are not accessed as analysis inputs.",
            "Isolated execution on the same host is not cross-platform or hermetic OS validation.",
            "Export includes installation locks but does not install or prove completeness of all inference environments."])


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
