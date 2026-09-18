"""Export a committed SwissTrees analysis and verify relocated reproduction."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

MODULES = ("bootstrap_qfo_swiss_comparators.py", "audit_qfo_swiss_comparators.py",
           "audit_qfo_swiss_counts.py", "prepare_ob_candidate_neighborhood.py",
           "report_qfo_recovered_stages.py", "validate_qfo_native_assessment.py",
           "qfo_summarize_scores.py", "bootstrap_qfo_swiss_stages.py", "plot_qfo_swiss_comparators.py")
DATA = ("qfo_swiss_comparator_counts_20260917.json", "QFO_SWISS_COMPARATOR_UNCERTAINTY_PROTOCOL_20260917.md",
        "qfo_swiss_comparator_intervals_20260917.json", "QFO_SWISS_COMPARATOR_INTERVALS_20260917.md")


def identity(path):
    content = path.read_bytes()
    return {"path": str(path.resolve()), "bytes": len(content), "sha256": hashlib.sha256(content).hexdigest()}


def compare_results(expected, observed):
    provenance = {"counts", "protocol", "source", "helpers"}
    if {k: v for k, v in expected.items() if k not in provenance} != {
            k: v for k, v in observed.items() if k not in provenance}:
        raise ValueError("Relocated scientific results differ")
    for key in ("counts", "protocol", "source"):
        if any(expected[key][k] != observed[key][k] for k in ("bytes", "sha256")):
            raise ValueError("Relocated source or input bytes differ")
    if [(r["bytes"], r["sha256"]) for r in expected["helpers"]] != [
            (r["bytes"], r["sha256"]) for r in observed["helpers"]]:
        raise ValueError("Relocated helper bytes differ")


def reproduce(repo, revision, python, output):
    if output.exists():
        raise FileExistsError(output)
    commit = subprocess.check_output(["git", "rev-parse", "--verify", revision + "^{commit}"], cwd=repo, text=True).strip()
    output.mkdir(parents=True)
    export = output / "source"
    paths = ["LICENSE.md", *["benchmark_tools/" + p for p in MODULES],
             *["benchmark_tools/results/" + p for p in DATA]]
    exported = []
    for name in paths:
        content = subprocess.check_output(["git", "show", commit + ":" + name], cwd=repo)
        path = export / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)
        exported.append(identity(path))
    env = {k: v for k, v in os.environ.items() if k not in ("PYTHONPATH", "PYTHONHOME", "LD_PRELOAD", "LD_LIBRARY_PATH", "LD_AUDIT")}
    env.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", PYTHONHASHSEED="0", PYTHONNOUSERSITE="1",
               MPLCONFIGDIR=str(output / "matplotlib_cache"))
    environment = json.loads(subprocess.check_output([str(python), "-I", "-c",
        "import sys, importlib.metadata as m, json; print(json.dumps({'python':sys.version,'prefix':sys.prefix,'base_prefix':sys.base_prefix,'packages':{d.metadata['Name']:d.version for d in m.distributions()}}))"], env=env, text=True))
    if environment["prefix"] == environment["base_prefix"]:
        raise ValueError("Require an isolated analysis virtual environment")
    results = export / "benchmark_tools/results"
    generated = output / "reproduced"
    generated.mkdir()
    commands = [[str(python), "-I", str(export / "benchmark_tools/bootstrap_qfo_swiss_comparators.py"),
                 "--counts", str(results / DATA[0]), "--protocol", str(results / DATA[1]),
                 "--output", str(generated / "intervals.json"), "--markdown", str(generated / "intervals.md")]]
    logs = []
    for command in commands:
        run = subprocess.run(command, cwd=export, env=env, capture_output=True, text=True, check=True)
        logs.append({"command": command, "stdout": run.stdout, "stderr": run.stderr, "returncode": run.returncode})
    expected = json.loads((results / DATA[2]).read_text())
    observed = json.loads((generated / "intervals.json").read_text())
    compare_results(expected, observed)
    if (results / DATA[3]).read_bytes() != (generated / "intervals.md").read_bytes():
        raise ValueError("Generated table differs")
    # Plot the pinned report only after its complete scientific content reproduced.
    command = [str(python), "-I", str(export / "benchmark_tools/plot_qfo_swiss_comparators.py"),
               "--results", str(results / DATA[2]), "--output", str(generated / "figures")]
    run = subprocess.run(command, cwd=export, env=env, capture_output=True, text=True, check=True)
    logs.append({"command": command, "stdout": run.stdout, "stderr": run.stderr, "returncode": run.returncode})
    for item in exported:
        if identity(Path(item["path"])) != item:
            raise ValueError("Export changed during reproduction")
    return {"status": "relocated_swiss_statistics_and_plot_workflow_reproduced", "source_commit": commit,
            "runner": identity(Path(__file__)), "exported": exported, "environment": environment, "execution": logs,
            "scientific_results_exact_match": True, "markdown_byte_match": True,
            "generated": [identity(p) for p in sorted(generated.rglob("*")) if p.is_file()],
            "limitations": ["Statistical/figure workflow only; native inference and raw-QfO scoring are not reproduced.",
                            "Historical absolute paths inside reports remain provenance, not accessed inputs.",
                            "Figure file bytes can vary with renderer metadata; no bitwise image-equivalence claim.",
                            "An isolated environment on the same host is not cross-platform validation.",
                            "This is not the complete publication archive, license clearance or a versioned public release."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument("--revision", required=True)
    parser.add_argument("--python", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    if args.report.exists():
        raise FileExistsError(args.report)
    result = reproduce(args.repo.resolve(), args.revision, args.python.absolute(), args.output.resolve())
    args.report.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
