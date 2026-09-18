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
DOMAIN_MODULES = ("analyze_swiss_domain_strata.py", "inventory_swiss_annotations.py",
                  "snapshot_orthohmm_input_order.py", "plot_swiss_domain_strata.py")
DOMAIN_DATA = ("swiss_domain_annotation_inventory_20260917.json", "SWISS_DOMAIN_STRATA_PROTOCOL_20260917.md",
               "swiss_domain_strata_results_20260917.json", "SWISS_DOMAIN_STRATA_RESULTS_20260917.md")


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


def compare_domain_results(expected, observed):
    provenance = {"inputs", "source", "helpers"}
    if {k: v for k, v in expected.items() if k not in provenance} != {
            k: v for k, v in observed.items() if k not in provenance}:
        raise ValueError("Relocated domain-stratified scientific results differ")
    for key in ("inputs", "helpers"):
        if [(r["bytes"], r["sha256"]) for r in expected[key]] != [
                (r["bytes"], r["sha256"]) for r in observed[key]]:
            raise ValueError("Relocated domain input or helper bytes differ")
    if any(expected["source"][k] != observed["source"][k] for k in ("bytes", "sha256")):
        raise ValueError("Relocated domain source differs")


def reproduce(repo, revision, python, output, include_domain_strata=False):
    if output.exists():
        raise FileExistsError(output)
    commit = subprocess.check_output(["git", "rev-parse", "--verify", revision + "^{commit}"], cwd=repo, text=True).strip()
    output.mkdir(parents=True)
    export = output / "source"
    modules = MODULES + (DOMAIN_MODULES if include_domain_strata else ())
    data = DATA + (DOMAIN_DATA if include_domain_strata else ())
    paths = ["LICENSE.md", *["benchmark_tools/" + p for p in modules],
             *["benchmark_tools/results/" + p for p in data]]
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
    if include_domain_strata:
        command = [str(python), "-I", str(export / "benchmark_tools/analyze_swiss_domain_strata.py"),
                   "--counts", str(results / DATA[0]), "--annotations", str(results / DOMAIN_DATA[0]),
                   "--protocol", str(results / DOMAIN_DATA[1]), "--output", str(generated / "domain_strata.json"),
                   "--markdown", str(generated / "domain_strata.md")]
        run = subprocess.run(command, cwd=export, env=env, capture_output=True, text=True, check=True)
        logs.append({"command": command, "stdout": run.stdout, "stderr": run.stderr, "returncode": run.returncode})
        compare_domain_results(json.loads((results / DOMAIN_DATA[2]).read_text()),
                               json.loads((generated / "domain_strata.json").read_text()))
        if (results / DOMAIN_DATA[3]).read_bytes() != (generated / "domain_strata.md").read_bytes():
            raise ValueError("Generated domain-stratified table differs")
        command = [str(python), "-I", str(export / "benchmark_tools/plot_swiss_domain_strata.py"),
                   "--results", str(results / DOMAIN_DATA[2]), "--output", str(generated / "domain_figures")]
        run = subprocess.run(command, cwd=export, env=env, capture_output=True, text=True, check=True)
        logs.append({"command": command, "stdout": run.stdout, "stderr": run.stderr, "returncode": run.returncode})
    for item in exported:
        if identity(Path(item["path"])) != item:
            raise ValueError("Export changed during reproduction")
    return {"status": "relocated_swiss_statistics_and_plot_workflow_reproduced", "source_commit": commit,
            "runner": identity(Path(__file__)), "exported": exported, "environment": environment, "execution": logs,
            "scientific_results_exact_match": True, "markdown_byte_match": True,
            "domain_strata": {"included": include_domain_strata,
                "scientific_results_exact_match": True if include_domain_strata else None,
                "markdown_byte_match": True if include_domain_strata else None,
                "annotation_regeneration": False},
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
    parser.add_argument("--include-domain-strata", action="store_true")
    args = parser.parse_args()
    if args.report.exists():
        raise FileExistsError(args.report)
    result = reproduce(args.repo.resolve(), args.revision, args.python.absolute(), args.output.resolve(), args.include_domain_strata)
    args.report.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
