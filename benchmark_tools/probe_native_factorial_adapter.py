"""Execute all eight bounded native diagnostics, preserving every attempt."""

import argparse
import itertools
import json
import os
from pathlib import Path
import random
import signal
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


def fixture(directory):
    directory.mkdir(exist_ok=False)
    generator = random.Random(317)
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    families = ["".join(generator.choice(alphabet) for _ in range(160)) for _ in range(6)]
    for species in range(4):
        sequences = []
        for family, ancestral in enumerate(families):
            sequence = list(ancestral)
            for index in range(species * 3, species * 3 + 3):
                sequence[index] = alphabet[(alphabet.index(sequence[index]) + species + 1) % 20]
            sequences.append(("s%d_f%d" % (species, family), "".join(sequence)))
        if species in (0, 1):
            name, sequence = sequences[-1]
            sequences.append((name + "_copy", sequence))
        (directory / ("s%d.fa" % species)).write_text("".join(
            ">%s\n%s\n" % pair for pair in sequences))


def probe(core, baseline, baseline_sha, python, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    python_pin = record(python)
    if record(baseline)["sha256"] != baseline_sha:
        raise ValueError("Wrong diagnostic baseline")
    output.mkdir(parents=True, exist_ok=False)
    fasta = output / "fixture"
    fixture(fasta)
    source = Path(__file__).with_name("diagnose_native_factorial.py")
    inputs = [record(p) for p in sorted(fasta.iterdir())]
    report = dict(schema="native_factorial_adapter_probe_v1", status="running_diagnostics",
        source=record(__file__), diagnostic_source=record(source),
        adapter_source=record(Path(__file__).with_name("native_factorial_adapter.py")),
        baseline=record(baseline), python_executable=python_pin, inputs=inputs, attempts=[],
        fixture_seed=317, diagnostic_only=True, full_dataset_execution_authorized=False,
        accuracy_evaluated=False, publication_ready=False)
    report_path = output / "probe.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    try:
        for values in itertools.product(range(2), repeat=3):
            cell = "p%d_c%d_r%d" % values
            target = output / cell
            env = {k: v for k, v in os.environ.items() if k not in
                   ("PYTHONHOME", "PYTHONPATH", "LD_PRELOAD", "LD_LIBRARY_PATH", "LD_AUDIT")}
            env.update(PYTHONHASHSEED="0", OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
                       PYTHONNOUSERSITE="1", PYTHONDONTWRITEBYTECODE="1", PYTHONPYCACHEPREFIX=str(output / ("unused_cache_" + cell)))
            argv = [str(python), "-B", str(source), "--diagnostic", "--core-root", str(core),
                "--baseline", str(baseline), "--baseline-sha256", baseline_sha,
                "--fasta-directory", str(fasta), "--output-directory", str(target), "--cell", cell]
            attempt = dict(cell=cell, argv=argv, exit_code=None, status="starting")
            report["attempts"].append(attempt)
            report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
            with (output / (cell + ".log")).open("x") as log:
                child = subprocess.Popen(argv, env=env, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
                try:
                    code = child.wait(timeout=300)
                except subprocess.TimeoutExpired:
                    attempt["status"] = "timed_out"
                    os.killpg(child.pid, signal.SIGTERM)
                    try:
                        child.wait(timeout=10)
                    except subprocess.TimeoutExpired:
                        os.killpg(child.pid, signal.SIGKILL)
                        child.wait()
                    raise
                finally:
                    attempt.update(exit_code=child.returncode, log=record(output / (cell + ".log")))
            if code:
                attempt["status"] = "native_failed"
                raise subprocess.CalledProcessError(code, argv)
            attempt["status"] = "native_process_completed"
            path = target / "diagnostic.json"
            diagnostic = json.loads(path.read_text())
            if diagnostic["status"] != "native_factorial_mechanical_diagnostic_complete":
                raise ValueError("Incomplete native diagnostic")
            if values[0] and diagnostic["counts"]["high_sensitivity_profiles"] <= 0:
                raise ValueError("Fixture did not exercise profile construction")
            if values[2] and diagnostic["counts"]["phylogeny_species_tree_families"] <= 0:
                raise ValueError("Fixture did not infer a species tree")
            attempt.update(diagnostic=record(path), stages=diagnostic["stages"], counts=diagnostic["counts"],
                           checkpoint=diagnostic["checkpoint"], native_groups=diagnostic["native_groups"],
                           status="diagnostic_checked")
        manifests = [json.loads(Path(row["checkpoint"]["path"]).read_text()) for row in report["attempts"]]
        if any(manifest != manifests[0] for manifest in manifests[1:]):
            raise ValueError("Factor changes altered initial search checkpoint bytes")
        for pin in inputs + [report[k] for k in ("source", "diagnostic_source", "adapter_source", "baseline", "python_executable")]:
            check(pin)
        report.update(status="eight_native_factorial_diagnostics_complete", initial_search_checkpoints_identical=True,
            genes=manifests[0]["genes"], hits=manifests[0]["hits"],
            limitations=["Mechanical fixtures, not accuracy simulations or full-dataset measurements.",
                "No frozen method/default changes or full production runs; missing configuration costs remain missing.",
                "Identical fixture checkpoints/stage presence do not prove full-dataset equivalence or causal efficiency.",
                "Native adapter handles explicit P/C/R experiments; a prospective resource protocol/launch handoff is still required."])
    except BaseException as error:
        report.update(status="diagnostic_probe_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--diagnostic", action="store_true", required=True)
    for name in ("core-root", "baseline", "python", "output-directory"):
        parser.add_argument("--" + name, required=True, type=Path)
    parser.add_argument("--baseline-sha256", required=True)
    args = parser.parse_args()
    probe(args.core_root.resolve(), args.baseline.resolve(), args.baseline_sha256,
          args.python.absolute(), args.output_directory.absolute())
