"""Prespecified eight-species length-adapter smoke test; never score accuracy."""

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.benchmark_production import file_record
from benchmark_tools.prepare_simulation_panel import DEFAULTS, parameter_sets
from benchmark_tools.run_zombi_seeded import family_lengths
from benchmark_tools.validate_simulation_outputs import validate_orthofinder
from benchmark_tools.zombi_truth import validate_run

SEED = 20260918


def run(source, output, orthofinder):
    output.mkdir(parents=True, exist_ok=False)
    sys.path.insert(0, str(source))
    import AuxiliarFunctions as af
    defaults = {mode: af.read_parameters(str(source / "Parameters" / name)) for mode, name in DEFAULTS.items()}
    parameters = parameter_sets(defaults, SEED, "baseline", allowed_seeds=(SEED,))
    lengths = family_lengths(SEED, [str(i) for i in range(1, 101)])
    length_path = output / "family_lengths.json"
    length_path.write_text(json.dumps({"schema_version": 1, "seed": SEED, "lengths": lengths}, indent=2, sort_keys=True) + "\n")
    driver = Path(__file__).with_name("run_zombi_seeded.py").resolve()
    env = os.environ.copy()
    env.update(PYTHONHASHSEED="0", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    report = {"schema_version": 1, "status": "running", "seed": SEED, "scientific_panel_member": False,
              "accuracy_computed": False, "runs": {}, "python": sys.version,
              "family_lengths": file_record(length_path, output),
              "protocol": file_record(Path(__file__).parent / "results/PUBLICATION_VARIABLE_LENGTH_PROTOCOL_20260916.md", Path(__file__).parent),
              "sources": [file_record(Path(__file__).with_name(n), Path(__file__).parent) for n in
                          ("smoke_zombi_variable_lengths.py", "run_zombi_seeded.py", "prepare_simulation_panel.py",
                           "zombi_truth.py", "validate_simulation_outputs.py")],
              "native_sources": [file_record(p, source) for p in sorted(source.glob("*.py"))]}
    report_path = output / "report.json"

    def save():
        report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

    try:
        for label in ("repeat_a", "repeat_b"):
            native = output / label
            param_dir = output / (label + "_parameters")
            param_dir.mkdir()
            record = report["runs"][label] = {"stages": []}
            for mode, values in parameters.items():
                path = param_dir / DEFAULTS[mode]
                path.write_text("".join(f"{k}\t{v}\n" for k, v in sorted(values.items())))
                argv = [sys.executable, str(driver), "--source", str(source), "--mode", mode,
                        "--parameters", str(path), "--output", str(native), "--seed", str(SEED)]
                if mode == "S":
                    argv.extend(["--family-lengths", str(length_path)])
                stage = {"stage": mode, "argv": argv, "parameters": file_record(path, output), "status": "running"}
                record["stages"].append(stage)
                save()
                with (output / f"{label}_{mode}.log").open("w") as log:
                    result = subprocess.run(argv, env=env, stdout=log, stderr=subprocess.STDOUT, check=False)
                stage.update(exit_code=result.returncode, status="complete" if result.returncode == 0 else "failed")
                if result.returncode:
                    raise RuntimeError(f"Native {label} {mode} failed")
            truth, sequences = validate_run(native)
            observed = {len(seq) for _, seq in sequences.values()}
            if len(observed) < 2 or any(len(seq) != lengths[gene.split("__", 1)[0][1:]]
                                         for gene, (_, seq) in sequences.items()):
                raise ValueError("Exported family sequence lengths differ from assignment")
            record.update(truth_validated=True, extant_genes=len(sequences),
                          unique_extant_lengths=len(observed), length_min=min(observed), length_max=max(observed),
                          native_products=[file_record(p, native) for p in sorted(native.rglob("*")) if p.is_file()])
            prepared = output / (label + "_prepared")
            argv = [sys.executable, str(driver.with_name("zombi_truth.py")), "--run", str(native), "--output", str(prepared)]
            with (output / f"{label}_truth.log").open("w") as log:
                subprocess.run(argv, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
            record["truth_export_command"] = argv
            save()
        if report["runs"]["repeat_a"]["native_products"] != report["runs"]["repeat_b"]["native_products"]:
            raise ValueError("Repeat native simulation products differ")
        report["same_seed_all_native_products_identical"] = True
        of_root = output / "orthofinder"
        copied = of_root / "input"
        copied.mkdir(parents=True)
        input_records = []
        for path in sorted((output / "repeat_a_prepared/input").glob("*.fasta")):
            shutil.copy2(path, copied / path.name)
            input_records.append(dict(file_record(path, path.parent), absolute_path=str(path)))
        argv = [str(orthofinder), "-f", str(copied), "-t", "4", "-a", "4", "-S", "diamond"]
        config = {"argv": argv, "output": str(of_root), "copy_inputs_to": str(copied)}
        of = report["orthofinder"] = {"argv": argv, "status": "running",
             "entrypoint": file_record(orthofinder, orthofinder.parent), "inputs": input_records}
        save()
        with (output / "orthofinder.log").open("w") as log:
            result = subprocess.run(["/usr/bin/time", "-v", "-o", str(output / "orthofinder.time.log"), *argv],
                                    env=env, stdout=log, stderr=subprocess.STDOUT, check=False)
        of.update(exit_code=result.returncode, status="process_succeeded" if result.returncode == 0 else "failed",
                  outputs=[dict(file_record(p, p.parent), absolute_path=str(p)) for p in sorted(of_root.rglob("*")) if p.is_file()])
        of["native_validation"] = validate_orthofinder(config, of, input_records)
        report["status"] = "complete"
    except Exception as error:
        report.update(status="failed", error=str(error))
        raise
    finally:
        save()
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--orthofinder", type=Path, required=True)
    args = parser.parse_args()
    report = run(args.source.resolve(), args.output.resolve(), args.orthofinder.absolute())
    print(json.dumps({"status": report["status"], "accuracy_computed": report["accuracy_computed"]}))


if __name__ == "__main__":
    main()
