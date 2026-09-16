"""Small end-to-end reproducibility smoke test, not an accuracy benchmark."""

import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.benchmark_production import file_record
from benchmark_tools.run_zombi_seeded import ZOMBI_COMMIT


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source, output = args.source.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    sys.path.insert(0, str(source))
    import AuxiliarFunctions as af
    driver = Path(__file__).with_name("run_zombi_seeded.py").resolve()
    report = {"schema_version": 1, "status": "running", "source_commit": ZOMBI_COMMIT,
              "accuracy_computed": False, "truth_validated": False,
              "python": sys.version, "packages": {p: importlib.metadata.version(p) for p in
                  ("Pyvolve", "ete3", "numpy", "scipy", "biopython", "networkx")},
              "scripts": [file_record(p, p.parent) for p in (driver, Path(__file__).resolve())],
              "runs": {}}
    report_path = output / "report.json"
    settings = {
        "T": ("SpeciesTreeParameters.tsv", {"TOTAL_LINEAGES": "4", "EXTINCTION": "f:0", "VERBOSE": "0"}),
        "G": ("GenomeParameters.tsv", {"INITIAL_GENOME_SIZE": "10", "MIN_GENOME_SIZE": "1",
            "DUPLICATION": "f:0.2", "LOSS": "f:0.1", "TRANSFER": "f:0", "ORIGINATION": "f:0",
            "INVERSION": "f:0", "TRANSPOSITION": "f:0", "VERBOSE": "0", "RECONCILED_TREES": "1"}),
        "S": ("SequenceParameters.tsv", {"SEQUENCE": "amino-acid", "SEQUENCE_SIZE": "100",
            "AA_MODEL": "WAG", "SCALING": "0.2", "VERBOSE": "0"}),
    }
    env = os.environ.copy()
    env.update(PYTHONHASHSEED="0", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    try:
        for label, seed in (("repeat_a", 20260916), ("repeat_b", 20260916), ("independent", 20260917)):
            run = output / label
            params = output / (label + "_parameters")
            params.mkdir()
            record = report["runs"][label] = {"seed": seed, "commands": [], "stages": {}}
            for mode, (filename, overrides) in settings.items():
                native = source / "Parameters" / filename
                parameters = af.read_parameters(str(native))
                parameters.update(overrides, SEED=str(seed))
                path = params / filename
                path.write_text("".join(f"{k}\t{v}\n" for k, v in sorted(parameters.items())))
                command = [sys.executable, str(driver), "--source", str(source), "--mode", mode,
                           "--parameters", str(path), "--output", str(run), "--seed", str(seed)]
                record["commands"].append(command)
                report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
                with (output / f"{label}_{mode}.log").open("w") as log:
                    subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, env=env, check=True)
                record["stages"][mode] = {"parameters": file_record(path, output),
                    "default_parameters": file_record(native, source),
                    "outputs": [file_record(p, run) for p in sorted((run / mode).rglob("*")) if p.is_file()]}
            record["sequence_files"] = [file_record(p, run) for p in sorted((run / "S").glob("*.fasta"))]
            if not record["sequence_files"]:
                raise ValueError("No simulated sequence files")
        def products(label):
            return {r["path"]: r["sha256"] for stage in report["runs"][label]["stages"].values() for r in stage["outputs"]}
        same = products("repeat_a") == products("repeat_b")
        different = report["runs"]["repeat_a"]["sequence_files"] != report["runs"]["independent"]["sequence_files"]
        report.update(same_seed_all_products_identical=same, different_seed_sequences_differ=different)
        if not same or not different:
            raise ValueError("Reproducibility smoke test failed")
        report["status"] = "complete"
    except Exception as error:
        report.update(status="failed", error=str(error))
        raise
    finally:
        report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
