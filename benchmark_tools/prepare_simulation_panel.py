"""Materialize the prespecified simulation generation panel without executing it."""

import argparse
import importlib.metadata
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.benchmark_production import file_record
from benchmark_tools.run_zombi_seeded import ZOMBI_COMMIT, family_lengths


SEEDS = tuple(range(20261001, 20261011))
VARIABLE_SEEDS = tuple(range(20261101, 20261111))
CONDITIONS = {"baseline": (2, 1, 0.2), "divergent": (2, 1, 0.8),
              "turnover": (10, 8, 0.2), "divergent_turnover": (10, 8, 0.8)}
DEFAULTS = {"T": "SpeciesTreeParameters.tsv", "G": "GenomeParameters.tsv", "S": "SequenceParameters.tsv"}


def parameter_sets(defaults, seed, condition, allowed_seeds=SEEDS):
    if seed not in allowed_seeds or condition not in CONDITIONS:
        raise ValueError("Seed or condition outside frozen panel")
    duplication, loss, scaling = CONDITIONS[condition]
    overrides = {
        "T": {"SPECIATION": "f:1", "EXTINCTION": "f:0", "STOPPING_RULE": "1",
              "TOTAL_LINEAGES": "8", "SCALE_TREE": "0", "VERBOSE": "0"},
        "G": {"INITIAL_GENOME_SIZE": "100", "MIN_GENOME_SIZE": "1",
              "DUPLICATION": f"f:{duplication}", "LOSS": f"f:{loss}", "TRANSFER": "f:0",
              "ORIGINATION": "f:0", "INVERSION": "f:0", "TRANSPOSITION": "f:0",
              **{event + "_EXTENSION": "g:1" for event in ("DUPLICATION", "TRANSFER", "LOSS", "INVERSION", "TRANSPOSITION")},
              "EVENTS_PER_BRANCH": "1", "PROFILES": "1", "GENE_TREES": "1",
              "RECONCILED_TREES": "1", "SCALE_TREE": "0", "VERBOSE": "0"},
        "S": {"SEQUENCE": "amino-acid", "SEQUENCE_SIZE": "300", "AA_MODEL": "WAG",
              "SCALING": str(scaling), "VERBOSE": "0"},
    }
    return {mode: {**defaults[mode], **overrides[mode], "SEED": str(seed)} for mode in DEFAULTS}


def prepare(source, output, variant="fixed_length_v1"):
    if variant not in {"fixed_length_v1", "variable_length_v2"}:
        raise ValueError("Unknown simulation panel variant")
    variable = variant == "variable_length_v2"
    seeds = VARIABLE_SEEDS if variable else SEEDS
    source, output = source.resolve(), output.resolve()
    if output.exists():
        raise FileExistsError("Refusing to overwrite a materialized panel")
    if any(c.isspace() for c in str(output)):
        raise ValueError("Upstream requires whitespace-free paths")
    commit = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
    if commit != ZOMBI_COMMIT:
        raise ValueError("Wrong simulator source revision")
    subprocess.run(["git", "-C", str(source), "diff", "--exit-code", "HEAD"], check=True, capture_output=True)
    sys.path.insert(0, str(source))
    import AuxiliarFunctions as af
    defaults = {mode: af.read_parameters(str(source / "Parameters" / filename)) for mode, filename in DEFAULTS.items()}
    workflow = Path(__file__).resolve().parent
    scripts = [workflow / name for name in ("prepare_simulation_panel.py", "run_zombi_seeded.py", "zombi_truth.py",
        "derive_simulation_conditions.py", "simulation_conditions.py", "benchmark_production.py")]
    packages = {p: importlib.metadata.version(p) for p in ("Pyvolve", "ete3", "numpy", "scipy", "biopython", "networkx")}
    if packages != {"Pyvolve": "1.1.0", "ete3": "3.1.3", "numpy": "2.2.6", "scipy": "1.15.3", "biopython": "1.86", "networkx": "2.8.8"}:
        raise ValueError("Simulation dependencies differ from validated smoke environment")
    output.mkdir(parents=True)
    report = {"schema_version": 1, "status": "materialized_not_executed", "source_commit": commit,
              "source": str(source), "python": sys.executable, "python_version": sys.version, "packages": packages,
              "protocol": file_record(workflow / "results" / ("PUBLICATION_VARIABLE_LENGTH_PROTOCOL_20260916.md" if variable else
                                                             "PUBLICATION_SIMULATION_PROTOCOL_20260916.md"), workflow),
              "panel_variant": variant, "seeds": list(seeds), "extra_inputs": [],
              "workflow_sources": [dict(file_record(p, workflow), absolute_path=str(p)) for p in scripts],
              "native_sources": [file_record(p, source) for p in sorted(source.glob("*.py"))],
              "native_defaults": [file_record(source / "Parameters" / name, source) for name in DEFAULTS.values()],
              "environment": {"PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"},
              "simulation_runs": [], "datasets": [], "history_equivalence_checks": [],
              "method_inference_launched": False,
              "remaining_gates": ["Generation runner integrity checks", "Seed-level aggregation tests",
                                  "Pinned method launch/conversion manifest", "Controlled resource allocation"]}
    names = sorted({d.metadata["Name"] for d in importlib.metadata.distributions() if d.metadata["Name"]})
    freeze = "".join(f"{name}=={importlib.metadata.version(name)}\n" for name in names)
    # Version inventory only: exclude direct URLs and environment credentials.
    freeze_path = output / "environment.freeze.txt"
    freeze_path.write_text(freeze)
    report["environment_freeze"] = file_record(freeze_path, output)
    for seed in seeds:
        length_path = None
        if variable:
            length_path = output / "family_lengths" / f"{seed}.json"
            length_path.parent.mkdir(exist_ok=True)
            length_path.write_text(json.dumps({"schema_version": 1, "seed": seed,
                "lengths": family_lengths(seed, [str(i) for i in range(1, 101)])}, indent=2, sort_keys=True) + "\n")
            report["extra_inputs"].append(file_record(length_path, output))
        for condition in CONDITIONS:
            label = f"{condition}_{seed}"
            run = output / "native" / label
            prepared = output / "prepared" / label
            params = output / "parameters" / label
            params.mkdir(parents=True)
            parameter_records, commands = {}, []
            for mode, values in parameter_sets(defaults, seed, condition, allowed_seeds=seeds).items():
                path = params / DEFAULTS[mode]
                path.write_text("".join(f"{k}\t{v}\n" for k, v in sorted(values.items())))
                parameter_records[mode] = {"values": values, "file": file_record(path, output)}
                commands.append({"stage": mode, "argv": [sys.executable, str(workflow / "run_zombi_seeded.py"),
                    "--source", str(source), "--mode", mode, "--parameters", str(path), "--output", str(run), "--seed", str(seed)]})
                if mode == "S" and length_path is not None:
                    commands[-1]["argv"].extend(["--family-lengths", str(length_path)])
            commands.append({"stage": "truth", "argv": [sys.executable, str(workflow / "zombi_truth.py"),
                "--run", str(run), "--output", str(prepared)]})
            if length_path is not None:
                commands[-1]["argv"].extend(["--family-lengths", str(length_path)])
            if condition == "baseline":
                derived = output / "derived" / str(seed)
                commands.append({"stage": "derive", "argv": [sys.executable, str(workflow / "derive_simulation_conditions.py"),
                    "--baseline-run", str(run), "--output", str(derived), "--seed", str(seed)]})
                for name in ("missing20", "uneven_taxa", "taxon_count_control"):
                    report["datasets"].append({"condition": name, "seed": seed, "input": str(derived / name / "input"),
                                                "truth": str(derived / name / "truth.json"), "parent": label})
            report["simulation_runs"].append({"label": label, "condition": condition, "seed": seed,
                "parameters": parameter_records, "commands": commands, "native_output": str(run)})
            report["datasets"].append({"condition": condition, "seed": seed,
                                       "input": str(prepared / "input"), "truth": str(prepared / "truth.json"), "parent": label})
        for first, second in (("baseline", "divergent"), ("turnover", "divergent_turnover")):
            report["history_equivalence_checks"].append({"first": f"{first}_{seed}", "second": f"{second}_{seed}",
                "scope": "T and G native biological files; exclude copied parameter files"})
    # Zombi creates its stage directories but not missing parent directories.
    (output / "native").mkdir()
    path = output / "manifest.json"
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--variant", choices=["fixed_length_v1", "variable_length_v2"], default="fixed_length_v1")
    args = parser.parse_args()
    report = prepare(args.source, args.output, args.variant)
    print(f"Materialized {len(report['simulation_runs'])} simulations and {len(report['datasets'])} datasets; no execution")


if __name__ == "__main__":
    main()
