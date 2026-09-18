"""Stage and independently audit corrected OrthoMCL inputs; do not run BLAST."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.prepare_qfo_corrected_proteinortho import PRIMARY_SHA, REGISTRY_SHA
from benchmark_tools.run_qfo_corrected_sonic import verify_corrected_inputs
from benchmark_tools.prepare_orthomcl_inputs import prepare_inputs
from benchmark_tools.audit_orthomcl_input_parity import compare_inputs

SOFTWARE = Path("/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE")


def commands(work, blastall, formatdb, threads=180):
    if type(threads) is not int or threads < 1:
        raise ValueError("Positive integer thread count required")
    fasta = str(Path(work) / "all.fa")
    return {"formatdb": [str(formatdb), "-i", fasta, "-p", "t"],
            "blast": [str(blastall), "-p", "blastp", "-i", fasta, "-d", fasta,
                      "-e", "1e-5", "-o", str(Path(work) / "all.blast.partial"),
                      "-m", "8", "-a", str(threads), "-v", "1000", "-b", "1000"]}


def require_parity(result):
    if (result["proteomes"] != 78 or result["total_sequences"] != 984137
            or result["identical_sequences"] != 984137 or result["sequence_difference_count"] != 0
            or result["differences"] or result["genome_map_complete_and_correct"] is not True):
        raise ValueError("Corrected OrthoMCL input parity failed")


def prepare(root, output, destination):
    if output.exists() or destination.exists():
        raise FileExistsError("Require fresh OrthoMCL output and manifest")
    results = root / "benchmark_tools/results"
    primary_path = results / "qfo_corrected_primary_commands_20260918.json"
    primary = read_frozen(primary_path, PRIMARY_SHA)
    directory = Path(primary["input_directory"])
    inputs = json.loads((directory / "staging_manifest.json").read_text())["input_fastas"]
    verify_corrected_inputs(primary, inputs)
    registry = results / "publication_comparison_orthomcl_complete_20260916.json"
    read_frozen(registry, REGISTRY_SHA)
    tools = {"blastall": SOFTWARE / "blast-2.2.13/bin/blastall",
             "formatdb": SOFTWARE / "blast-2.2.13/bin/formatdb",
             "mcl": SOFTWARE / "mcl-02-063/bin/mcl",
             "perl": SOFTWARE / "ORTHOMCLV1.4/venv_orthomcl/bin/perl",
             "orthomcl": SOFTWARE / "ORTHOMCLV1.4/orthomcl.pl",
             "module": SOFTWARE / "ORTHOMCLV1.4/orthomcl_module.pm"}
    sources = [record(root / "benchmark_tools" / name) for name in (
        "prepare_orthomcl_inputs.py", "audit_orthomcl_input_parity.py", "audit_qfo_input_sequences.py",
        "configure_orthomcl_1_4.py", "parallelize_orthomcl_pairs.py", "convert_orthomcl_blast.py",
        "run_orthomcl_qfo.slurm", "run_orthomcl_final_groups_qfo.slurm", "audit_orthomcl_blast.py",
        "run_qfo_corrected_sonic.py")]
    checked = [record(primary_path), record(registry), *sources, *[record(p) for p in tools.values()], *inputs]
    output.mkdir(parents=True, exist_ok=False)
    work = output / "work"
    work.mkdir()
    failure = output / "preparation_status.json"
    try:
        summary = prepare_inputs(directory, work / "all.fa", work / "all.gg", output / "input_summary.json")
        if summary["taxon_count"] != 78 or summary["sequence_count"] != 984137:
            raise ValueError("Unexpected prepared dimensions")
        parity = compare_inputs([Path(item["path"]) for item in inputs], work / "all.fa", work / "all.gg")
        require_parity(parity)
        for item in checked:
            check(item)
        verify_corrected_inputs(primary, inputs)
        report = {"status": "corrected_orthomcl_inputs_prepared_pending_execution_freeze",
                  "execution_authorized": False, "accuracy_admitted": False,
                  "source": record(__file__), "checked_records": checked, "input_fastas": inputs,
                  "primary_manifest": record(primary_path), "output_root": str(output),
                  "prepared_inputs": [record(work / "all.fa"), record(work / "all.gg"), record(output / "input_summary.json")],
                  "summary": summary, "independent_parity": parity,
                  "tools": {name: record(path) for name, path in tools.items()},
                  "search_commands": commands(work, tools["blastall"], tools["formatdb"]),
                  "search_reuse": False,
                  "resources_proposed": {"cpus": 180, "pair_workers": 64, "memory_gib": 900,
                                         "time_limit_days": 14, "node": "bizon"},
                  "pair_semantics": "Final native OrthoMCL groups expanded to cross-species pairs, not pre-clustering matrix edges.",
                  "remaining_gates": [
                      "Freeze configured native modules, runtime environment and dependencies, scheduler allocation and guarded launcher.",
                      "Fresh formatdb and legacy BLAST searches against the full corrected database; no old-release search reuse.",
                      "Native BLAST log/query-failure audit, database and BPO/index validation before clustering admission.",
                      "Preserve legacy masking/scoring defaults and failed queries; do not silently replace the search engine.",
                      "Freeze final-group conversion and corrected scorer participant before scoring; never invoke old matrix-edge scoring.",
                  ],
                  "limitations": ["Prepared sequence/genome-map parity does not verify an unbuilt BLAST database.",
                                  "Proposed shared-host resources are not matched timing evidence.",
                                  "Original-release failed-query counts cannot be transferred to this corrected run."]}
        with destination.open("x") as stream:
            json.dump(report, stream, indent=2, sort_keys=True)
            stream.write("\n")
        failure.write_text(json.dumps({"status": report["status"], "manifest": record(destination)}, indent=2) + "\n")
        return report
    except Exception as exc:
        failure.write_text(json.dumps({"status": "preparation_failed", "error": str(exc)}, indent=2) + "\n")
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "output-root", "manifest"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    prepare(args.root.resolve(), args.output_root.resolve(), args.manifest.resolve())
