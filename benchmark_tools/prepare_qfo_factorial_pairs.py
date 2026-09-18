"""Convert one frozen QfO cell without conflating native pairs and group cliques."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import SeqIO
from benchmark_tools.admit_qfo_factorial_cell import PREPARED_SHA, admit
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_qfo_factorial_cell import select_cell, verify_prepared
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.simulation_method_outputs import orthohmm_pairs
from benchmark_tools.qfo_filter_pairs import load_mapping, filter_pairs
from qfo_benchmark.og_to_pairwise import _strip_to_uniprot

RECOVERED_SHA = "ce6f19cd005b886a91dc3aee63cb7b41e7f448d8cad5e65ff7cb10d92a636100"


def select_conversion(prepared, index):
    if type(index) is not int or not 0 <= index < 8:
        raise ValueError("Require factorial cell index0-7")
    reconciliation, output, executor = select_cell(prepared, index // 2)
    cell = prepared["cells"][index]
    if cell["reconciliation"] is not bool(index % 2):
        raise ValueError("Unexpected factorial cell order")
    return cell, reconciliation, output, executor


def reusable_stage(cell, arm, recovered, fastas, mapping):
    if cell["reconciliation"] or cell["candidate_expansion"]:
        raise ValueError("Only unexpanded R-off cells can reuse recovered pairs")
    if (recovered["status"] != "four_stage_pairs_prepared_unscored" or recovered.get("accuracy_evaluated") is not False
            or len(recovered["stages"]) != 4 or recovered["input_fastas"] != fastas or recovered["mapping"] != mapping):
        raise ValueError("Recovered conversion provenance differs")
    index = 3 if cell["profile_expansion"] else 1
    stage = recovered["stages"][index]
    expected_label = "strict_profiles_refined" if cell["profile_expansion"] else "multipass_refined"
    if stage["index"] != index or stage["stage"] != expected_label:
        raise ValueError("Recovered stage identity differs")
    a, b = arm["candidate_partition"], stage["partition"]
    if (a["sha256"], a["bytes"]) != (b["sha256"], b["bytes"]):
        raise ValueError("Candidate partition differs from recovered conversion")
    if not 0 < stage["retained_pairs"] <= stage["total_pairs"] or stage["removed_mapping_pairs"] != stage["total_pairs"] - stage["retained_pairs"]:
        raise ValueError("Invalid recovered conversion counts")
    return stage


def write_native_pairs(native, output, owners, expected_count):
    # Require injective accession normalization before stripping composite IDs.
    normalized = {}
    for gene in owners:
        accession = _strip_to_uniprot(gene)
        if accession in normalized:
            raise ValueError("Ambiguous normalized accession")
        normalized[accession] = gene
    count, previous = 0, None
    with output.open("x") as stream:
        for pair in orthohmm_pairs(native, owners):
            if pair[0] >= pair[1] or (previous is not None and pair <= previous):
                raise ValueError("Native pair order or uniqueness changed")
            previous = pair
            a, b = sorted(_strip_to_uniprot(g) for g in pair)
            stream.write(f"{a}\t{b}\n")
            count += 1
    if count != expected_count:
        raise ValueError("Native conversion count differs from admitted count")
    return count


def prepare(root, index):
    if Path.cwd().resolve() != root:
        raise ValueError("Run from original repository verification directory")
    path = root / "benchmark_tools/results/qfo_factorial_prepared_20260917.json"
    prepared = read_frozen(path, PREPARED_SHA)
    cell, reconciliation, _, _ = select_conversion(prepared, index)
    scheduler = verify_prepared(prepared, reconciliation, 21670)
    output = root / "benchmarks/results/qfo_factorial_pairs_v1" / cell["label"]
    if output.exists():
        raise FileExistsError(output)
    recovered_path = root / "benchmark_tools/results/qfo_recovered_stage_pairs_20260917.json"
    recovered = read_frozen(recovered_path, RECOVERED_SHA)
    mapping = recovered["mapping"]
    check(mapping)
    fastas = prepared["input_fastas"]
    fasta = Path(fastas[0]["path"]).parent
    arm = prepared["candidate_arms"][f"p{int(cell['profile_expansion'])}_c{int(cell['candidate_expansion'])}"]
    converter = root / "qfo_benchmark/og_to_pairwise.py"
    source = record(__file__)
    output.mkdir(parents=True)
    report = {"status": "preparing", "accuracy_evaluated": False, "index": index, "cell": cell["label"],
              "participant": f"ohmm_qfo_factorial_{cell['label']}", "job_id": os.environ.get("SLURM_JOB_ID"),
              "array_job_id": os.environ.get("SLURM_ARRAY_JOB_ID"), "array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"),
              "prepared": record(path), "preparation_scheduler": scheduler, "mapping": mapping, "source": source,
              "input_fastas": fastas, "candidate_partition": arm["candidate_partition"],
              "semantics": "native phylogenetically inferred pairs" if cell["reconciliation"] else "cross-species group-derived clique pairs",
              "sources": [record(converter), record(sys.modules["qfo_benchmark.og_to_pairwise"].__file__), record(Path(__file__).with_name("qfo_filter_pairs.py")),
                          record(Path(__file__).with_name("simulation_method_outputs.py")),
                          record(Path(__file__).with_name("admit_qfo_factorial_cell.py"))]}
    checked = [source, report["prepared"], mapping, *fastas, arm["candidate_partition"], *report["sources"]]
    try:
        if not cell["reconciliation"] and not cell["candidate_expansion"]:
            stage = reusable_stage(cell, arm, recovered, fastas, mapping)
            checked.extend([record(recovered_path), *recovered["sources"], stage["partition"],
                            stage["pairs"], stage["filtered_pairs"], stage["conversion_log"]])
            for item in checked:
                check(item)
            report.update(reused_conversion={"manifest": record(recovered_path), "stage": stage["stage"]},
                          **{k: stage[k] for k in ("pairs", "filtered_pairs", "total_pairs", "retained_pairs", "removed_mapping_pairs")})
        else:
            pairs, filtered = output / "pairs.tsv", output / "pairs.qfo.tsv"
            if cell["reconciliation"]:
                # Independently revalidate the terminal native run; do not trust
                # a stale admission flag or a RootHOG path in the group plan.
                native = admit(root, index // 2)
                native_path = output / "native_admission_recheck.json"
                native_path.write_text(json.dumps(native, indent=2, sort_keys=True) + "\n")
                report["native_admission"] = record(native_path)
                report["native_input"] = native["native_pairs"]
                checked.extend([report["native_admission"], native["native_pairs"]])
                native_manifest = json.loads(Path(native["native_group_integrity"]["native_manifest"]["path"]).read_text())
                taxa = {r["filename"]: r["taxon"] for r in native_manifest["input_proteomes"]}
                owners = {}
                for item in fastas:
                    for sequence in SeqIO.parse(item["path"], "fasta"):
                        if sequence.id in owners:
                            raise ValueError("Duplicate input gene")
                        owners[sequence.id] = taxa[Path(item["path"]).name]
                write_native_pairs(Path(native["native_pairs"]["path"]), pairs, owners, native["native_pair_count"])
            else:
                command = [sys.executable, str(converter), arm["candidate_partition"]["path"], str(fasta)]
                with pairs.open("x") as stream, (output / "conversion.log").open("x") as log:
                    subprocess.run(command, stdout=stream, stderr=log, check=True)
                report["command"] = command
                report["conversion_log"] = record(output / "conversion.log")
            total, retained = filter_pairs(pairs, filtered, load_mapping(Path(mapping["path"])))
            if not 0 < retained <= total:
                raise ValueError("No valid converted prediction pairs")
            report.update(pairs=record(pairs), filtered_pairs=record(filtered), total_pairs=total,
                          retained_pairs=retained, removed_mapping_pairs=total - retained)
        for item in checked:
            check(item)
        report["checked_records"] = checked
        report["status"] = "cell_pairs_prepared_unscored"
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        (output / "results.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--index", type=int, choices=range(8), required=True)
    args = parser.parse_args()
    prepare(args.root.resolve(), args.index)
