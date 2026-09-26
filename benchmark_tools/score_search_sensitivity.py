"""Independently validate native search hits and score the frozen simulation panel."""

import argparse
from collections import Counter
import csv
from fractions import Fraction
import json
import math
from pathlib import Path
import subprocess

from Bio import SeqIO

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_search_sensitivity_cell import INSTALL_SHA, PLAN_SHA, diamond_commands


def reference(row):
    check(row["truth"])
    truth = json.loads(Path(row["truth"]["path"]).read_text())
    genes = {}
    for item in row["inputs"]:
        check(item)
        for seq in SeqIO.parse(item["path"], "fasta"):
            if seq.id in genes:
                raise ValueError("Duplicate reference gene")
            genes[seq.id] = (Path(item["path"]).name, len(seq))
    families, denominators = {}, {}
    for family, members in truth["families"].items():
        for gene in members:
            if gene not in genes or gene in families:
                raise ValueError("Invalid reference family membership")
            families[gene] = family
        counts = Counter(genes[g][0] for g in members)
        denominators[family] = sum(n * (len(members) - n) for n in counts.values())
    if set(families) != set(genes) or denominators != row["family_directed_homology_pairs"]:
        raise ValueError("Reference denominators differ from frozen plan")
    return genes, families, denominators


def read_hits(path, genes, engine, target=None):
    hits, seen = [], set()
    with path.open() as stream:
        reader = csv.reader(stream, delimiter="\t")
        if engine == "hmm" and next(reader, None) != ["query_species", "target_species", "query_id", "target_id", "score", "evalue"]:
            raise ValueError("Invalid HMM header")
        for row in reader:
            if len(row) != (6 if engine == "hmm" else 7):
                raise ValueError("Invalid hit row width")
            if engine == "hmm":
                qs, ts, q, t, score, e = row
                if q not in genes or t not in genes or (qs, ts) != (genes[q][0], genes[t][0]):
                    raise ValueError("Foreign gene/species in HMM hits")
                numeric = [float(score), float(e)]
                if not 0 <= numeric[-1] < 1e-4:
                    raise ValueError("HMM E-value outside native cutoff")
            else:
                q, t, qlen, tlen, score, bits, e = row
                if q not in genes or t not in genes or genes[t][0] != target:
                    raise ValueError("Foreign gene/target in DIAMOND hits")
                if (int(qlen), int(tlen)) != (genes[q][1], genes[t][1]):
                    raise ValueError("DIAMOND sequence length mismatch")
                numeric = [float(score), float(bits), float(e)]
                if not 0 <= numeric[-1] <= 1:
                    raise ValueError("DIAMOND E-value outside native cutoff")
            if not all(map(math.isfinite, numeric)) or (q, t) in seen:
                raise ValueError("Nonfinite or duplicate directed hit")
            seen.add((q, t))
            hits.append((q, t, numeric[-1]))
    return hits


def score(hits, genes, families, denominators, cutoff):
    recovered = Counter()
    nonhomologs, total = 0, 0
    for q, t, evalue in hits:
        if evalue > cutoff or genes[q][0] == genes[t][0]:
            continue
        total += 1
        if families[q] == families[t]:
            recovered[families[q]] += 1
        else:
            nonhomologs += 1
    if any(recovered[f] > denominators[f] for f in recovered):
        raise ValueError("Recall numerator exceeds denominator")
    eligible = [f for f in denominators if denominators[f]]
    if not eligible:
        raise ValueError("No eligible homology family")
    recall = sum((Fraction(recovered[f], denominators[f]) for f in eligible), Fraction()) / len(eligible)
    return dict(primary_recall=float(recall), exact_recall=str(recall), total_hits=total,
                homolog_hits=sum(recovered.values()), nonhomolog_hits=nonhomologs,
                homology_denominator=sum(denominators.values()), eligible_families=len(eligible))


def mean(rows, engine, index=None):
    values = [r[engine] if index is None else r[engine][index] for r in rows]
    return sum((Fraction(v["exact_recall"]) for v in values), Fraction()) / len(values)


def select_cutoff(rows, grid):
    calibration = [r for r in rows if r["split"] == "calibration"]
    if not calibration:
        raise ValueError("Missing calibration rows")
    baseline = mean(calibration, "hmm")
    return min(range(len(grid)), key=lambda i: (abs(mean(calibration, "diamond", i) - baseline), grid[i]))


def aggregate(rows, index):
    hmm, diamond = mean(rows, "hmm"), mean(rows, "diamond", index)
    return dict(datasets=len(rows), hmm_recall=float(hmm), diamond_recall=float(diamond),
                difference=float(diamond - hmm), exact_difference=str(diamond - hmm),
                hmm_hits=sum(r["hmm"]["total_hits"] for r in rows),
                diamond_hits=sum(r["diamond"][index]["total_hits"] for r in rows),
                hmm_nonhomolog_hits=sum(r["hmm"]["nonhomolog_hits"] for r in rows),
                diamond_nonhomolog_hits=sum(r["diamond"][index]["nonhomolog_hits"] for r in rows),
                hmm_homolog_hits=sum(r["hmm"]["homolog_hits"] for r in rows),
                diamond_homolog_hits=sum(r["diamond"][index]["homolog_hits"] for r in rows),
                homology_denominator=sum(r["hmm"]["homology_denominator"] for r in rows))


def validate_cell(path, row, index, submission, plan):
    receipt = json.loads((path / "execution.json").read_text())
    if (receipt["status"] != "native_completed_pending_independent_readback"
            or receipt["dataset_index"] != index or receipt["attempt"] != 1
            or receipt["slurm_array_job_id"] != submission["job_id"]
            or receipt["slurm_array_task_id"] != str(index)
            or any(receipt[k] != row[k] for k in ("condition", "seed", "split"))
            or receipt["plan"] != submission["plan"]):
        raise ValueError("Incomplete or inconsistent cell receipt")
    for item in receipt["checked_records"] + receipt["private_inputs"] + receipt["outputs"] + receipt["logs"]:
        check(item)
    install_record = next(r for r in plan["checked_evidence"] if r["sha256"] == INSTALL_SHA)
    check(install_record)
    installed = json.loads(Path(install_record["path"]).read_text())
    required_checks = [submission["plan"], plan["diamond"], install_record, *installed["checked_records"],
                       *[s for s in submission["sources"] if Path(s["path"]).name in
                         {"run_search_sensitivity_cell.py", "search_sensitivity_hmm.py"}]]
    if any(r not in receipt["checked_records"] for r in required_checks):
        raise ValueError("Missing pinned runtime/source checks")
    staging = next(Path(r["path"]) for r in installed["checked_records"] if Path(r["path"]).name == "staging.json")
    if receipt["python_command_path"] != str(staging.parent / "venv_clean/bin/python"):
        raise ValueError("Unexpected interpreter")
    expected_env = dict(PATH="/usr/bin:/bin", HOME=str(path), LC_ALL="C", OMP_NUM_THREADS="1",
                        OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", PYTHONHASHSEED="0")
    if receipt["environment"] != expected_env:
        raise ValueError("Unexpected execution environment")
    expected_inputs = {(Path(r["path"]).name, r["bytes"], r["sha256"]) for r in row["inputs"]}
    actual_inputs = { (p.name, (r := record(p))["bytes"], r["sha256"]) for p in (path / "input").iterdir() }
    if actual_inputs != expected_inputs:
        raise ValueError("Private input inventory mismatch")
    helper = Path(submission["executor"]) / "benchmark_tools/search_sensitivity_hmm.py"
    commands = [[receipt["python_command_path"], "-I", str(helper), "--input", str(path / "input"), "--output", str(path / "hmm")]]
    for target in sorted((path / "input").iterdir()):
        commands.extend(diamond_commands(plan["diamond"]["path"], path / "queries.fasta", target, path / "diamond"))
    if [s["command"] for s in receipt["stages"]] != commands or any(s["returncode"] != 0 for s in receipt["stages"]):
        raise ValueError("Incomplete or altered stage commands")
    hmm = json.loads((path / "hmm/receipt.json").read_text())
    if (hmm["source"] != record(helper) or hmm["engine"] not in receipt["checked_records"]
            or hmm["settings"] != plan["hmm_settings"] or hmm["truth_labels_loaded"]
            or not hmm["isolated"] or hmm["inputs"] != receipt["private_inputs"]):
        raise ValueError("HMM identity/settings mismatch")
    output_records = {r["path"]: r for r in receipt["outputs"]}
    required = [path / "hmm/hits.tsv", path / "hmm/receipt.json", path / "queries.fasta"]
    for target in sorted((path / "input").iterdir()):
        required.extend([path / "diamond" / (target.stem + suffix) for suffix in (".tsv", ".dmnd")])
    if set(output_records) != {str(p) for p in required} or hmm["hits"] != output_records[str(path / "hmm/hits.tsv")]:
        raise ValueError("Native output inventory mismatch")
    return receipt, hmm


def run(submission_path, output):
    if output.exists():
        raise FileExistsError(output)
    submission = json.loads(submission_path.read_text())
    check(submission["plan"])
    if submission["plan"]["sha256"] != PLAN_SHA:
        raise ValueError("Unrecognized plan")
    for source in submission["sources"]:
        check(source)
    plan = json.loads(Path(submission["plan"]["path"]).read_text())
    accounting = subprocess.check_output(["sacct", "-j", submission["job_id"], "--format=JobID,State,ExitCode", "-n", "-P"], text=True)
    states = {r[0]: r[1:3] for line in accounting.splitlines() if len(r := line.split("|")) >= 3}
    if any(states.get(submission["job_id"] + "_" + str(i)) != ["COMPLETED", "0:0"] for i in range(70)):
        raise ValueError("All 70 native cells must be terminal and successful")
    results, receipts = [], []
    root = Path(submission["output"])
    for index, row in enumerate(plan["datasets"]):
        path = root / f"cell_{index}"
        execution, native = validate_cell(path, row, index, submission, plan)
        genes, families, denominators = reference(row)
        hmm = read_hits(path / "hmm/hits.tsv", genes, "hmm")
        counts = Counter((genes[q][0], genes[t][0]) for q, t, _ in hmm)
        expected_pairs = {(a, b) for a in {g[0] for g in genes.values()} for b in {g[0] for g in genes.values()}}
        pair_rows = {(r["query"], r["target"]): r for r in native["pairs"]}
        if len(pair_rows) != len(native["pairs"]) or set(pair_rows) != expected_pairs:
            raise ValueError("HMM species-pair inventory mismatch")
        if any(r["hits"] != counts[p] or r["candidates"] < r["hits"] for p, r in pair_rows.items()):
            raise ValueError("HMM count readback mismatch")
        diamond = []
        for item in row["inputs"]:
            target = Path(item["path"])
            diamond.extend(read_hits(path / "diamond" / (target.stem + ".tsv"), genes, "diamond", target.name))
        results.append(dict(condition=row["condition"], seed=row["seed"], split=row["split"],
                            hmm=score(hmm, genes, families, denominators, 1e-4),
                            diamond=[score(diamond, genes, families, denominators, e) for e in plan["evalue_grid"]]))
        receipts.append(record(path / "execution.json"))
    chosen = select_cutoff(results, plan["evalue_grid"])
    reporting = [r for r in results if r["split"] == "reporting"]
    overall = aggregate(reporting, chosen)
    conditions = {c: aggregate([r for r in reporting if r["condition"] == c], chosen) for c in sorted({r["condition"] for r in reporting})}
    passed = (abs(Fraction(overall["exact_difference"])) <= Fraction(2, 100)
              and all(abs(Fraction(r["exact_difference"])) <= Fraction(5, 100) for r in conditions.values()))
    result = dict(status="complete_panel_scored", source=record(__file__), submission=record(submission_path),
                  plan=submission["plan"], receipts=receipts, accounting=accounting,
                  evalue_grid=plan["evalue_grid"], datasets=results,
                  calibration_grid=[aggregate([r for r in results if r["split"] == "calibration"], i) for i in range(len(plan["evalue_grid"]))],
                  selected_cutoff=plan["evalue_grid"][chosen], reporting_overall=overall,
                  reporting_conditions=conditions, reporting_match_gate_passed=passed,
                  diamond_cutoff_convention="reported E-value <= cutoff; HMM native strict < 1e-4",
                  independent_validation=False, scientific_defaults_changed=False,
                  computational_effort_matched=False, publication_ready=False)
    for item in receipts:
        check(item)
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--submission", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.submission.resolve(), args.output.absolute())
