"""Retrospective FastOMA OrthoBench input, tree, task and output identity audit."""

import argparse
from collections import Counter
import csv
import json
from pathlib import Path

from benchmark_tools.admit_qfo_corrected_fastoma import tree_clades
from benchmark_tools.audit_ob_orthofinder_provenance import fasta
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


def match_input(genes, sources):
    matches = [name for name, value in sources.items() if genes.keys() == value.keys()]
    if len(matches) != 1:
        raise ValueError("No unique frozen species gene-set match")
    name = matches[0]
    return name, sorted(g for g in genes if genes[g] != sources[name][g])


def table_counts(path, universe, root_hogs=False):
    header = ["RootHOG", "Protein", "OMAmerRootHOG"] if root_hogs else ["Group", "Protein"]
    groups, seen = set(), set()
    with path.open(newline="") as stream:
        rows = csv.reader(stream, delimiter="\t")
        if next(rows) != header:
            raise ValueError("Unexpected FastOMA table header")
        for row in rows:
            if len(row) != len(header) or not row[0] or row[1] not in universe or row[1] in seen:
                raise ValueError("Invalid or repeated native member")
            groups.add(row[0])
            seen.add(row[1])
    if not groups:
        raise ValueError("Empty FastOMA table")
    return dict(groups=len(groups), assigned_genes=len(seen))


def task_records(work):
    rows = []
    for path in sorted(work.glob("*/*/.exitcode")):
        value = path.read_text().strip()
        if not value.isdigit():
            raise ValueError("Malformed task exit code")
        directory = path.parent
        rows.append(dict(directory=str(directory), exit_code=int(value),
                         records=[record(p) for p in (path, directory / ".command.sh", directory / ".command.run", directory / ".command.log")]))
    if not rows:
        raise ValueError("No retained task exit records")
    return rows


def audit(root):
    reports = root / "benchmark_tools/results"
    inventory_path, admission_path = [reports / name for name in
        ("installed_ob_input_inventory_20260926.json", "retained_ob_comparator_readback_20260926.json")]
    checked = [record(inventory_path), record(admission_path)]
    if [r["sha256"] for r in checked] != ["8b429398e6c381fcdec8289c641a7f20969edcad5595e3d3f6930335601571bd", "5dcb4e65f277eb9debf2bbacb1d4e298c678ac13fc981e75c4ba78e4851eeb55"]:
        raise ValueError("Changed frozen input/admission reports")
    inventory, admission = json.loads(inventory_path.read_text()), json.loads(admission_path.read_text())
    sources, source_records = {}, {}
    for item in inventory["inputs"]:
        check(item)
        checked.append(item)
        name = Path(item["path"]).name
        sources[name], source_records[name] = fasta(Path(item["path"])), item
    base = root / "benchmarks/results/fastoma_run"
    rows, universe = [], set()
    for path in sorted((base / "in_folder/proteome").glob("*.fa")):
        item = record(path)
        checked.append(item)
        genes = fasta(path)
        name, mismatches = match_input(genes, sources)
        if universe & genes.keys():
            raise ValueError("Duplicate staged gene IDs")
        universe.update(genes)
        rows.append(dict(species=path.stem, staged=item, frozen_input=source_records[name], genes=len(genes), sequence_mismatches=mismatches))
    if len(rows) != 12 or len({r["frozen_input"]["path"] for r in rows}) != 12 or len(universe) != inventory["genes"]:
        raise ValueError("Incomplete staged universe")
    species = [r["species"] for r in rows]
    supplied, checked_tree = base / "in_folder/species_tree.nwk", base / "output/species_tree_checked.nwk"
    checked.extend(record(p) for p in (supplied, checked_tree))
    topology_same = tree_clades(supplied, species) == tree_clades(checked_tree, species)
    prediction = next(r for r in admission["rows"] if r["key"] == "fastoma_0_3_5")["prediction"]
    check(prediction)
    checked.append(prediction)
    published = record(base / "output/OrthologousGroups.tsv")
    checked.append(published)
    if any(prediction[k] != published[k] for k in ("sha256", "bytes")):
        raise ValueError("Published final groups differ from scored native groups")
    final_counts = table_counts(Path(prediction["path"]), universe)
    root_table = base / "output/RootHOGs.tsv"
    checked.append(record(root_table))
    root_counts = table_counts(root_table, universe, True)
    tasks = task_records(base / "work")
    checked.extend(r for task in tasks for r in task["records"])
    collector = next(t for t in tasks if t["directory"] == str(Path(prediction["path"]).parent))
    if collector["exit_code"] != 0:
        raise ValueError("Scored collection task did not exit successfully")
    log_path = base.with_name("fastoma_run.log")
    checked.append(record(log_path))
    log = log_path.read_text()
    relevant = [line for line in log.splitlines() if line.strip().startswith(
        ("Cmd line:", "Manifest's pipeline version:", "Completed at", "Duration", "CPU hours", "Succeeded", "Failed", "Processes"))]
    if "Manifest's pipeline version: 0.3.5" not in log or "Succeeded   : 90" not in log or "Failed      : 10" not in log:
        raise ValueError("Unexpected retained workflow summary")
    counts = dict(Counter(str(t["exit_code"]) for t in tasks))
    if counts != {"0": 90, "137": 10}:
        raise ValueError("Native task counts disagree with reviewed run")
    checked.extend(record(Path(__file__).with_name(name)) for name in
                   ("audit_ob_fastoma_provenance.py", "audit_ob_orthofinder_provenance.py", "admit_qfo_corrected_fastoma.py"))
    for item in checked:
        check(item)
    return dict(status="retained_fastoma_ob_provenance_readback", checked_records=checked, inputs=rows,
                all_protein_sequences_match=all(not r["sequence_mismatches"] for r in rows),
                supplied_tree_topology_preserved=topology_same, final_groups=final_counts, root_hogs_diagnostic=root_counts,
                scored_prediction=prediction, published_prediction=published, task_exit_counts=counts, tasks=tasks,
                retained_run_summary=relevant, cpu_time_precision="rounded Nextflow report", peak_rss_kib=None,
                publication_ready=False, historical_consumption_proven=False, inference_rerun=False,
                limitations=["Supplied species tree, not independent FastOMA tree inference; tree origin is not established by this audit.",
                             "Ninety successful and ten exit-137 tasks are retained; exit 137 alone does not establish out-of-memory as the cause.",
                             "Workflow-reported elapsed/CPU time is descriptive, not independently verified controlled runtime or memory evidence.",
                             "Staged sequence and output consistency does not prove historical binary, database or input-consumption identity.",
                             "No inference/scoring rerun or substitution of root HOGs for final groups."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.root.resolve())
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    raise SystemExit(0 if result["all_protein_sequences_match"] and result["supplied_tree_topology_preserved"] else 1)
