"""Prepare exact-byte, prespecified legacy BLAST diagnostic queries; do not run."""

import argparse
import json
from pathlib import Path
import sys

import Bio
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.run_qfo_corrected_blast import verify, PLAN_SHA

FIRST = "sp|H2VFI5|SIAA_NEIMB"
FAILED = "sp|P62945|RL41_HUMAN"
SELENO = "sp|P63302|SELW_HUMAN"
BOUNDARY = "tr|F1MRE1|F1MRE1_BOVIN"


def panel_ids(ids):
    if not ids or ids[0] != FIRST or len(ids) != len(set(ids)):
        raise ValueError("Unexpected first query or duplicate input IDs")
    if not {FIRST, FAILED, SELENO, BOUNDARY}.issubset(ids):
        raise ValueError("Required diagnostic query missing")
    boundary = ids.index(BOUNDARY)
    if boundary == 0:
        raise ValueError("Boundary has no predecessor")
    chosen = {FIRST, FAILED, SELENO, BOUNDARY, ids[boundary - 1]}
    if len(chosen) != 5:
        raise ValueError("Diagnostic roles overlap")
    return [gene for gene in ids if gene in chosen]


def prepare(plan_path, runtime_path, output):
    if output.exists():
        raise FileExistsError(output)
    plan = verify(plan_path, runtime_path)
    source = next(item for item in plan["prepared_inputs"] if item["path"].endswith("/all.fa"))
    status = json.loads((Path(plan["output_root"]) / "search_execution/status.json").read_text())
    database = status["database_files"]
    for item in database:
        check(item)
    index = SeqIO.index(source["path"], "fasta")
    try:
        ids = list(index)
        selected = panel_ids(ids)
        if len(ids) != 984137:
            raise ValueError("Unexpected corrected input size")
        raw = {gene: index.get_raw(gene) for gene in selected}
    finally:
        index.close()
    check(source)
    output.mkdir(parents=True, exist_ok=False)
    queries = []
    for number, gene in enumerate(selected):
        path = output / f"query_{number:02d}.fa"
        with path.open("xb") as stream:
            stream.write(raw[gene])
        queries.append({"id": gene, "input_ordinal_0based": ids.index(gene), **record(path)})
    combined = output / "combined.fa"
    with combined.open("xb") as stream:
        for gene in selected:
            stream.write(raw[gene])
    jobs = []
    for name, query in [("combined", combined), *[(f"single_{i:02d}", Path(q["path"]))
                                                for i, q in enumerate(queries)]]:
        command = list(plan["search_commands"]["blast"])
        command[command.index("-i") + 1] = str(query)
        command[command.index("-o") + 1] = str(output / (name + ".blast"))
        jobs.append({"name": name, "argv": command})
    report = {"status": "diagnostic_panel_prepared_not_executed", "source": record(__file__),
              "biopython_version": Bio.__version__, "plan": record(plan_path),
              "runtime": record(runtime_path), "input": source, "database": database,
              "queries": queries, "combined": record(combined), "commands": jobs,
              "cpus_per_task": 180, "search_admitted": False, "reuse_authorized": False,
              "limitations": ["No replay executed; no output equivalence established.",
                              "Small diagnostic panel does not validate the full retained prefix."]}
    for item in database:
        check(item)
    read_frozen(plan_path, PLAN_SHA)
    with (output / "panel.json").open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("plan", "runtime", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    report = prepare(args.plan.resolve(), args.runtime.resolve(), args.output.resolve())
    print(json.dumps({"status": report["status"], "queries": report["queries"]}))
