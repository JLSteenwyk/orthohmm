"""Compare terminal legacy BLAST diagnostics without authorizing prefix reuse."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess

from benchmark_tools.audit_orthomcl_blast import parse_diagnostics
from benchmark_tools.audit_orthomcl_search_table import audit_table
from benchmark_tools.prepare_blast_recovery_partition import AUDIT_SHA
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_blast_replay_panel import PANEL_SHA, query_file_records, validate_commands
from benchmark_tools.run_qfo_corrected_blast import verify
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.verify_ygob_validation import require_completed_job


def hsp_rows(content, allowed):
    rows = {gene: Counter() for gene in allowed}
    if content and not content.endswith(b"\n"):
        raise ValueError("Incomplete replay table")
    for line in content.splitlines():
        fields = line.split(b"\t")
        if len(fields) != 12 or b"\x00" in line:
            raise ValueError("Invalid HSP row")
        gene = fields[0].decode("ascii")
        if gene not in rows:
            raise ValueError("Unexpected replay query")
        rows[gene][tuple(field.decode("ascii") for field in fields)] += 1
    return rows


def messages(diagnostics, gene):
    value = diagnostics.get(gene)
    return Counter((r["level"], r["category"], r["message"]) for r in value["messages"]) if value else Counter()


def compare_query(gene, combined, single, original, incomplete, combined_log, single_log, original_log):
    combined_equal = combined == single
    original_relation = ("unobserved" if original is None else
        "partial_subset" if incomplete and not (original - combined) else
        "equal" if not incomplete and original == combined else "mismatch")
    log_equal = messages(combined_log, gene) == messages(single_log, gene)
    original_log_equal = messages(original_log, gene) == messages(combined_log, gene)
    return dict(query=gene, combined_hsp_rows=sum(combined.values()), single_hsp_rows=sum(single.values()),
        original_observed_hsp_rows=None if original is None else sum(original.values()),
        combined_single_hsp_multisets_equal=combined_equal,
        combined_single_diagnostics_equal=log_equal, original_hsp_relation=original_relation,
        original_diagnostics_equal=original_log_equal,
        diagnostic_compatible=combined_equal and log_equal and original_log_equal and original_relation != "mismatch")


def compare(root, output):
    root = root.resolve()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    results = root / "benchmark_tools/results"
    panel_path = results / "qfo_blast_replay_panel_20260923.json"
    panel = read_frozen(panel_path, PANEL_SHA)
    audit_path = root / "benchmarks/work/qfo_blast_prefix_audit_20260923/status.json"
    audit = read_frozen(audit_path, AUDIT_SHA)
    plan = verify(Path(panel["plan"]["path"]), Path(panel["runtime"]["path"]))
    directory = validate_commands(panel, plan["search_commands"]["blast"])
    status_path = directory / "execution/status.json"
    status = json.loads(status_path.read_text())
    if (status["status"] != "native_replays_finished_pending_comparison" or status["job_id"] != "22055"
            or status["panel"] != record(panel_path) or status["commands"] != panel["commands"]
            or status["search_admitted"] is not False or status["reuse_authorized"] is not False):
        raise ValueError("Require terminal frozen replay panel")
    accounting = subprocess.check_output(["sacct", "-j", "22055", "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed"], text=True)
    scheduler = require_completed_job(accounting, 22055)
    if [s["name"] for s in status["stages"]] != [c["name"] for c in panel["commands"]]:
        raise ValueError("Incomplete replay stage inventory")
    checked = [record(panel_path), record(audit_path), record(status_path), record(__file__),
        status["source"], panel["source"], panel["plan"], panel["runtime"], panel["input"],
        panel["combined"], *query_file_records(panel["queries"]), *panel["database"], audit["blocks"]]
    checked.extend(record(Path(__file__).with_name(name)) for name in (
        "audit_orthomcl_blast.py", "audit_orthomcl_search_table.py", "convert_orthomcl_blast.py",
        "run_blast_replay_panel.py", "run_qfo_corrected_blast.py", "verify_ygob_validation.py"))
    tables, logs = {}, {}
    ids = [q["id"] for q in panel["queries"]]
    for n, stage in enumerate(status["stages"]):
        name = stage["name"]
        if (stage["exit_code"] != 0 or stage["finished_epoch"] < stage["started_epoch"]
                or stage["output"]["path"] != str(directory / (name + ".blast"))
                or stage["log"]["path"] != str(directory / "execution" / (name + ".log"))):
            raise ValueError("Failed replay stage or unexpected paths")
        checked.extend([stage["output"], stage["log"]])
    for item in checked:
        check(item)
    for n, stage in enumerate(status["stages"]):
        name = stage["name"]
        blast, log = Path(stage["output"]["path"]), Path(stage["log"]["path"])
        allowed = ids if n == 0 else [ids[n-1]]
        tables[name] = hsp_rows(blast.read_bytes(), allowed)
        logs[name] = parse_diagnostics(log)
        if not set(logs[name]) <= set(allowed):
            raise ValueError("Unexpected diagnostic query")
        structural = audit_table(blast, Path(panel["input"]["path"]), log) if blast.stat().st_size else None
        if structural and structural["hsp_rows_above_1e_minus_5"]:
            raise ValueError("Replay HSP exceeds frozen cutoff")
    blocks = {}
    with Path(audit["blocks"]["path"]).open() as stream:
        for line in stream:
            row = json.loads(line)
            if row["query"] in ids:
                if row["query"] in blocks:
                    raise ValueError("Repeated retained diagnostic query")
                blocks[row["query"]] = row
    partial = Path(plan["output_root"]) / "work/all.blast.partial"
    partial_record = record(partial)
    if partial_record["sha256"] != audit["boundary"]["sha256"]:
        raise ValueError("Interrupted BLAST bytes changed")
    checked.append(partial_record)
    original = {}
    with partial.open("rb") as stream:
        for gene, block in blocks.items():
            stream.seek(block["start"])
            content = stream.read(block["end"] - block["start"])
            if hashlib.sha256(content).hexdigest() != block["sha256"]:
                raise ValueError("Retained query block changed")
            original[gene] = hsp_rows(content, [gene])[gene]
            if sum(original[gene].values()) != block["rows"]:
                raise ValueError("Retained row count differs")
    original_log = {r["gene"]: r for r in audit["content"]["diagnostics"]}
    comparisons = [compare_query(gene, tables["combined"][gene], tables[f"single_{i:02d}"][gene],
        original.get(gene), gene == audit["boundary"]["last_query"], logs["combined"],
        logs[f"single_{i:02d}"], original_log) for i, gene in enumerate(ids)]
    for item in checked:
        check(item)
    report = dict(status="diagnostic_replay_comparison_complete", scheduler=scheduler, records=checked,
        comparisons=comparisons, all_diagnostics_compatible=all(r["diagnostic_compatible"] for r in comparisons),
        search_admitted=False, reuse_authorized=False, publication_ready=False,
        limitations=["Exact printed HSP-field multisets retain duplicates but ignore row order; structural checks separately require contiguous blocks.",
            "A partial-boundary subset match never establishes complete boundary output.",
            "Five-query agreement cannot establish full-prefix durable completion or authorize recovery.",
            "Any mismatch requires investigation; no automatic replay, merge, or downstream admission is performed."])
    with output.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    compare(args.root, args.output)
