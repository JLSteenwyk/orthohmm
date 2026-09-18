"""Extract a legacy protein database and compare every record with its input."""

import argparse
from itertools import zip_longest
import json
from pathlib import Path
import subprocess
import sys

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_qfo_input_sequences import sequence_identity
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.prepare_qfo_corrected_orthomcl import SOFTWARE
from benchmark_tools.run_qfo_corrected_blast import RUNTIME_SHA, environment
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.snapshot_runtime_trees import verify as verify_tree

FASTACMD = SOFTWARE / "blast-2.2.13/bin/fastacmd"


def compare_dump(fasta, dump):
    seen, differences = set(), []
    count = input_residues = database_residues = 0
    with fasta.open() as original, dump.open() as extracted:
        for ordinal, (source, native) in enumerate(zip_longest(
                SeqIO.parse(original, "fasta"), SeqIO.parse(extracted, "fasta"))):
            if source is None or native is None:
                raise ValueError("Database and input sequence counts differ")
            if source.id in seen:
                raise ValueError("Duplicate input sequence ID")
            seen.add(source.id)
            # The frozen formatdb command does not request parsed sequence IDs.
            expected_title = f"gnl|BL_ORD_ID|{ordinal} {source.description}"
            if native.description != expected_title:
                raise ValueError(f"Database header/order differs at ordinal {ordinal}")
            expected, observed = sequence_identity(str(source.seq)), sequence_identity(str(native.seq))
            count += 1
            input_residues += expected["length"]
            database_residues += observed["length"]
            if expected != observed:
                differences.append({"id": source.id, "ordinal": ordinal,
                                    "input": expected, "database": observed})
    if count == 0:
        raise ValueError("Empty input/database")
    return {"input_sequences": count, "database_sequences": count,
            "input_residues": input_residues, "database_residues": database_residues,
            "exact_sequence_matches": count - len(differences),
            "sequence_difference_count": len(differences), "differences": differences,
            "exact_sequence_parity": not differences, "header_order_verified": True}


def database_files(fasta):
    paths = sorted(fasta.parent.glob(fasta.name + ".*"))
    expected = {fasta.name + "." + suffix for suffix in ("phr", "pin", "psq")}
    if {p.name for p in paths} != expected or any(
            p.is_symlink() or not p.is_file() or not p.stat().st_size for p in paths):
        raise ValueError("Require the three nonempty unsplit native database files")
    return paths


def audit(fasta, runtime_path, output):
    if output.exists():
        raise FileExistsError(output)
    runtime = read_frozen(runtime_path, RUNTIME_SHA)
    verify_tree(runtime)
    tool = record(FASTACMD)
    if not any(all(item.get(k) == tool[k] for k in ("path", "bytes", "sha256"))
               for item in runtime["records"] if item["kind"] == "file"):
        raise ValueError("Extractor not bound to frozen runtime")
    for config in (Path.home() / ".ncbirc", Path("/etc/.ncbirc"), Path("/etc/ncbi.ini"),
                   Path("/etc/ncbi"), Path("/etc/ld.so.preload")):
        if config.exists() or config.is_symlink():
            raise ValueError("Unreviewed configuration or loader override: " + str(config))
    db = database_files(fasta)
    checked = [record(fasta), record(runtime_path), tool, *[record(p) for p in db],
               record(__file__), record(Path(__file__).with_name("audit_qfo_input_sequences.py")),
               record(Path(__file__).with_name("run_qfo_corrected_blast.py")),
               record(Path(__file__).with_name("snapshot_runtime_trees.py"))]
    output.mkdir(parents=True, exist_ok=False)
    dump, log = output / "database.fasta", output / "fastacmd.log"
    argv = [str(FASTACMD), "-d", str(fasta), "-p", "T", "-D", "1"]
    report = {"status": "running", "command": argv, "cwd": str(output),
              "environment": environment(), "checked_records": checked,
              "database_admitted": False, "search_admitted": False, "accuracy_admitted": False,
              "publication_ready": False}
    try:
        with dump.open("xb") as stream, log.open("xb") as errors:
            done = subprocess.run(argv, cwd=output, env=report["environment"],
                                  stdout=stream, stderr=errors)
        report["exit_code"] = done.returncode
        if done.returncode or log.stat().st_size:
            raise ValueError("Database extraction failed or emitted unreviewed diagnostics")
        report["content"] = compare_dump(fasta, dump)
        if database_files(fasta) != db:
            raise ValueError("Database file membership changed")
        for item in checked:
            check(item)
        verify_tree(runtime)
        report["status"] = ("database_exact_sequence_parity_verified"
                            if report["content"]["exact_sequence_parity"]
                            else "database_sequence_differences_require_review")
        report["limitations"] = [
            "Extraction checks current database bytes, not their use by a particular BLAST job.",
            "Terminal search provenance and before/after database identity must be bound separately.",
            "Sequence differences are retained without case conversion or residue normalization.",
            "Only the frozen default, unsplit, unparsed-ID protein database layout is supported.",
            "Extraction success does not admit search completeness, downstream BPO/clustering or accuracy."]
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        report["outputs"] = [record(p) for p in (dump, log) if p.is_file()]
        with (output / "report.json").open("x") as stream:
            json.dump(report, stream, indent=2, sort_keys=True)
            stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("fasta", "runtime", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.fasta.resolve(), args.runtime.resolve(), args.output.resolve())
    print(json.dumps({"status": result["status"], "content": result["content"]}))
