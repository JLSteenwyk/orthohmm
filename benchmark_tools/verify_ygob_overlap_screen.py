"""Recheck retained YGOB homology-screen evidence without reading predictions."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_historical_profile_ablation import verify_file
from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.screen_ygob_homology import search_command, summarize_hits
from benchmark_tools.verify_ygob_validation import require_completed_job


def verify(root):
    root = root.resolve()
    report_path = root / "benchmark_tools/results/ygob_homology_screen_20260916.json"
    report = json.loads(report_path.read_text())
    if report["status"] != "complete":
        raise ValueError("Incomplete homology screen")
    provenance = report["provenance"]
    labels = ["qfo_2020", "orthobench", "three_kingdoms"]
    if set(provenance["development_inputs"]) != set(labels) or provenance["target_counts"] != [976504, 251378, 443217]:
        raise ValueError("Unexpected development panel")
    for item in [report["hits"], report["target_mapping"], provenance["source"],
                 provenance["manifest"], provenance["diamond_binary"]]:
        verify_file(Path(item["path"]), item)
    for records in provenance["development_inputs"].values():
        for item in records:
            verify_file(Path(item["path"]), item)
    if provenance["diamond_version"] != "diamond version 2.1.11":
        raise ValueError("Wrong screen version")
    work = Path("benchmarks/work/ygob_homology_screen_v1")
    expected = search_command(provenance["diamond_binary"]["path"], work / "queries.fasta",
                              work / "development", work / "hits.tsv", 32)
    if provenance["commands"][1] != expected:
        raise ValueError("Screen command differs from frozen protocol")
    prepared = root / "benchmarks/work/ygob_validation_v1"
    manifest = json.loads((prepared / "manifest.json").read_text())
    query_species = {}
    for item in manifest["inputs"]:
        path = Path(item["path"])
        verify_file(path, item)
        for record in SeqIO.parse(path, "fasta"):
            if record.id in query_species:
                raise ValueError("Duplicate query ID")
            query_species[record.id] = path.stem
    verify_file(Path(manifest["reference"]["path"]), manifest["reference"])
    references = json.loads(Path(manifest["reference"]["path"]).read_text())
    recomputed = summarize_hits(Path(report["hits"]["path"]), query_species, references,
                                provenance["target_counts"], labels)
    if any(report[key] != value for key, value in recomputed.items()):
        raise ValueError("Retained overlap summary does not reproduce")
    accounting = subprocess.check_output(["sacct", "-j", "20918", "--parsable2",
                                         "--format=JobIDRaw,State,ExitCode,Elapsed"], text=True)
    scheduler = require_completed_job(accounting, 20918)
    frozen = root / "benchmarks/work/publication_method_native_v2"
    installed = root.parent.parent / "SOFTWARE/orthofinder_3.1.5/lib/python3.12/site-packages/orthofinder"
    sources = {
        frozen / "orthohmm/search/profile_expansion.py": "1025cfab4f80e3a85d6eaa3e768ae5e0e25dedb48fdba9deb08cb6e6aa7bb42a",
        frozen / "orthohmm/phylogeny_pipeline.py": "44e00316546b5df78354badd2b1a6bb595b685e98b86f686f1a32543d7f15e4f",
        installed / "run/config.json": "10a10a93262c9676865f03ba05b2b77d48285c312dca435d2458939894a474b3",
        installed / "run/run_commands.py": "10aeb00c3affb31175a337478957807b79b153b9b91597b29e4dccb782d8d387"}
    source_records = [file_provenance(path) for path in sources]
    if any(item["sha256"] != sources[Path(item["path"])] for item in source_records):
        raise ValueError("Reference-resource audited source changed")
    return {"schema_version": 1, "status": "overlap_evidence_verified", "accuracy_evaluated": False,
            "scheduler": scheduler, "screen": file_provenance(report_path),
            "matching_proteins": recomputed["matching_query_proteins"], "query_proteins": len(query_species),
            "matching_pillars": recomputed["reference_pillars_with_hit"], "reference_pillars": len(references),
            "reference_resource_audit": file_provenance(root / "benchmark_tools/results/YGOB_REFERENCE_RESOURCE_AUDIT_20260916.md"),
            "reviewed_sources": source_records, "verifier": file_provenance(Path(__file__)),
            "permitted_interpretation": "Frozen novel-taxon curated-group transfer; NOT family-disjoint validation.",
            "limitations": [*recomputed["limitations"],
                            "The source/resource audit is not a syscall trace or exhaustive annotation-ancestry audit.",
                            "Raw database redistribution permissions remain unresolved.",
                            "This rechecks retained inputs/hits, not a fresh search or historical proof of intermediate database bytes."]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = verify(args.root)
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")


if __name__ == "__main__":
    main()
