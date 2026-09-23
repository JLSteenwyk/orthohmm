"""Check historical annotation evidence and freeze descriptive family bins."""

import argparse
from datetime import datetime
import hashlib
import json
from pathlib import Path
import subprocess

from Bio import SeqIO, SwissProt
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.verify_ygob_validation import require_completed_job

INVENTORY_SHA = "912363276fa5456f11aeff9f398fce71aec8d377c8c8eda7b144105945e554f1"
HELPER_SHA = "0642beeee11199853518b328b15bdb12d297d3a1324bcd3327815ccbc7355680"
COLLECTOR_SHA = "1c672aa338f22070bb8863297dd8c54fe4609b335f9f5b264d008bdfe7cdcac9"
PROTOCOL_SHA = "447714557c357e752152ab4d2c95e624d7250195a7f0a59a1a1995c19f7f8359"


def selected_history(history, accession, sv):
    rows = history["results"]
    if any(r["accession"] != accession for r in rows) or len({r["entryVersion"] for r in rows}) != len(rows):
        raise ValueError("Ambiguous accession history")
    cutoff = datetime(2020, 8, 12)
    active, later = [], []
    for row in rows:
        first = datetime.strptime(row["firstReleaseDate"], "%d-%b-%Y")
        last = datetime.strptime(row["lastReleaseDate"], "%d-%b-%Y")
        if last < first:
            raise ValueError("Reversed release interval")
        if row["sequenceVersion"] == sv:
            if first <= cutoff <= last:
                active.append(row)
            elif first > cutoff:
                later.append((first, row["entryVersion"], row))
    if len(active) > 1:
        raise ValueError("Overlapping baseline entry versions")
    if active:
        return active[0], "baseline_release"
    # A version that existed before baseline but not at baseline is not a later update.
    if any(r["sequenceVersion"] == sv and datetime.strptime(r["firstReleaseDate"], "%d-%b-%Y") <= cutoff for r in rows):
        raise ValueError("Sequence version absent at baseline")
    if later:
        return sorted(later, key=lambda r: r[:2])[0][2], "later_sequence_version"
    raise ValueError("No applicable annotation version")


def family_bins(families, annotations, baseline_only=False):
    wanted = {g for members in families.values() for g in members}
    if set(annotations) != wanted or any(not m or len(m) != len(set(m)) for m in families.values()):
        raise ValueError("Incomplete family annotation universe")
    bins = dict(annotation_positive=[], all_matched_unflagged=[], missing_without_positive=[])
    for family, members in sorted(families.items()):
        states = []
        for gene in members:
            row = annotations[gene]
            if row is None or (baseline_only and row["selection_class"] != "baseline_release"):
                states.append(None)
            else:
                states.append(row["fragment_flag"] or bool(row["incomplete_sequence_features"]))
        key = ("annotation_positive" if any(v is True for v in states) else
               "missing_without_positive" if any(v is None for v in states) else "all_matched_unflagged")
        bins[key].append(family)
    return bins


def verify(root, output):
    root = root.resolve()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    directory = root / "benchmarks/work/swiss_historical_fragment_panel_20260923"
    status_path = directory / "status.json"
    report = json.loads(status_path.read_text())
    if (report["status"] != "collection_complete_pending_independent_validation" or report["job_id"] != "22116"
            or any(report[k] is not False for k in ("prediction_statistics_evaluated", "annotation_panel_admitted", "publication_ready"))):
        raise ValueError("Incomplete or incorrectly labeled collection")
    accounting = subprocess.check_output(["sacct", "-j", "22116", "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed"], text=True)
    scheduler = require_completed_job(accounting, 22116)
    inventory_path = root / "benchmark_tools/results/corrected_swiss_sequence_strata_20260918.json"
    inventory = read_frozen(inventory_path, INVENTORY_SHA)
    preflight = report["preflight"]
    families, genes = inventory["family_memberships"], preflight["genes"]
    wanted = set(inventory["genes"])
    if (len(wanted) != 563 or len(families) != 18 or families != preflight["families"]
            or set(genes) != wanted or set(report["entries"]) != wanted):
        raise ValueError("Changed full panel membership")
    tracked = [record(status_path), record(inventory_path), record(__file__), *preflight["records"]]
    if ([r["sha256"] for r in preflight["records"][:4]] !=
            [INVENTORY_SHA, PROTOCOL_SHA, HELPER_SHA, COLLECTOR_SHA]
            or preflight["records"][4:] != inventory["fasta_inputs"]):
        raise ValueError("Changed input protocol or executable identities")
    for item in tracked:
        check(item)
    observed = set()
    for item in inventory["fasta_inputs"]:
        for entry in SeqIO.parse(item["path"], "fasta"):
            gene = entry.id.split("|")[1]
            if gene not in wanted:
                continue
            row = genes[gene]
            descriptor = inventory["genes"][gene]
            tags = dict(token.split("=", 1) for token in entry.description.split() if "=" in token)
            if (gene in observed or row["accession"] != gene or row["length"] != len(entry.seq)
                    or row["sequence_sha256"] != hashlib.sha256(str(entry.seq).encode("ascii")).hexdigest()
                    or row["taxid"] != tags["OX"] or row["sequence_version"] != int(tags["SV"])
                    or descriptor["description"] != entry.description):
                raise ValueError("Historical query identity differs from staged FASTA")
            observed.add(gene)
    if observed != wanted:
        raise ValueError("Missing staged sequences")
    annotations = {}
    for gene, state in sorted(report["entries"].items()):
        if state["status"] == "missing":
            if not state.get("error_type") or not state.get("error"):
                raise ValueError("Unexplained missing annotation")
            annotations[gene] = None
            continue
        if state["status"] != "sequence_matched" or state["audit"]["path"] != str(directory / gene / "audit.json"):
            raise ValueError("Unknown annotation state or path")
        check(state["audit"])
        audit = json.loads(Path(state["audit"]["path"]).read_text())
        expected_paths = [directory / gene / name for name in ("history.json", "entry.txt")]
        if audit["records"] != [record(p) for p in expected_paths] or audit["source"]["sha256"] != HELPER_SHA:
            raise ValueError("Annotation raw source identities changed")
        tracked.extend([state["audit"], *audit["records"]])
        chosen, kind = selected_history(json.loads(expected_paths[0].read_text()), gene, genes[gene]["sequence_version"])
        with expected_paths[1].open() as stream:
            entry = SwissProt.read(stream)
        flags = []
        for field in entry.description.split(";"):
            if field.strip().startswith("Flags: "):
                flags.append(field.strip()[7:])
        features = [dict(type=f.type, location=str(f.location), qualifiers=f.qualifiers)
                    for f in entry.features if f.type in {"NON_TER", "NON_CONS"}]
        expected = dict(**genes[gene], entry_version=chosen["entryVersion"], selection=chosen,
            selection_class=kind, annotation_date=entry.annotation_update[0], data_class=entry.data_class,
            flags=flags, fragment_flag=any(f in ("Fragment", "Fragments") for f in flags),
            incomplete_sequence_features=features, outcome_analysis_performed=False)
        if (entry.accessions[0] != gene or entry.taxonomy_id != [genes[gene]["taxid"]]
                or entry.sequence_update[1] != genes[gene]["sequence_version"]
                or entry.annotation_update[1] != chosen["entryVersion"]
                or hashlib.sha256(entry.sequence.encode("ascii")).hexdigest() != genes[gene]["sequence_sha256"]
                or any(audit[k] != v for k, v in expected.items()) or state["selection_class"] != kind):
            raise ValueError("Raw annotation/receipt disagreement")
        annotations[gene] = expected
    matched = sum(r is not None for r in annotations.values())
    if (matched, len(wanted)-matched) != (report["matched"], report["missing"]):
        raise ValueError("Collection accounting mismatch")
    for item in tracked:
        check(item)
    result = dict(status="historical_annotation_panel_checked_with_explicit_missingness", scheduler=scheduler,
        records=tracked, families=families, annotations=annotations, matched=matched, missing=len(wanted)-matched,
        strata=family_bins(families, annotations), baseline_only_strata=family_bins(families, annotations, True),
        annotation_panel_admitted=True, prediction_statistics_evaluated=False, publication_ready=False,
        limitations=["Shares Bio.SwissProt parser with acquisition; this is not independent parser validation.",
            "Missing acquisition records remain missing; their external cause is not independently established.",
            "Unflagged annotations do not prove complete proteins; historical labels may contain errors."])
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    verify(args.root, args.output)
