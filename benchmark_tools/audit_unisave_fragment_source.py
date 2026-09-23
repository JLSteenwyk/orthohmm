"""Acquire a sequence-matched historical annotation, without reading predictions."""

import argparse
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import re
from urllib.request import urlopen

from Bio import SwissProt
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def choose_version(history, accession, sequence_version):
    cutoff = datetime(2020, 8, 12)
    rows = history["results"]
    if not rows or any(r["accession"] != accession for r in rows):
        raise ValueError("Missing or mixed accession history")
    if len({r["entryVersion"] for r in rows}) != len(rows):
        raise ValueError("Duplicate entry versions")
    date = lambda r: datetime.strptime(r["firstReleaseDate"], "%d-%b-%Y")
    eligible = [r for r in rows if r["sequenceVersion"] == sequence_version]
    before = [r for r in eligible if date(r) <= cutoff]
    if before:
        chosen = max(before, key=lambda r: (date(r), r["entryVersion"]))
        if datetime.strptime(chosen["lastReleaseDate"], "%d-%b-%Y") < cutoff:
            raise ValueError("Matching sequence version not present at baseline release")
        return chosen, "baseline_release"
    if eligible:
        return min(eligible, key=lambda r: (date(r), r["entryVersion"])), "later_sequence_version"
    raise ValueError("No matching sequence version")


def inspect_entry(text, accession, version, sequence_version, sequence_sha256, taxid):
    entry = SwissProt.read(io.StringIO(text))
    digest = hashlib.sha256(entry.sequence.encode("ascii")).hexdigest()
    if (entry.accessions[0] != accession or entry.annotation_update[1] != version
            or entry.sequence_update[1] != sequence_version or digest != sequence_sha256
            or entry.taxonomy_id != [str(taxid)]):
        raise ValueError("Historical annotation does not match exact input identity")
    flags = re.findall(r"Flags: ([^;]+);", entry.description)
    features = [dict(type=f.type, location=str(f.location), qualifiers=f.qualifiers)
                for f in entry.features if f.type in {"NON_TER", "NON_CONS"}]
    return dict(accession=accession, entry_version=version, sequence_version=sequence_version,
        sequence_sha256=digest, taxid=str(taxid), length=len(entry.sequence),
        annotation_date=entry.annotation_update[0], data_class=entry.data_class,
        flags=flags, fragment_flag=any(f in {"Fragment", "Fragments"} for f in flags),
        incomplete_sequence_features=features,
        limitation="Unflagged does not prove completeness; these annotations are not experimental ground truth.")


def acquire(accession, sequence_version, sequence_sha256, taxid, output):
    if not re.fullmatch(r"[A-Z0-9]+", accession) or not re.fullmatch(r"[0-9a-f]{64}", sequence_sha256):
        raise ValueError("Invalid accession or sequence digest")
    output.mkdir(parents=True, exist_ok=False)
    urls = [f"https://rest.uniprot.org/unisave/{accession}?format=json"]

    def download(url, name):
        with urlopen(url, timeout=60) as response:
            payload = response.read()
        path = output / name
        with path.open("xb") as stream:
            stream.write(payload)
        return payload

    history = json.loads(download(urls[0], "history.json"))
    selected, provenance = choose_version(history, accession, sequence_version)
    urls.append(f"https://rest.uniprot.org/unisave/{accession}?format=txt&versions={selected['entryVersion']}")
    content = download(urls[1], "entry.txt").decode("utf-8")
    result = inspect_entry(content, accession, selected["entryVersion"], sequence_version, sequence_sha256, taxid)
    result.update(selection=selected, selection_class=provenance, urls=urls,
        acquired_utc=datetime.now(timezone.utc).isoformat(),
        records=[record(output / n) for n in ("history.json", "entry.txt")],
        source=record(__file__), outcome_analysis_performed=False)
    with (output / "audit.json").open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accession", required=True)
    parser.add_argument("--sequence-version", required=True, type=int)
    parser.add_argument("--sequence-sha256", required=True)
    parser.add_argument("--taxid", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(acquire(args.accession, args.sequence_version, args.sequence_sha256, args.taxid, args.output)))
