"""Combine explicitly selected, checksum-bound CSL exports without altering records."""

import argparse
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def assemble(manifest, output, provenance):
    if output.resolve() == provenance.resolve():
        raise ValueError("Output and provenance must be distinct")
    if output.exists() or provenance.exists():
        raise FileExistsError("Output or provenance already exists")
    selection = json.loads(manifest.read_text())
    if not isinstance(selection, list) or not selection:
        raise ValueError("Expected a nonempty source selection")
    citations, sources, identities = [], [], []
    files, ids, dois = set(), set(), set()
    for item in selection:
        name = item["file"]
        if not isinstance(name, str) or Path(name).name != name or name in ("", ".", ".."):
            raise ValueError("Source must be a local filename")
        if name in files:
            raise ValueError("Duplicate source file")
        files.add(name)
        path = manifest.parent / name
        raw = path.read_bytes()
        digest = hashlib.sha256(raw).hexdigest()
        if digest != item["sha256"]:
            raise ValueError("Source checksum differs: " + name)
        rows = json.loads(raw)
        if not isinstance(rows, list) or not rows:
            raise ValueError("Expected a nonempty CSL list")
        for row in rows:
            if not isinstance(row, dict) or any(
                not isinstance(row.get(k), str) or not row[k].strip()
                for k in ("id", "type", "title")
            ):
                raise ValueError("Missing citation identity/type/title")
            if row["id"] in ids:
                raise ValueError("Duplicate citation ID: " + row["id"])
            ids.add(row["id"])
            if "DOI" in row:
                doi = row["DOI"]
                if not isinstance(doi, str) or not doi.strip():
                    raise ValueError("Invalid DOI")
                normalized = doi.strip().casefold()
                if normalized in dois:
                    raise ValueError("Duplicate DOI: " + doi)
                dois.add(normalized)
            citations.append(row)
            identities.append({"id": row["id"], "source": name})
        sources.append({"file": name, "sha256": digest, "bytes": len(raw), "records": len(rows)})
    with output.open("x") as handle:
        handle.write(json.dumps(citations, indent=2, ensure_ascii=True) + "\n")
    result = {
        "status": "selected_citations_assembled",
        "manifest": record(manifest), "assembler": record(__file__),
        "output": record(output), "sources": sources, "records": identities,
        "count": len(citations),
        "limitations": [
            "Explicit selected exports only; not a claim of complete bibliography coverage.",
            "No CSL fields, bylines, dates or article identifiers were corrected or normalized.",
            "Not journal-style rendering, scientific source review, version provenance or rights clearance.",
        ],
    }
    with provenance.open("x") as handle:
        handle.write(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("manifest", "output", "provenance"):
        parser.add_argument("--" + name, required=True, type=Path)
    args = parser.parse_args()
    print(f"Assembled {assemble(args.manifest, args.output, args.provenance)['count']} citations")
