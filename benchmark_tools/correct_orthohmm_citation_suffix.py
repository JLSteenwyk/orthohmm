"""Apply one manually reviewed author-source correction without changing Crossref."""

import argparse
from copy import deepcopy
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

RAW_SHA = "c44a7eb5b6eec9e4fa5d0c13146d194fe932888931fd8adb62fa66095f45f3e8"
AUTHOR_SHA = "d713be57a6457ba281f9a275f290df253ab095692de61c9aba4c58431b1d5a25"
API_SHA = "fbbf4f22c8b2f98bada38e48ffdd7e1eb1efb5864edb6c03ffad4737bcf0d933"
DOI = "10.1101/2024.12.07.627370"


def corrected(records):
    result = deepcopy(records)
    targets = [r for r in result if r.get("id") == "orthohmm2024preprint"]
    if len(targets) != 1:
        raise ValueError("Require one selected OrthoHMM citation")
    target = targets[0]
    if (target.get("DOI") != DOI or target.get("genre") != "preprint"
            or target.get("author") != [
                {"given": "Jacob L", "family": "Steenwyk"},
                {"given": "Thomas J.", "family": "Buida"},
                {"given": "Antonis", "family": "Rokas"},
                {"given": "Nicole", "family": "King"}]):
        raise ValueError("Citation DOI, type or deposited byline differs")
    target["author"][1]["suffix"] = "III"
    return result


def export(raw, author_page, api, output, provenance):
    if (output.resolve() == provenance.resolve() or output.exists() or provenance.exists()
            or output.is_symlink() or provenance.is_symlink()):
        raise FileExistsError("Require distinct fresh output paths")
    evidence = [record(raw), record(author_page), record(api)]
    if [e["sha256"] for e in evidence] != [RAW_SHA, AUTHOR_SHA, API_SHA]:
        raise ValueError("Reviewed citation or source bytes changed")
    before = json.loads(raw.read_text())
    after = corrected(before)
    source = record(__file__)
    for item in [source, *evidence]:
        check(item)
    with output.open("x") as stream:
        json.dump(after, stream, indent=2, ensure_ascii=True)
        stream.write("\n")
    report = dict(status="author_source_supported_orthohmm_suffix_exported", raw=evidence[0],
        author_page=evidence[1], publisher_api=evidence[2], output=record(output), source=source,
        citation_id="orthohmm2024preprint", doi=DOI,
        change={"author_index": 1, "before": {"given": "Thomas J.", "family": "Buida"},
                "after": {"given": "Thomas J.", "family": "Buida", "suffix": "III"}},
        source_urls=["https://jlsteenwyk.com/publications.html",
                     "https://api.biorxiv.org/details/biorxiv/" + DOI],
        limitations=["Manual review of the author publication list supports this suffix; not a parsed publisher byline.",
            "Crossref and bioRxiv API omit the suffix; both original snapshots remain unchanged.",
            "Only one suffix field changes. No title, date, author order, other record, or preprint status changes.",
            "The bioRxiv article page returned HTTP 403 during review; no access restriction was bypassed.",
            "Not journal rendering, full-text scientific validation or executed-method provenance."])
    with provenance.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("raw", "author-page", "api", "output", "provenance"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    export(args.raw, args.author_page, args.api, args.output, args.provenance)
