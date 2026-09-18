"""Create a separate, source-bound CSL byline correction; preserve raw Crossref."""

import argparse
from copy import deepcopy
import json
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check

RAW_SHA = "8c844fbef9a2694da284a2e98e22baf78a5b865465bcabd7768b321e877fc53a"
SOURCES = {
    "qfo2020": ("PMC7319555", "50162a498e935d807943160e706d3ed3f6f6f5769313277e9de4eea0fe4074ec", 23),
    "qfo2022": ("PMC9252809", "d6b9ae310a6fd7bf617d6011e4649313dbe1992b4df5dd76af1ff57231f50255", 31),
}


def text(element):
    return " ".join("".join(element.itertext()).split()) if element is not None else ""


def byline(xml, doi):
    meta = ET.fromstring(xml).find("./front/article-meta")
    if meta is None or [text(e) for e in meta.findall("article-id[@pub-id-type='doi']")] != [doi]:
        raise ValueError("Article DOI differs or is missing")
    groups = meta.findall("contrib-group")
    if len(groups) != 1:
        raise ValueError("Require one explicit top-level contributor group")
    authors = []
    for contributor in groups[0].findall("contrib"):
        if contributor.get("contrib-type", "author") != "author":
            raise ValueError("Unexpected non-author contributor")
        names, collectives = contributor.findall("name"), contributor.findall("collab")
        if len(names) + len(collectives) != 1:
            raise ValueError("Ambiguous contributor name")
        if collectives:
            if list(collectives[0]):
                raise ValueError("Unexpected nested collective markup")
            author = {"literal": text(collectives[0])}
        else:
            author = {"given": text(names[0].find("given-names")),
                      "family": text(names[0].find("surname"))}
            # Reviewed 2020 source encodes this collective as a personal name.
            if doi == "10.1093/nar/gkaa308" and author == {
                    "given": "the Quest", "family": "for Orthologs Consortium"}:
                author = {"literal": "the Quest for Orthologs Consortium"}
        if not all(author.values()):
            raise ValueError("Empty contributor name")
        authors.append(author)
    if not authors:
        raise ValueError("Empty byline")
    return authors


def export(raw, xml_directory, output, provenance):
    if output.exists() or provenance.exists() or output.resolve() == provenance.resolve():
        raise FileExistsError("Require distinct fresh output paths")
    original = record(raw)
    if original["sha256"] != RAW_SHA:
        raise ValueError("Unexpected raw CSL export")
    records = json.loads(raw.read_text())
    if len(records) != 2 or {r["id"] for r in records} != set(SOURCES):
        raise ValueError("Wrong citation selection")
    corrected, evidence = deepcopy(records), []
    for item in corrected:
        pmc, sha, expected_count = SOURCES[item["id"]]
        xml_path = xml_directory / (pmc + ".xml")
        source = record(xml_path)
        if source["sha256"] != sha:
            raise ValueError("Article XML changed")
        previous = item["author"]
        item["author"] = byline(xml_path.read_bytes(), item["DOI"])
        if len(item["author"]) != expected_count:
            raise ValueError("Reviewed byline count differs")
        evidence.append(dict(id=item["id"], source=source,
            url=f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/fullTextXML",
            original_author_count=len(previous), corrected_author_count=len(item["author"]),
            original_authors=previous, corrected_authors=item["author"]))
    for source in [original, *(e["source"] for e in evidence)]:
        check(source)
    output.write_text(json.dumps(corrected, indent=2, ensure_ascii=True) + "\n")
    provenance.write_text(json.dumps(dict(status="reviewed_qfo_service_bylines_exported",
        raw=original, output=record(output), source=record(__file__), evidence=evidence,
        changes=["Only author fields change; raw Crossref remains untouched.",
                 "Use top-level article-meta contributor order, not consortium membership expansion.",
                 "2020 collective converted from given/family to CSL literal.",
                 "2022 combined collective label preserved verbatim apart from whitespace.",
                 "Personal names are not deduplicated or spelling-normalized."],
        limitations=["Not a complete manuscript bibliography or journal-style rendering.",
                     "XML snapshots are retained locally, not redistributed in this export."]),
        indent=2, ensure_ascii=True) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("raw", "xml-directory", "output", "provenance"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    export(args.raw, args.xml_directory, args.output, args.provenance)
