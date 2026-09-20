"""Inventory explicit manuscript citations, not the adequacy of scientific attribution."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess
from urllib.parse import unquote, urlsplit

from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.render_publication_bibliography import text


def doi_url(url):
    parsed = urlsplit(url)
    if parsed.scheme.lower() in {"http", "https"} and parsed.netloc.lower() in {"doi.org", "dx.doi.org"}:
        return unquote(parsed.path).lstrip("/").lower()
    return None


def inventory(document, bibliography):
    by_id, by_doi, by_url = {}, {}, {}
    for row in bibliography:
        identity = row["id"]
        if not isinstance(identity, str) or not identity or identity in by_id:
            raise ValueError("Missing or duplicate bibliography identity")
        by_id[identity] = row
        if row.get("DOI"):
            doi = row["DOI"].lower()
            if doi in by_doi:
                raise ValueError("Ambiguous bibliography DOI")
            by_doi[doi] = identity
        if row.get("URL"):
            by_url.setdefault(row["URL"], []).append(identity)
    occurrences, external, internal, used = [], [], 0, set()

    def visit(node, section):
        nonlocal internal
        if isinstance(node, list):
            for child in node:
                visit(child, section)
        elif isinstance(node, dict):
            if node.get("t") == "Link":
                _, label, target = node["c"]
                url = target[0]
                doi = doi_url(url)
                matches = [by_doi[doi]] if doi in by_doi else by_url.get(url, [])
                row = dict(section=section, label=text(label), url=url, doi=doi, bibliography_ids=matches)
                if doi is not None or matches:
                    occurrences.append(dict(row, kind="link", status="matched" if len(matches) == 1 else "ambiguous" if matches else "unmatched_doi"))
                    if len(matches) == 1:
                        used.add(matches[0])
                elif urlsplit(url).scheme in {"http", "https"}:
                    external.append(row)
                else:
                    internal += 1
            if node.get("t") == "Cite":
                for citation in node["c"][0]:
                    identity = citation["citationId"]
                    matched = identity in by_id
                    occurrences.append(dict(section=section, kind="cite", citation_id=identity,
                        bibliography_ids=[identity] if matched else [], status="matched" if matched else "unknown_citation_id"))
                    if matched:
                        used.add(identity)
            for child in node.values():
                visit(child, section)

    section = "Preamble"
    for block in document["blocks"]:
        if block.get("t") == "Header":
            section = text(block["c"][2])
        visit(block, section)
    return dict(status="explicit_citation_inventory", bibliography_entries=len(by_id),
        matched_bibliography_ids=sorted(used), bibliography_not_explicitly_matched=sorted(set(by_id)-used),
        occurrences=occurrences, unresolved_citations=[r for r in occurrences if r["status"] != "matched"],
        other_external_links=external, internal_link_occurrences=internal,
        citation_adequacy_established=False, publication_ready=False,
        limitations=["Matches explicit DOI resolver links, exact bibliography URLs and Pandoc citation IDs only.",
            "Other external links may be citations requiring manual review; unmatched bibliography entries are not automatically irrelevant.",
            "Does not detect uncited claims or validate attribution, source semantics, versions, rights or linked local evidence.",
            "Section names locate occurrences; this does not supply exact Markdown line numbers."])


def audit(manuscript, bibliography):
    binary = shutil.which("pandoc")
    if binary is None:
        raise FileNotFoundError("pandoc")
    sources = [record(p) for p in (manuscript, bibliography, binary, __file__)]
    command = [binary, "--from=markdown", "--to=json", str(manuscript)]
    parsed = subprocess.run(command, capture_output=True, text=True, check=True, timeout=60)
    result = inventory(json.loads(parsed.stdout), json.loads(Path(bibliography).read_text()))
    for item in sources:
        check(item)
    return dict(result, sources=sources, command=command, stderr=parsed.stderr,
        pandoc_version=subprocess.check_output([binary, "--version"], text=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("manuscript", "bibliography", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.manuscript, args.bibliography)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
