"""Export explicitly selected Crossref citation metadata as CSL-JSON."""

import argparse
from datetime import datetime, timezone
from html.parser import HTMLParser
import json
from pathlib import Path
import sys
from urllib.parse import quote
from urllib.request import Request, urlopen

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record


class PlainText(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.parts = []

    def handle_data(self, data):
        self.parts.append(data)


def citation(metadata, item):
    if metadata["DOI"].casefold() != item["doi"].casefold():
        raise ValueError("Returned DOI differs from requested DOI")
    types = {"journal-article": "article-journal", "posted-content": "article"}
    if metadata["type"] not in types or not isinstance(metadata.get("title"), str):
        raise ValueError("Unsupported record type or title")
    title = PlainText()
    title.feed(metadata["title"])
    title.close()
    authors = []
    for author in metadata.get("author", []):
        name = {k: author[k] for k in ("given", "family", "suffix", "literal") if author.get(k)}
        if not name.get("family") and not name.get("literal"):
            if author.get("name"):
                name = {"literal": author["name"]}
            else:
                raise ValueError("Incomplete author name")
        authors.append(name)
    issued = metadata["issued"]
    if not authors or not "".join(title.parts).strip() or issued["date-parts"][0][0] != item["year"]:
        raise ValueError("Missing citation fields or publication year disagreement")
    result = {"id": item["id"], "type": types[metadata["type"]], "title": "".join(title.parts),
              "DOI": item["doi"], "URL": "https://doi.org/" + item["doi"],
              "author": authors, "issued": issued}
    for field in ("container-title", "volume", "issue", "page", "publisher"):
        value = metadata.get(field)
        if value:
            if not isinstance(value, str):
                raise ValueError("Unexpected bibliographic field type: " + field)
            result[field] = value
    if metadata["type"] == "posted-content":
        result["genre"] = "preprint"
    article_number = metadata.get("article-number")
    if article_number:
        if not isinstance(article_number, str):
            raise ValueError("Unexpected article identifier type")
        result["number"] = article_number
    return result


def export(manifest, raw_directory, output, provenance, cached_provenance=None):
    if output.exists() or provenance.exists():
        raise FileExistsError("Citation output or provenance already exists")
    selected = json.loads(manifest.read_text())
    if not selected or len({r["id"] for r in selected}) != len(selected) or len({r["doi"].casefold() for r in selected}) != len(selected):
        raise ValueError("Empty or duplicate citation selection")
    if any(Path(r["id"]).name != r["id"] or r["id"] in (".", "..") for r in selected):
        raise ValueError("Unsafe citation identifier")
    cached = None
    if cached_provenance is not None:
        cached = json.loads(cached_provenance.read_text())
        if cached["selection"]["sha256"] != record(manifest)["sha256"]:
            raise ValueError("Cached selection checksum differs")
        if [(r["id"], r["doi"]) for r in cached["records"]] != [(r["id"], r["doi"]) for r in selected]:
            raise ValueError("Cached source inventory differs")
    else:
        raw_directory.mkdir(parents=True, exist_ok=False)
    citations, sources = [], []
    for index, item in enumerate(selected):
        url = "https://api.crossref.org/works/" + quote(item["doi"], safe="") + "/transform/application/vnd.citationstyles.csl+json"
        raw = raw_directory / (item["id"] + ".json")
        if cached is not None:
            source = cached["records"][index]
            observed = record(raw)
            if source["url"] != url or any(observed[k] != source["response"][k] for k in ("sha256", "bytes")):
                raise ValueError("Cached response provenance differs")
            data = raw.read_bytes()
            retrieved = source["retrieved_utc"]
        else:
            request = Request(url, headers={"User-Agent": "OrthoHMM-publication-reference-audit/1.0"})
            with urlopen(request, timeout=30) as response:
                data = response.read()
            raw.write_bytes(data)
            retrieved = datetime.now(timezone.utc).isoformat()
        entry = citation(json.loads(data), item)
        citations.append(entry)
        sources.append({"id": item["id"], "doi": item["doi"], "url": url, "response": record(raw),
                        "retrieved_utc": retrieved, "authors": len(entry["author"])})
    with output.open("x") as stream:
        json.dump(citations, stream, indent=2, ensure_ascii=True)
        stream.write("\n")
    result = {"status": "selected_citation_metadata_exported", "selection": record(manifest), "source": record(__file__),
              "output": record(output), "records": sources,
              "article_identifier_mapping": "Crossref article-number is CSL number, not page.",
              "limitations": ["Crossref bibliographic metadata, not full-text scientific review or software-version provenance.",
                              "Author names reflect deposited metadata; known suffix/name omissions need explicit review, not silent invention.",
                              "Only selected references exported; dependency/resource bibliography and journal-specific style remain separate.",
                              "Raw responses retained locally; abstracts and cited-reference lists are omitted from committed CSL-JSON."]}
    with provenance.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("manifest", "raw-directory", "output", "provenance"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--cached-provenance", type=Path, help="Replay checksum-verified raw responses without network access")
    args = parser.parse_args()
    result = export(args.manifest, args.raw_directory, args.output, args.provenance, args.cached_provenance)
    print(f"Exported {len(result['records'])} citation records")
