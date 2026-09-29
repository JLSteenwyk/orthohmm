"""Check PDF page bounds and render selected pages without asserting visual review."""

import argparse
import json
from pathlib import Path

import fitz

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


def review(pdf, assets, output, phrases):
    if output.exists():
        raise FileExistsError(output)
    if not phrases or any(not p.strip() for p in phrases):
        raise ValueError("Require explicit nonempty text selectors")
    refs = [record(pdf.resolve()), record(assets.resolve()), record(Path(__file__).resolve())]
    receipt = json.loads(assets.read_text())
    if receipt.get("status") != "manuscript_review_rendered":
        raise ValueError("Require successful HTML render receipt")
    refs.extend([receipt["html"], *receipt["sources"], *receipt["targets"]])
    for ref in refs:
        check(ref)
    violations, selected = [], set()
    matches = {p: [] for p in phrases}
    with fitz.open(pdf) as document:
        for i, page in enumerate(document):
            text = " ".join(page.get_text().split())
            for phrase in phrases:
                if " ".join(phrase.split()) in text:
                    matches[phrase].append(i + 1)
                    selected.update(range(max(0, i - 1), min(len(document), i + 2)))
            for block in page.get_text("dict")["blocks"]:
                box = fitz.Rect(block["bbox"])
                bounds = page.rect + (-1, -1, 1, 1)
                if not bounds.contains(box):
                    violations.append(dict(page=i + 1, type=block["type"], bounds=list(box)))
        if any(not pages for pages in matches.values()):
            raise ValueError("One or more required text selectors are absent")
        output.mkdir(parents=True, exist_ok=False)
        images = []
        for i in sorted(selected):
            path = output / f"page_{i + 1:03d}.png"
            document[i].get_pixmap(matrix=fitz.Matrix(1.5, 1.5)).save(path)
            images.append(record(path.resolve()))
        result = dict(status="pdf_bounds_checked", page_count=len(document), bounds_tolerance_points=1,
            bounds_violations=violations, phrase_matches=matches, rendered_pages=images,
            checked_records=refs, pymupdf_version=fitz.__version__,
            local_occurrences=receipt["local_occurrences"], unique_targets=receipt["unique_targets"],
            untracked_targets=receipt["untracked_targets"], visual_review_complete=False,
            publication_ready=False, limitations=[
                "Page bounds do not establish readability, nonoverlap, image fidelity or scientific correctness.",
                "Text selectors and adjacent pages guide manual review, not full-document visual certification.",
                "PDF correspondence to HTML is not proven solely by these post-print file checks."])
    for ref in refs:
        check(ref)
    (output / "report.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("pdf", "assets", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--phrase", action="append", required=True)
    args = parser.parse_args()
    result = review(args.pdf, args.assets, args.output, args.phrase)
    print(json.dumps({k: result[k] for k in ("page_count", "bounds_violations", "phrase_matches")}))
