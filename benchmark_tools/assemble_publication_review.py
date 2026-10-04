"""Assemble pinned manuscript/figure PDFs without rerunning scientific analyses."""

import argparse
import hashlib
import json
from pathlib import Path

import fitz


def require(condition, message):
    if not condition:
        raise ValueError(message)


def record(path):
    path = Path(path).resolve()
    content = path.read_bytes()
    return {"path": str(path), "bytes": len(content),
            "sha256": hashlib.sha256(content).hexdigest()}


def pinned(root, reference):
    relative = Path(reference["path"])
    require(not relative.is_absolute() and ".." not in relative.parts,
            "Unsafe input path")
    path = (root / relative).resolve()
    require(path.is_relative_to(root), "Input escapes root")
    observed = record(path)
    require((observed["bytes"], observed["sha256"]) ==
            (reference["bytes"], reference["sha256"]),
            "Input checksum mismatch: " + str(relative))
    return observed


def box(page, rectangle, text, size, font="helv"):
    require(page.insert_textbox(rectangle, text, fontsize=size, fontname=font) >= 0,
            "Guide text does not fit: " + text)


def assemble(root, selection, selection_sha256, output):
    root, selection, output = map(lambda path: Path(path).resolve(),
                                  (root, selection, output))
    if output.exists():
        raise FileExistsError("Refusing existing output directory: " + str(output))
    selection_ref = record(selection)
    require(selection_ref["sha256"] == selection_sha256,
            "Selection checksum mismatch")
    chosen = json.loads(selection.read_text())
    require(chosen["schema"] == "publication_figure_selection_v1", "Unknown selection schema")
    main_ref = pinned(root, chosen["main"])
    inputs = [selection_ref, main_ref, record(__file__)]
    inputs.extend(pinned(root, ref) for ref in chosen["provenance"])
    figures = []
    aliases = {}
    for number, item in enumerate(chosen["figures"], 1):
        require(item["number"] == number, "Nonsequential figure inventory")
        require(all(isinstance(item[key], str) and item[key].strip()
                    for key in ("title", "caption")), "Missing title or caption")
        ref = pinned(root, item["pdf"])
        require(ref["path"] not in [figure["pdf"]["path"] for figure in figures],
                "Duplicate figure PDF")
        inputs.append(ref)
        figure = {"number": number, "title": item["title"], "caption": item["caption"],
                  "pdf": ref, "relative_path": item["pdf"]["path"]}
        figures.append(figure)
        for alias in [ref["path"], *item.get("link_aliases", [])]:
            require(isinstance(alias, str) and alias, "Invalid link alias")
            require(alias not in aliases or aliases[alias] == number, "Ambiguous link alias")
            aliases[alias] = number
    require(figures, "Empty figure inventory")
    # Validate every PDF before creating the output directory.
    fitz.TOOLS.mupdf_warnings(reset=True)
    for ref in [main_ref, *[item["pdf"] for item in figures]]:
        with fitz.open(ref["path"]) as source:
            expected = chosen["main_pages"] if ref == main_ref else 1
            require(type(expected) is int and expected > 0 and len(source) == expected,
                    "Unexpected source page count")
            require(not source.is_repaired, "Repaired source PDF")
    document = fitz.open()
    with fitz.open(main_ref["path"]) as source:
        document.insert_pdf(source)
    main_pages = len(document)
    guide_pages = (len(figures) + 3) // 4
    first_figure = main_pages + guide_pages
    guide_links = []
    for group in range(guide_pages):
        page = document.new_page(width=595.276, height=841.89)
        box(page, fitz.Rect(40, 38, 555, 72), "Figure Appendix - Working Review", 18, "hebo")
        box(page, fitz.Rect(40, 78, 555, 120),
            "Main text: pages 1-%d. Original vector pages and dimensions are preserved. "
            "Shared-host timings have unknown contention effects; no isolated speed ranking "
            "or publication-readiness claim is made." % main_pages, 10)
        for slot, item in enumerate(figures[group * 4:(group + 1) * 4]):
            top = 135 + slot * 160
            box(page, fitz.Rect(40, top, 555, top + 30),
                "A%d. %s" % (item["number"], item["title"]), 12, "hebo")
            box(page, fitz.Rect(40, top + 35, 555, top + 110), item["caption"], 10.5)
            box(page, fitz.Rect(40, top + 115, 555, top + 145),
                "Source: " + item["relative_path"], 8)
            guide_links.append((page.number, fitz.Rect(40, top, 555, top + 145),
                                first_figure + item["number"] - 1))
        box(page, fitz.Rect(40, 790, 555, 817),
            "Guide %d of %d - click an entry or use PDF bookmarks." % (group + 1, guide_pages), 9)
    for item in figures:
        item["combined_page"] = len(document) + 1
        with fitz.open(item["pdf"]["path"]) as source:
            document.insert_pdf(source)
    for page, rectangle, target in guide_links:
        document[page].insert_link({"kind": fitz.LINK_GOTO, "from": rectangle, "page": target})
    destinations = {alias: first_figure + number - 1 for alias, number in aliases.items()}
    redirected, restored = [], []
    # insert_pdf can reinterpret Chrome file-URI actions as launch actions.
    # Use the original link actions, not guessed path substitutions.
    with fitz.open(main_ref["path"]) as source:
        for index in range(main_pages):
            page = document[index]
            for original in source[index].get_links():
                filename = original.get("file")
                if filename not in destinations and original["kind"] != fitz.LINK_LAUNCH:
                    continue
                copied = [link for link in page.get_links() if link["from"] == original["from"]]
                require(len(copied) == 1, "Ambiguous copied link")
                action = {"xref": copied[0]["xref"], "from": original["from"]}
                if filename in destinations:
                    action.update(kind=fitz.LINK_GOTO, page=destinations[filename])
                    redirected.append({"main_page": index + 1, "source": filename,
                                       "combined_page": destinations[filename] + 1})
                else:
                    uri_type, uri = source.xref_get_key(original["xref"], "A/URI")
                    require(uri_type == "string" and uri.startswith("file:"),
                            "Unexpected original launch action")
                    action.update(kind=fitz.LINK_URI, uri=uri)
                    restored.append({"main_page": index + 1, "uri": uri})
                page.update_link(action)
    document.set_toc([[1, "Main Text", 1], [1, "Figure Guide", main_pages + 1],
                      *[[1, "A%d. %s" % (item["number"], item["title"]), item["combined_page"]]
                        for item in figures]])
    document.set_metadata({"title": "OrthoHMM Working Manuscript With Figure Appendix",
                           "subject": "Shared-host observations; scientific limitations retained"})
    require(not fitz.TOOLS.mupdf_warnings(reset=True), "PDF backend warnings")
    output.mkdir(parents=True)
    output_pdf = output / "document.pdf"
    document.save(output_pdf, garbage=3, deflate=True)
    document.close()
    checks, renders = [], []
    with fitz.open(output_pdf) as combined:
        require(len(combined) == main_pages + guide_pages + len(figures), "Wrong combined page count")
        source_pages = [(main_ref, page, page) for page in range(main_pages)] + [
            (item["pdf"], 0, item["combined_page"] - 1) for item in figures]
        for ref, original_index, combined_index in source_pages:
            with fitz.open(ref["path"]) as original:
                before, after = original[original_index], combined[combined_index]
                require(before.rect == after.rect and before.get_text("words") == after.get_text("words"),
                        "Source text or geometry changed")
                a, b = before.get_pixmap(alpha=False), after.get_pixmap(alpha=False)
                require((a.width, a.height, a.n, a.samples) == (b.width, b.height, b.n, b.samples),
                        "Source pixels changed")
                checks.append({"source": ref, "original_page": original_index + 1,
                               "combined_page": combined_index + 1,
                               "pixel_sha256": hashlib.sha256(a.samples).hexdigest(),
                               "exact_text_geometry_and_pixels": True})
        for index in range(main_pages, len(combined)):
            path = output / ("page_%02d.png" % (index + 1))
            combined[index].get_pixmap(matrix=fitz.Matrix(1.25, 1.25), alpha=False).save(path)
            renders.append(record(path))
        expected = [(index + 1, rectangle, target + 1)
                    for index, rectangle, target in guide_links]
        observed = [(index + 1, link["from"], link["page"] + 1)
                    for index in range(main_pages, first_figure)
                    for link in combined[index].get_links() if link["kind"] == fitz.LINK_GOTO]
        require(observed == expected, "Guide link targets changed")
        for index in range(main_pages):
            actual = [link for link in combined[index].get_links() if link["kind"] == fitz.LINK_GOTO]
            require(sorted(link["page"] + 1 for link in actual) ==
                    sorted(link["combined_page"] for link in redirected if link["main_page"] == index + 1),
                    "Manuscript figure link targets changed")
        outlines = combined.get_toc()
        require(len(outlines) == len(figures) + 2, "Wrong bookmarks")
    for ref in inputs:
        require(record(ref["path"]) == ref, "Input changed during assembly")
    warnings = fitz.TOOLS.mupdf_warnings(reset=True)
    require(not warnings, warnings)
    report = {"schema": "publication_figure_assembly_v1", "selection": selection_ref,
              "pdf": record(output_pdf), "inputs": inputs, "figures": figures,
              "main_pages": main_pages, "guide_pages": guide_pages, "figure_pages": len(figures),
              "total_pages": main_pages + guide_pages + len(figures), "bookmarks": outlines,
              "redirected_main_figure_links": redirected, "restored_original_file_uri_actions": restored,
              "preserved_source_pages": checks, "rendered_pages": renders,
              "original_text_and_figure_pages_pixel_identical": True,
              "pymupdf_version": fitz.__version__, "mupdf_warnings": warnings,
              "visual_review_complete": False, "scientific_settings_or_results_changed": False,
              "publication_ready": False,
              "limitations": ["Presentation only; scientific inference and statistical inputs are not replayed.",
                              "External nonfigure file links retain historical locations, not portable data access.",
                              "Shared-host resource evidence does not establish isolated tool speed.",
                              "Source page sizes differ; journal typesetting and publication readiness are not certified."]}
    with (output / "assembly.json").open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--selection", type=Path, required=True)
    parser.add_argument("--selection-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = assemble(args.root, args.selection, args.selection_sha256, args.output)
    print(json.dumps({"pdf": report["pdf"], "pages": report["total_pages"],
                      "preserved_pages": len(report["preserved_source_pages"]),
                      "redirected_links": len(report["redirected_main_figure_links"])}, indent=2))


if __name__ == "__main__":
    main()
