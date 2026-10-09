"""Recover one clipped manuscript table without changing scientific HTML content."""

import argparse
from html.parser import HTMLParser
import json
from pathlib import Path

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


PINS = {
    "html": ("native_qfo_four_cell_strata_manuscript_20261009_v1_review.html",
             "d457907836a1d739bb1fcaeee5751c6a942f244f355004d2ff786277a613d1c8"),
    "assets": ("native_qfo_four_cell_strata_manuscript_20261009_v1_assets.json",
               "d430d54836357d2790f742102f78a9e31bc5d940a80ebacebe6303f7182418d2"),
}
ANCHOR = "profile-refinement-across-all-fixed-bins"
SELECTOR = "#" + ANCHOR + " + p + table"
STYLE = '\n<style id="profile-stratum-print-recovery">\n@media print {\n' + "\n".join([
    SELECTOR + " { display: table; width: 100%; table-layout: fixed; overflow: visible; font-size: 10pt; line-height: 1.25; }",
    SELECTOR + " col:nth-child(1) { width: 12% !important; }",
    SELECTOR + " col:nth-child(2) { width: 40% !important; }",
    SELECTOR + " col:nth-child(3) { width: 9% !important; }",
    SELECTOR + " col:nth-child(n+4) { width: 13% !important; }",
    SELECTOR + " th, " + SELECTOR + " td { overflow-wrap: anywhere; padding: 4px; vertical-align: top; }",
    SELECTOR + " td:nth-child(n+3) { white-space: nowrap; }",
    SELECTOR + " tr { break-inside: avoid; }",
]) + "\n}\n</style>\n"


def require(condition, message):
    if not condition:
        raise ValueError(message)


class ProfileTable(HTMLParser):
    def __init__(self):
        super().__init__()
        self.anchors = 0
        self.state = "before"
        self.rows = []
        self.row = []
        self.cell = None

    def handle_starttag(self, tag, attrs):
        if dict(attrs).get("id") == ANCHOR:
            require(tag == "h3", "Wrong profile anchor element")
            self.anchors += 1
            self.state = "heading"
        elif self.state == "after_heading":
            require(tag == "p", "Profile table is not after an introductory paragraph")
            self.state = "paragraph"
        elif self.state == "after_paragraph":
            require(tag == "table", "Scoped selector does not match profile table")
            self.state = "table"
        elif self.state == "table" and tag in ("th", "td"):
            require(self.cell is None, "Nested table cell")
            self.cell = []

    def handle_data(self, text):
        if self.cell is not None:
            self.cell.append(text)

    def handle_endtag(self, tag):
        if self.state == "heading" and tag == "h3": self.state = "after_heading"
        elif self.state == "paragraph" and tag == "p": self.state = "after_paragraph"
        elif self.state == "table":
            if tag in ("th", "td"):
                require(self.cell is not None, "Missing table cell start")
                self.row.append(" ".join("".join(self.cell).split()))
                self.cell = None
            elif tag == "tr":
                self.rows.append(self.row)
                self.row = []
            elif tag == "table": self.state = "done"


def styled(html):
    parser = ProfileTable()
    parser.feed(html)
    parser.close()
    require(parser.anchors == 1 and parser.state == "done" and len(parser.rows) == 24
            and parser.rows[0] == ["Suite", "Bin", "Families", "F1 Difference", "PPV Difference", "TPR Difference"]
            and all(len(row) == 6 for row in parser.rows), "Incomplete profile table DOM")
    require(html.count("</head>") == 1 and "profile-stratum-print-recovery" not in html, "Already styled or ambiguous head")
    updated = html.replace("</head>", STYLE + "</head>", 1)
    require(updated.replace(STYLE, "", 1) == html, "Changed scientific HTML content")
    return updated


def run(root, output, assets):
    base = Path(root).resolve() / "benchmark_tools/results"
    output, assets = Path(output).absolute(), Path(assets).absolute()
    require(output != assets and output.parent.resolve() == assets.parent.resolve() == base, "Preserve relative assets")
    for path in (output, assets):
        if path.exists() or path.is_symlink():
            raise FileExistsError(path)
    refs = {key:record(base / name) for key, (name, _) in PINS.items()}
    require(all(refs[key]["sha256"] == digest for key, (_, digest) in PINS.items()), "Changed original render input")
    receipt = json.loads(Path(refs["assets"]["path"]).read_text())
    require(receipt["status"] == "manuscript_review_rendered" and receipt["html"] == refs["html"]
            and receipt["publication_ready"] is False, "Wrong original render binding")
    source = record(__file__)
    checked = [*refs.values(), *receipt["sources"], *receipt["targets"], source]
    for ref in checked:
        check(ref)
    html = styled(Path(refs["html"]["path"]).read_text())
    with output.open("x") as stream:
        stream.write(html)
    revised = dict(receipt, html=record(output), sources=[*receipt["sources"], refs["assets"], refs["html"], source],
                   render_command_scope="inherited executed Pandoc invocation; recovery only inserts scoped print CSS",
                   rendering_recovery=dict(operation="insert_scoped_print_css_only", original_inputs=refs,
                                           source=source, scientific_html_unchanged=True, profile_rows=23,
                                           all_columns_required=6, selector=SELECTOR))
    for ref in checked:
        check(ref)
    with assets.open("x") as stream:
        json.dump(revised, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return revised


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "output", "assets"):
        parser.add_argument("--" + name, required=True, type=Path)
    args = parser.parse_args()
    result = run(args.root, args.output, args.assets)
    print(json.dumps(dict(html=result["html"], rendering_recovery=result["rendering_recovery"])))
