"""Inventory published family-level duplication summaries, without outcome joins."""

import argparse
from html.parser import HTMLParser
import json
from pathlib import Path

from benchmark_tools.prepare_ob_candidate_neighborhood import record, check

URL = "https://swisstree.sib.swiss/cgi-bin/swisst?page=gold_standard"
FAMILIES = dict(zip([f"ST{i:03d}" for i in range(1, 20)],
    ["POP", "NOX", "VATB", "SERC", "SUMF", "HOX", "RPS", "BAMBI", "ASTER", "CITE",
     "GH14", None, "TRFE", "CASP", "BAR", "PSEN", "Clusterin", "APP", "MAPT"]))


class TableParser(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.rows, self.row, self.cell = [], None, None

    def handle_starttag(self, tag, attrs):
        if tag == "tr":
            self.row = []
        elif tag in {"td", "th"} and self.row is not None:
            self.cell = []

    def handle_data(self, data):
        if self.cell is not None:
            self.cell.append(data)

    def handle_endtag(self, tag):
        if tag in {"td", "th"} and self.cell is not None:
            self.row.append(" ".join("".join(self.cell).split()))
            self.cell = None
        elif tag == "tr" and self.row is not None:
            self.rows.append(self.row)
            self.row = None


def parse(text):
    parser = TableParser()
    parser.feed(text)
    headers = [r for r in parser.rows if r and r[0] == "ID" and "Duplications" in r]
    if len(headers) != 1 or len(set(headers[0])) != len(headers[0]):
        raise ValueError("Missing or ambiguous duplication table header")
    header, output = headers[0], {}
    for cells in parser.rows:
        if not cells or not cells[0].startswith("ST"):
            continue
        if cells[0] not in FAMILIES or cells[0] in output or len(cells) != len(header):
            raise ValueError("Unexpected or duplicate SwissTree family row")
        row = dict(zip(header, cells))
        if not row["Duplications"].isascii() or not row["Duplications"].isdigit():
            raise ValueError("Missing or nonnumeric published duplication count")
        output[cells[0]] = dict(source_fields=row, qfo_family=FAMILIES[cells[0]],
            published_duplications=int(row["Duplications"]),
            updated_marker=row["Family"].endswith("*"))
    if set(output) != set(FAMILIES):
        raise ValueError("Incomplete published SwissTree family inventory")
    return output


def inventory(source, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    identity = record(source)
    rows = parse(source.read_text())
    check(identity)
    report = dict(status="published_swisstree_duplication_summary_inventory", url=URL,
        raw_source=identity, source=record(__file__), families=rows, prediction_statistics_evaluated=False,
        exact_benchmark_duplication_history_validated=False, publication_ready=False,
        limitations=["Current published summary counts may use different release/topology/taxon coverage than QfO 2020.",
            "Family-name crosswalk is explicit; exact tree and leaf correspondence has not been validated.",
            "ST012 ant transformer family is absent from the retained 18-family benchmark panel.",
            "Counts are curated-resource summaries, not experimental ground truth or independent test data.",
            "No outcome strata, normalization by tree size, or accuracy comparisons computed."])
    with output.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    inventory(args.source, args.output)
