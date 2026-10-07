"""Append the already-reviewed native figure to a new pinned presentation inventory."""

import argparse
import copy
import json
from pathlib import Path

from benchmark_tools.assemble_publication_review import pinned, record


ROOT = Path(__file__).resolve().parents[1]
BASE = "benchmark_tools/results/"
PRIOR = dict(path=BASE + "publication_figure_selection_20261004_v3.json", bytes=0,
             sha256="a9ece39a341947250dfe662ef9649c3abd8289342d5e7c4fbdb7d5cc7b5bfdf2")
MAIN = dict(path=BASE + "publication_main_print_20261006_v2/document.pdf", bytes=349569,
            sha256="6f288300c73352e26cb485cc22b870d583510be1c660ea5f6f674e574711b872")
PRINT = dict(path=BASE + "publication_main_print_20261006_v2/print.json", bytes=0,
             sha256="3342e5cf442508ae45967e3edef726dd6eb30e7e9888d841ab46705ed6e7640d")
FIGURE = dict(path=BASE + "native_qfo_three_cell_figure_20261006_v1/native_qfo_three_cell.pdf", bytes=24109,
              sha256="c3e6ca725b1a0e48653179956f87d41c8f957093517d9610324843387f8613f1")


def checked(root, ref):
    observed = record(root / ref["path"])
    if observed["sha256"] != ref["sha256"] or ref["bytes"] and observed["bytes"] != ref["bytes"]:
        raise ValueError("Presentation input checksum mismatch")
    return observed


def select(previous, printed):
    if (previous["schema"] != "publication_figure_selection_v1"
            or previous["publication_ready"] is not False
            or [f["number"] for f in previous["figures"]] != list(range(1, 17))
            or printed["status"] != "verified_html_printed"
            or printed["publication_ready"] is not False or printed["page_count"] != 19
            or (printed["pdf"]["bytes"], printed["pdf"]["sha256"]) != (MAIN["bytes"], MAIN["sha256"])):
        raise ValueError("Unexpected historical inventory or new print receipt")
    result = copy.deepcopy(previous)
    result.update(main=MAIN, main_pages=printed["page_count"], publication_ready=False,
                  scientific_settings_or_results_changed=False)
    result["figures"].append(dict(number=17, title="Three admitted native QfO cells",
        caption="Initial sensitive HMM search is on in all three P0 cells. F1, GO/EC/FAS, prediction coverage "
            "and SwissTrees precision/recall remain separate; two conditional adjusted F1 intervals include zero. "
            "Four fresh cells and twelve contrasts remain unavailable. Recovered R1 timing remains failed; "
            "this is not initial-HMM-off, complete-factorial or isolated-efficiency evidence.",
        pdf=FIGURE, link_aliases=[FIGURE["path"]]))
    return result


def build(output, root=ROOT):
    output, root = Path(output).absolute(), Path(root).resolve()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    refs = [checked(root, r) for r in (PRIOR, PRINT, MAIN, FIGURE)]
    previous, printed = [json.loads((root / r["path"]).read_text()) for r in (PRIOR, PRINT)]
    result = select(previous, printed)
    result["provenance"] = [dict(path=r["path"], bytes=observed["bytes"], sha256=observed["sha256"])
                            for r, observed in zip((PRIOR, PRINT, MAIN, FIGURE), refs)]
    source = record(__file__)
    result["provenance"].append(dict(source, path=str(Path(__file__).relative_to(ROOT))))
    for item in result["figures"]:
        pinned(root, item["pdf"])
    for ref in refs:
        observed = record(ref["path"])
        if observed != ref:
            raise ValueError("Input changed during selection")
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return record(output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(build(args.output)))
