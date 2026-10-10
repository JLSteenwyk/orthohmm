"""Insert a pinned completed-results supplement without changing parent bytes."""

import json
from pathlib import Path

from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.export_native_factorial_progress import require


def insert(parent, supplement, anchor):
    require(parent.count(anchor) == 1 and supplement.strip(), "Missing/ambiguous insertion")
    piece = supplement.rstrip("\n") + "\n\n"
    require(piece not in parent, "Supplement already integrated")
    revised = parent.replace(anchor, piece + anchor, 1)
    require(revised.count(piece) == 1 and revised.replace(piece, "", 1) == parent,
            "Unscoped parent alteration")
    return revised, piece


def run(parent_ref, supplement_ref, output, receipt, anchor):
    output, receipt = Path(output).absolute(), Path(receipt).absolute()
    require(output != receipt and output.suffix == ".md" and receipt.suffix == ".json",
            "Distinct manuscript/record paths required")
    require(output.parent.resolve() == receipt.parent.resolve() == Path(parent_ref["path"]).resolve().parent
            == Path(supplement_ref["path"]).resolve().parent, "Relative links would change")
    require(not any(p.exists() or p.is_symlink() for p in (output, receipt)), "Outputs must be fresh")
    for ref in (parent_ref, supplement_ref):
        check(ref)
    parent = Path(parent_ref["path"]).read_bytes()
    supplement = Path(supplement_ref["path"]).read_bytes()
    revised, piece = insert(parent.decode("utf-8"), supplement.decode("utf-8"), anchor)
    data = revised.encode("utf-8")
    require(data.replace(piece.encode("utf-8"), b"", 1) == parent, "Parent bytes changed")
    for ref in (parent_ref, supplement_ref):
        check(ref)
    with output.open("xb") as stream:
        stream.write(data)
    result = dict(schema="pinned_manuscript_supplement_v1", source=record(__file__),
                  inputs=dict(parent=parent_ref, supplement=supplement_ref), manuscript=record(output),
                  inserted_section=piece, parent_unchanged_except_insertion=True,
                  new_inference_or_scoring=False, publication_ready=False, manuscript_rendered=False)
    with receipt.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return result
