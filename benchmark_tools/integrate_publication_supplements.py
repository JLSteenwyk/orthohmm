"""Append completed, bounded results without altering prior manuscript content."""

import json
import math
from pathlib import Path

from benchmark_tools.export_native_qfo_three_cell_strata import record, require


PINS = {
    "parent": ("controlled_fragment_trace_manuscript_20261010_v1.md", "df932f2299f47494851235fce68989a33f2bc24390bb4511597f445239b24b3c"),
    "poisson": ("paired_poisson_f1_validation_20261010_v1.json", "314ddd68b690fb8a8ef210f49fa5be85e3a1a6a17b8f6362cefbf9d091ebee07"),
    "relocation": ("controlled_fragment_relocated_readback_20261010_v2.json", "0128cd57399dbb4d1bdb7b72660d3be046b098f802059b81eb30162e6d25257f"),
}
ANCHORS = ("### Controlled Fragment Scope And Limitations\n", "## References\n")


def sections(docs):
    poisson, relocation = docs["poisson"], docs["relocation"]
    require(poisson["status"] == "conditional_exact_boundary_enumeration_complete"
            and poisson["all_conditional_cells_pass"] is True, "Incomplete conditional panel")
    require(all(poisson[k] is False for k in (
        "native_intervals_admitted", "overall_uncertainty_method_admitted",
        "publication_ready", "independent_biological_confirmation", "new_inference_or_benchmark_scoring"))
        and poisson["new_bootstrap_draws"] == 0, "Changed conditional scope")
    rows, controls = poisson["rows"], poisson["model_violation_controls"]
    means = ((0, 0, 0), (2, 0, 0), (1, 1, 0), (1, 0, 1), (.5, .3, .7), (10, 20, 12))
    require(len(rows) == 12 and {(r["truth_mass"], tuple(r["means"])) for r in rows}
            == {(t, m) for t in (512, 23934) for m in means}, "Changed conditional cells")
    require(all(r["numerical_check_pass"] is True and r["numerical_rounding_certified"] is False
                and type(r["covered_enumerated_mass"]) in (int, float)
                and math.isfinite(r["covered_enumerated_mass"])
                and .95 <= r["covered_enumerated_mass"] <= 1
                and type(r["omitted_tail_mass"]) in (int, float)
                and math.isfinite(r["omitted_tail_mass"])
                and 0 <= r["omitted_tail_mass"] <= 1e-12 for r in rows), "Invalid panel numbers")
    require(len(controls) == 2 and {c["truth_mass"] for c in controls} == {512, 23934}
            and all(c["coverage"] == 0 and c["poisson_marginal_assumption_satisfied"] is False
                    and c["native_dependence_model"] is False for c in controls), "Changed failure controls")
    require(relocation["status"] == "relocated_native_readback_matches_retained_result"
            and relocation["artifact_files"] == 1783 and relocation["artifact_payload_bytes"] == 39428185,
            "Changed relocation result")
    require(all(relocation[k] is False for k in (
        "historical_inference_reproduced", "new_inference_or_scoring", "original_artifacts_modified",
        "original_path_fallback", "original_reader_kernels_changed", "publication_ready")), "Changed replay scope")
    native = relocation["native_readback"]
    require(native["status"] == "independent_native_stages_and_na_table_verified"
            and (native["stage_rows"], native["selected_cases"], native["native_contexts"], native["checked_input_files"])
            == (106, 53, 44, 415), "Incomplete relocated observations")
    minimum = min(r["covered_enumerated_mass"] for r in rows)
    tail = max(r["omitted_tail_mass"] for r in rows)
    uncertainty = f"""### Conditional Count Model Does Not Establish Native Uncertainty

A restricted mathematical experiment examines a fixed perfect-recall truth
mass and Poisson shared/method-only whole-panel false-positive counts. Its
target is the difference of F1 ratios at expected count vectors, not expected
sample F1 or generalization from native families. Simultaneous component
Poisson mean limits are sharply projected through that ratio difference.
Zero observed errors retain a positive rate bound, rather than a zero-width
plug-in interval. The model is not fitted to native errors.

All 12 prespecified independent-Poisson enumeration cells passed; minimum
covered probability mass was {minimum:.12f}, with maximum omitted tail
{tail:.12e}. This finite panel is conservative, not proof of a native law
or floating-point certification. Both non-Poisson common-shock controls had
zero coverage despite having the same mean-count target form. Imperfect
recall, unequal family sizes and shared-clade dependence remain outside the
construction. Earlier failed dyadic screens are not replaced. No native VGNC,
TreeFam, GO/EC/FAS or secondary-mean confidence interval is admitted, and no
benchmark statistic, resampling unit or default is changed.
[Frozen model and derivation](PAIRED_POISSON_F1_PROTOCOL_20261010.md),
[all results and failure controls](paired_poisson_f1_validation_20261010_v1.json),
[execution and independent tail/corner readback](paired_poisson_f1_execution_20261010_v1.json).

"""
    replay = """### Relocated Fragment Raw Stage Replay

The selected fragment-stage postprocessing was executed from a separately
located component containing 1,783 pinned artifact files (39,428,185 payload
bytes) and unchanged native-reader kernels. Its complete result matches the
retained native readback: 106 observations, 53 selected cases, 44 contexts and
415 original input references. The first ownership-map failure is preserved;
the corrected map includes all configured ownership FASTAs. Python-open
auditing observed no original-artifact fallback, not operating-system
containment. This is raw-stage postprocessing portability, not a repeat of
historical inference, new biological evidence or whole-study restoration.

The component remains local. The Git repository contains its map, result and
reproduction instructions, not the raw payload; no public deposition or
redistribution clearance is inferred. The existing 30-page review PDF remains
the rendered copy of the prior manuscript, not this extended text.
[Replay result and preserved failure](CONTROLLED_FRAGMENT_RELOCATION_RESULT_20261010.md),
[complete copied readback](controlled_fragment_relocated_readback_20261010_v2.json),
[reproduction workflow](../CONTROLLED_FRAGMENT_REPRODUCTION.md).

"""
    return uncertainty, replay


def manuscript(parent, docs):
    require(all(parent.count(anchor) == 1 for anchor in ANCHORS), "Ambiguous manuscript anchors")
    inserts = sections(docs)
    revised = parent
    for anchor, section in zip(ANCHORS, inserts):
        require(section.splitlines()[0] not in parent, "Already integrated")
        revised = revised.replace(anchor, section + anchor, 1)
    restored = revised
    for section in inserts:
        restored = restored.replace(section, "", 1)
    require(restored == parent, "Unscoped manuscript alteration")
    return revised, inserts


def run(root, output, receipt):
    root = Path(root).resolve()
    directory = root / "benchmark_tools/results"
    output, receipt = Path(output).absolute(), Path(receipt).absolute()
    for path in (output, receipt):
        if path.exists() or path.is_symlink():
            raise FileExistsError(path)
    require(output != receipt and output.suffix == ".md" and receipt.suffix == ".json"
            and output.parent.resolve() == receipt.parent.resolve() == directory, "Preserve relative links and destinations")
    refs, docs = {}, {}
    for key, (name, digest) in PINS.items():
        path = directory / name
        refs[key] = record(path)
        require(refs[key]["sha256"] == digest, "Changed input " + key)
        text = path.read_text()
        docs[key] = text if key == "parent" else json.loads(text)
    revised, inserts = manuscript(docs["parent"], docs)
    for ref in refs.values():
        require(record(ref["path"]) == ref, "Input changed during integration")
    with output.open("x") as stream:
        stream.write(revised)
    result = dict(schema="publication_completed_supplements_v1", inputs=refs,
                  source=record(__file__), manuscript=record(output), inserted_sections=list(inserts),
                  parent_unchanged_except_insertions=True, new_inference_or_scoring=False,
                  new_bootstrap_draws=0, native_intervals_admitted=False, publication_ready=False,
                  manuscript_rendered=False, visual_reviewed=False)
    with receipt.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return result
