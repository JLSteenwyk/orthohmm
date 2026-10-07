"""Add verified error-stratum evidence to a new, explicitly provisional manuscript."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess

PINS = {
    "parent": ("PUBLICATION_MAIN_TEXT_20261006_v2.md",
               "2804448cbd99ffa3164fcb9e997640f966103a827c9da761bfb7c3732d79a9c1"),
    "fixed": ("native_qfo_three_cell_strata_20261007_v1/report.json",
              "55ed72ce771c756424877cfa412da329ff104efb7dca4b35113b0d80703a59e5"),
    "fixed_reader": ("native_qfo_three_cell_strata_readback_20261007_v2.json",
                     "6f3eef63c0e96e1d85805fcfb70ed9b66fbb1db76a08ca2a105483422a0f6969"),
    "figure": ("native_qfo_three_cell_strata_figure_20261007_v1/manifest.json",
               "e0157657d24ed39da8f8e39a41590311ec480537ac5aa35e06483d758ef85b02"),
    "figure_review": ("native_qfo_three_cell_strata_figure_review_20261007_v1.json",
                      "9f40fc01d48c92c2d8e2d6bc42d67c395e9130a5929f43bb94d306797abe2556"),
    "distance": ("swiss_model_divergence_strata_20261007_v1/report.json",
                 "ddfa27757c4c256df9c881079ac352c37d7e17646f75fbd417b624e7980be78f"),
    "distance_reader": ("swiss_model_divergence_strata_readback_20261007_v1.json",
                        "186c0151137fd825513ef36c27df595b1a0f1cd4b2b4bf68de4b9665ad328d0e"),
    "feature_reader": ("swiss_model_divergence_readback_23932_v1.json",
                       "45cb12b3d201daf56e3fee5c5833153a69dd1732a2fee1c28d59261130eeb312"),
    "bibliography": ("publication_bibliography_20260920_v5.csl.json",
                     "bc2669a4c291a86473c88db22cc575b46eca55b6ec6b049c9c2cf860cb372aef"),
    "citation": ("iqtree3_citation_crossref_20261007.json",
                 "628a6ab8ac1272f6215112a15784c02a1f975762eb0f2d0465ac9fd1b5391df1"),
}
OLD_SELECTOR = "This source is [mechanically generated](../prepare_native_main_text_v2.py)\n"
NEW_SELECTOR = "This source is [mechanically extended](../prepare_error_strata_main_text.py)\n"
ANCHORS = dict(methods="### Shared-Host Resource Measurement\n",
               results="### Synthetic Null Tails Depend On Composition\n",
               discussion="Original TreeFam-A family mappings and complete source trees remain unavailable.\n",
               availability="## References\n")
HEADER = """Error-stratum evidence revision of 7 October 2026. Not submission-ready.
This version preserves the [reviewed 6 October v2 source](PUBLICATION_MAIN_TEXT_20261006_v2.md)
and its earlier review/rc5 archive unchanged. The complete three-cell fixed-bin
figure and model-based divergence analysis are later additions, not retroactive
members of rc5. This source-generation step does not establish a new HTML/PDF
review, complete executable study release, independent validation or publication
readiness. Dated execution and review receipts govern those separate claims.
The earlier opening below is retained as chronological history.

"""


def require(condition, message):
    if not condition:
        raise ValueError(message)


def record(path):
    path = Path(path).resolve(strict=True)
    return dict(path=str(path), bytes=path.stat().st_size,
                sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def validate(docs):
    fixed, distance = docs["fixed"], docs["distance"]
    require(fixed["memberships"] == distance["memberships"] and len(fixed["memberships"]) == 18
            and sum(map(len, fixed["memberships"].values())) == 563
            and fixed["family_rows"] == distance["family_rows"], "Changed native family/count universe")
    require(len(fixed["rows"]) == 60 and len(fixed["differences"]) == 40
            and docs["fixed_reader"]["family_rows_checked"] == 54
            and docs["fixed_reader"]["score_rows_checked"] == 60
            and docs["fixed_reader"]["differences_checked"] == 40, "Incomplete fixed-stratum admission")
    require(len(distance["rows"]) == 9 and len(distance["differences"]) == 6
            and docs["distance_reader"]["status"] == "projection_verified"
            and docs["distance_reader"]["score_rows_checked"] == 9
            and docs["distance_reader"]["differences_checked"] == 6, "Incomplete distance admission")
    require(docs["feature_reader"]["status"] == "features_verified"
            and docs["feature_reader"]["families_checked"] == 18
            and docs["feature_reader"]["proteins_checked"] == 563
            and docs["feature_reader"]["pairs_checked"] == 10765
            and docs["feature_reader"]["strata"] == distance["bins"]
            and docs["feature_reader"]["features"] == distance["features"], "Incomplete feature admission")
    require(distance["model"] == "WAG+G4" and distance["seed"] == 20261007
            and distance["unit"] == "model_estimated_expected_amino_acid_substitutions_per_site"
            and distance["cells"][1]["timing_eligible"] is False
            and distance["cells"][1]["timing_admitted"] is False, "Changed model/native timing scope")
    for key in ("fixed", "fixed_reader", "figure", "distance", "distance_reader", "feature_reader"):
        require(all(docs[key][k] is False for k in ("new_uncertainty", "independent_confirmation",
                    "new_accuracy_or_resource_admission", "publication_ready"))
                and docs[key]["new_bootstrap_draws"] == 0, "Inflated reporting scope")
    figure, review = docs["figure"], docs["figure_review"]
    require(len(figure["plotted_rows"]) == 66 and len(figure["excluded_bins"]) == 9
            and review["full_png_actually_viewed"] is True
            and review["full_pdf_page_raster_actually_viewed"] is True
            and review["publication_ready"] is False
            and review["old_figures_or_archives_rebuilt"] is False, "Incomplete existing figure review")
    candidate = [r for r in fixed["differences"] if r["contrast"] == "C_at_P0_R0"
                 and r["status"] == "descriptive"]
    require(len(candidate) == 15 and all(r["PPV"] < 0 and r["TPR"] > 0 for r in candidate),
            "Unsupported all-bin precision/recall statement")


def blocks(docs):
    validate(docs)
    distance = docs["distance"]
    cutoff = distance["median_family_distance"]
    methods = """### Fixed-Stratum SwissTrees Error Analyses

A complete three-cell diagnostic projects already-admitted family counts into
the previously frozen sequence, Pfam and reference-duplication descriptors.
All18families/563canonical members and all54count rows, including TN, are
retained. Original prior-adjusted macro-family precision/recall followed by
harmonic F1 is preserved; no raw scorer, bootstrap or original tree traversal
is rerun. All60score and40conditional-difference records are retained, including
empty bins as NA. Overlapping bins and repeated all-family records are not
independent tests. [Fixed-bin protocol](NATIVE_QFO_THREE_CELL_STRATA_PROTOCOL_20261007.md),
[complete independent readback](native_qfo_three_cell_strata_readback_20261007_v2.json).

A separately prespecified feature procedure reused every retained reference
alignment, without realignment, masking or filtering taxa/sites. Installed
IQ-TREE3.0.1 [@iqtree3_2026] inferred one tree per family using fixed WAG+G4,
seed20261007, one thread and identical-tip retention, without ModelFinder,
supplied/species trees, support resampling or a clock. The primary family
descriptor is the median of all unordered tip-pair patristic distances in
model-estimated expected amino-acid substitutions per site; paralog,
within-species and zero-distance pairs remain. A median-of-family-medians
cutoff assigns <=median to the lower bin and >median to the higher bin, with
no outcome-selected alternatives. Every selected family must succeed before
any cutoff/projection. The feature stage never reads prediction outcomes.
Independent edge-split summation checks all10765distances using shared
Bio.Phylo parsing, not a second inference/model method. A separate exact-rational
reader verifies every derived score and percentage-point contrast.
[Prospective model-distance protocol](SWISS_MODEL_DIVERGENCE_PROTOCOL_20261007.md),
[complete feature readback](swiss_model_divergence_readback_23932_v1.json),
[independent score readback](swiss_model_divergence_strata_readback_20261007_v1.json).

"""
    rows = []
    for row in distance["differences"]:
        label = {"all": "All families", "lower_or_equal_median": "Lower/equal distance",
                 "higher_than_median": "Higher distance"}[row["stratum"]]
        contrast = {"R_at_P0_C0": "R at P0/C0", "C_at_P0_R0": "C at P0/R0"}[row["contrast"]]
        values = ["NA" if row[m + "_pp"] is None else f"{row[m + '_pp']:+.3f}" for m in ("F1", "PPV", "TPR")]
        rows.append("| " + " | ".join([label, str(row["families"]), contrast, *values]) + " |")
    results = """### Fixed-Stratum Errors And Estimated Sequence Divergence

Candidate expansion lowers precision and increases recall in every nonempty
fixed sequence/domain/reference-duplication bin. These records overlap and
include repeated all-family summaries; they are not multiple independent
confirmations. The [complete fixed-bin figure](native_qfo_three_cell_strata_figure_20261007_v1/fixed_stratum_contrasts.pdf)
shows both conditional contrasts and all three metrics in all11nonempty,
nonredundant bins (66points). Nine excluded records are five empty bins and
four redundant all-family copies, not outcome-based omissions. Complete
tables retain every record. Dots are descriptive, not confidence intervals,
and the two panels use explicitly different percentage-point axis ranges.
[Full table/readback](NATIVE_QFO_THREE_CELL_STRATA_RESULT_20261007.md),
[actual figure review](NATIVE_QFO_THREE_CELL_STRATA_FIGURE_RESULT_20261007.md).

The fixed-model distance procedure succeeded for all18families, without a
timeout or retry. The median family-distance cutoff was approximately
""" + f"{cutoff:.11f}" + """ substitutions/site, giving nine families in each bin.
The complete conditional differences were:

| Distance Stratum | Families | Conditional Contrast | Delta F1 (pp) | Delta Precision (pp) | Delta Recall (pp) |
| --- | ---: | --- | ---: | ---: | ---: |
""" + "\n".join(rows) + """

Expansion trades precision for recall in both distance bins, with a small
higher-distance F1 increase but a lower-distance and overall decrease.
Reconciliation raises F1/precision and lowers recall in both. These exposed
subgroup point patterns have no new intervals or significance claim; neither
a causal mechanism nor a distance-dependent default follows. Initial HMM
search is on and downstream profile refinement is off in all three cells.
C and R contrasts are conditional, not an identified interaction.
[Complete model-distance result](SWISS_MODEL_DIVERGENCE_RESULT_20261007.md),
[complete scores and differences](swiss_model_divergence_strata_20261007_v1/TABLE.md).

"""
    discussion = """The new error-stratum analyses describe development-exposed families, not
independent confirmation. Fixed-model branch-path distances are an explicit
estimated sequence-divergence descriptor, not calibrated biological time,
known ancestral history, verified model adequacy or a causal explanation.
Paralogy, sampling, composition, domains/length and gaps/alignment error may
affect them. Independent arithmetic does not eliminate shared parsing,
alignment/model assumptions or heuristic topology uncertainty. Length and
composition remain proxies rather than literal fragment truth; Pfam summaries
are not complete validated architectures, and informative reference-node
duplication fractions are not complete ancestral histories. None of these
results promotes a new default or closes the remaining descriptor and
generalization gaps.

"""
    availability = """The new three-cell fixed-bin figure and model-distance supplements retain
all positive, negative, neutral and empty outcomes. A small
[model-distance evidence archive](swiss_model_divergence_evidence_23932_v1.tar.gz)
contains all18selected alignments and the complete new tree/log/pair/receipt
payloads; [execution receipt](swiss_model_divergence_execution_23932_v1.json)
and independent readers state the exact scope. Auxiliary inference Slurm
TotalCPU and MaxRSS were unavailable, not measured zeros or native tool timing.
The IQ-TREE3citation uses retained publisher-deposited
[Crossref metadata](iqtree3_citation_crossref_20261007.json) for the current
Wong et al.2026 paper, rather than the installed binary's old submitted-paper
notice. The executed binary remains3.0.1; no software update or rerun follows.
The new bibliography preserves every prior entry and adds this reference only.
Earlier manuscript reviews and rc5 remain unchanged; their dated receipts do
not establish review or package inclusion of this later source revision.
This is reporting integration, not transitive raw/dependency closure, a hermetic
whole-study release, redistribution clearance, independent confirmation or DOI
deposition. Future render/package receipts must cover their actual new inputs.

"""
    return dict(methods=methods, results=results, discussion=discussion, availability=availability)


def extend(parent, docs):
    additions = blocks(docs)
    require(parent.startswith("# OrthoHMM:") and parent.count("\n\n") > 1, "Invalid parent draft")
    heading, rest = parent.split("\n\n", 1)
    output = heading + "\n\n" + HEADER + rest
    require(output.count(OLD_SELECTOR) == 1, "Changed parent generator selector")
    output = output.replace(OLD_SELECTOR, NEW_SELECTOR)
    for key, anchor in ANCHORS.items():
        require(output.count(anchor) == 1, "Changed insertion anchor: " + key)
        require(additions[key] not in output, "Already-integrated block")
        output = output.replace(anchor, additions[key] + anchor)
    return output, additions


def citation_entry(document):
    data = document["message"]
    require(document["status"] == "ok" and data["DOI"] == "10.1093/molbev/msag117"
            and data["type"] == "journal-article" and data["issued"]["date-parts"][0][0] == 2026,
            "Unexpected publisher-deposited citation")
    return dict(id="iqtree3_2026", type="article-journal", title=data["title"][0],
                author=[{k: a[k] for k in ("given", "family")} for a in data["author"]],
                issued=data["issued"], DOI=data["DOI"], URL="https://doi.org/" + data["DOI"],
                volume=data["volume"], issue=data["issue"], page=data["article-number"],
                **{"container-title": data["container-title"][0]})


def generate(repo, output, bibliography, receipt, source_commit):
    repo = Path(repo).resolve()
    for path in (output, bibliography, receipt):
        require(not Path(path).exists() and not Path(path).is_symlink(), "Existing generation output")
    require(len({Path(p).resolve() for p in (output, bibliography, receipt)}) == 3, "Distinct outputs required")
    docs, refs = {}, {}
    for key, (name, sha) in PINS.items():
        ref = record(repo / "benchmark_tools/results" / name)
        require(ref["sha256"] == sha, "Changed pinned input: " + key)
        refs[key] = ref
        if name.endswith(".json"):
            docs[key] = json.loads(Path(ref["path"]).read_text())
    require(docs["fixed_reader"]["report"] == refs["fixed"]
            and docs["distance_reader"]["report"] == refs["distance"]
            and docs["figure"]["inputs"]["report"] == refs["fixed"]
            and docs["figure_review"]["manifest_sha256"] == refs["figure"]["sha256"], "Broken direct admission link")
    figure_outputs = docs["figure"]["outputs"]
    for ref in figure_outputs:
        require(record(ref["path"]) == ref, "Changed reviewed figure output")
    require(Path(__file__).read_bytes() == subprocess.check_output(["git", "show", source_commit
            + ":benchmark_tools/prepare_error_strata_main_text.py"], cwd=repo), "Uncommitted integration source")
    text, additions = extend(Path(refs["parent"]["path"]).read_text(), docs)
    entries = docs["bibliography"]
    addition = citation_entry(docs["citation"])
    require(addition["id"] not in [e["id"] for e in entries], "Duplicate new bibliography identity")
    with Path(output).open("x", encoding="utf-8") as stream:
        stream.write(text)
    with Path(bibliography).open("x", encoding="ascii") as stream:
        json.dump([*entries, addition], stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    result = dict(schema="error_strata_main_text_generation_v1", status="provisional_source_generated",
                  source=record(__file__), source_commit=source_commit, inputs=refs,
                  figure_outputs=figure_outputs, output=record(output), bibliography=record(bibliography),
                  blocks_sha256={k: hashlib.sha256(v.encode()).hexdigest() for k, v in additions.items()},
                  old_evidence_modified=False, new_scoring_or_admission=False, native_inference_reexecuted=False,
                  new_bootstrap_draws=0, publication_ready=False, new_render_or_package_proved=False,
                  bibliography_entries_preserved=len(entries), bibliography_entries_added=1,
                  limitations=["Source integration only, not browser/PDF review or package restoration.",
                               "No new scientific validation, inference, interval, tuning or default change."])
    with Path(receipt).open("x", encoding="ascii") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("repo", "output", "bibliography", "receipt"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--source-commit", required=True)
    args = parser.parse_args()
    result = generate(args.repo, args.output, args.bibliography, args.receipt, args.source_commit)
    print(json.dumps(dict(status=result["status"], output=result["output"], bibliography=result["bibliography"])))
