"""Integrate verified representatives without replacing frozen scientific results."""

from collections import Counter
import json
from pathlib import Path

from benchmark_tools.export_native_qfo_three_cell_strata import record, require, write_tsv


PINS = {
    "selection": ("benchmark_tools/results/controlled_fragment_trace_selection_20261009_v1/selection.json", "e072dbff2c570f1926bf46ba2ad1dcdbf31ce5197d8c3461dc3576edb41b8594"),
    "stages": ("benchmark_tools/results/controlled_fragment_stage_trace_20261010_v1/report.json", "0d5f5777effd9e33a5fe14a188e07af118f0c783b91bc0655264dfd707751b8c"),
    "readback": ("benchmark_tools/results/controlled_fragment_stage_readback_20261010_v2.json", "2d37e635fb9d9f08d1b1857d829e0798c19452936785b92e5dd40a058658da35"),
    "execution": ("benchmark_tools/results/controlled_fragment_trace_execution_20261010_v1.json", "221ac09bcb21ee9ef1192f9fff7d948d1c1ea0a7f846cfc2e9dbb7a2aa7dad32"),
    "parent": ("benchmark_tools/results/controlled_fragment_manuscript_20261009_v1.md", "71b499efab6a54b2d8557ed7c04241001abf08b436e4d56344e6fdc14f93e41a"),
}
METHODS = ("orthohmm_high_sensitivity", "orthohmm_satellite_v2", "orthofinder_full", "orthofinder_sequence_only")
LABELS = dict(zip(METHODS, ("HS", "PHY", "OF", "OF checkpoint")))
ANCHORS = ("### Fresh Native Ablation Protocol\n", "### Gene-Tree And Event-History Controls Localize Errors\n", "## Reproducibility And Availability\n")
OLD_STATUS = "Representative search/candidate/group/reconciliation tracing\nremains separate work."
NEW_STATUS = "Representative search/candidate/group/reconciliation tracing\nis reported in the next section, with unobserved stages retained explicitly."
FIELDS = ("method", "observed_location", "baseline_representatives", "fragment_representatives")


def tables(docs):
    selected, stages, reader, execution = (docs[k] for k in ("selection", "stages", "readback", "execution"))
    require(selected["planned_bins"] == 84 and selected["selected_records"] == selected["unique_method_pair_records"] == 53
            and Counter(b["status"] for b in selected["bins"]) == {"selected":53, "empty_bin":31}, "Changed representative selection")
    require(stages["selection"]["sha256"] == PINS["selection"][1] and reader["report"]["sha256"] == PINS["stages"][1]
            and reader["status"] == "independent_native_stages_and_na_table_verified"
            and reader["stage_rows"] == stages["stage_rows"] == 106 and reader["selected_cases"] == 53
            and reader["summary"] == stages["summary"], "Unverified stage evidence")
    terminal = execution["independent_stage_readback_recovery"]["terminal"]
    observed = json.loads(terminal["output"])
    require(terminal["exit_code"] == 0 and observed["status"] == reader["status"] and observed["stage_rows"] == 106
            and execution["independent_stage_readback_recovery"]["readback_sha256"] == PINS["readback"][1]
            and execution["independent_stage_readback_failed"]["terminal_output_transcription"]["exit_code"] == 1,
            "Mixed actual readback history")
    require(all(d[k] is False for d in (selected, stages, reader, execution) for k in ("publication_ready", "new_inference_or_scoring")), "Changed scientific scope")
    counts = Counter()
    require(len(stages["cases"]) == len(selected["cases"]) == 53, "Incomplete cases")
    for case, original in zip(stages["cases"], selected["cases"]):
        require({k:v for k,v in case.items() if k != "stages"} == original, "Changed case identity")
        for arm in ("baseline", "fragment"):
            counts[case["method"], case["stages"][arm]["observed_location"], arm] += 1
    rows = [dict(method=m, observed_location=l, baseline_representatives=counts[m,l,"baseline"], fragment_representatives=counts[m,l,"fragment"])
            for m in METHODS for l in sorted({l for method,l,arm in counts if method == m})]
    require(sum(r["baseline_representatives"] for r in rows) == sum(r["fragment_representatives"] for r in rows) == 53, "Incomplete location table")
    cases = {c["case_id"]:c for c in stages["cases"]}
    for case_id, method, endpoints, fields in (
        ("Case0000", METHODS[0], 0, dict(hit_forward=False, hit_reverse=False, graph_connected=True, same_seed_group=False, observed_location="group_separation")),
        ("Case0026", METHODS[1], 2, dict(hit_forward=True, hit_reverse=True, graph_connected=True, pair_event="duplication", observed_location="observed_duplication_exclusion")),
        ("Case0030", METHODS[2], 0, dict(same_seed_group=True, pair_event=None, search_status="unavailable_matched_adapter", observed_location="within_mcl_native_pair_exclusion")),
    ):
        case = cases[case_id]
        require(case["method"] == method and case["fragment_endpoints"] == endpoints and case["category"] == "fragment_fn" and case["truth"] is True
                and case["stages"]["fragment"]["predicted"] is False and all(case["stages"]["fragment"][k] == v for k,v in fields.items()), "Changed illustrative observation")
    return rows


def sections(docs):
    rows = tables(docs)
    methods = """### Retained Fragment Representative Protocol

Representative selection is retrospective: aggregate results had already been
inspected. Four methods, three truncated-endpoint strata and seven categories
define 84 bins. Categories are fragment FN/FP, new FN/FP relative to baseline,
recovered FN, removed FP and retained TP controls. Within each nonempty bin,
the minimum SHA256 of UTF-8 `method:category:seed:gene_a:gene_b`, with lexical
tie-breaking, selects one eligible pair-record across the ten retained seeds.
Thirty-one bins are empty, without replacement. Fifty-three selected records
are unique in this realization; the protocol permits category overlap and
does not make these observations independent or representative of prevalence.

Read-only tracing retains both arms for every selected case. Significant-hit
NumPy checkpoints provide exact directed rows/scores; the weighted graph and
groups provide direct edges, connectivity and seed membership. PHY additionally
uses saved candidate families, observed reconciliation-node ancestors, logged
merge constraints and final RootHOGs. Positive-paralogy excludes observed
duplications but permits mapping-only uncertain calls; membership filtering
is active only when saved constraints are nonempty. Final RootHOG separation
alone is not proof of native pair exclusion. Unambiguous families without tree
inference are separately marked as bypasses. OF is traced only through its
MCL checkpoint and native pair output, not a matched search or event adapter.

Independent native-file readers reproduce all 106 observations in 44 contexts,
including graph connectivity by BFS and bottom-up root/constraint reconstruction.
The initial readback refused the TSV's null serialization after native checks;
a separate tested correction accepts the original `NA` cells, not empty or
false fields, while preserving the failed reader and original report/table.
No inference, scoring, bootstrap or tuning is repeated.

"""
    lines = ["### Retained Fragment Pipeline Observations", "",
             "The table counts the fixed representatives (HS 12, PHY 18, OF 11,",
             "OF checkpoint 12) in each arm, including errors, recoveries and retained",
             "TP controls. These are not rates, prevalence estimates or new accuracy",
             "endpoints. All 53 cases and empty bins remain available; no preferred",
             "stage explanation was used in selection.", "",
             "| Method | Recorded Decision | Baseline Representatives | Fragment Representatives |",
             "| --- | --- | ---: | ---: |"]
    lines.extend(f"| {LABELS[r['method']]} | {r['observed_location'].replace('_', ' ')} | {r['baseline_representatives']} | {r['fragment_representatives']} |" for r in rows)
    lines += ["", "Three illustrative records show distinct observable boundaries, not causal",
              "or biological validation. HS Case0000 (unflagged FN) lacks direct",
              "significant hits in either direction, yet its graph is connected and its",
              "groups are separate. PHY Case0026 (two truncated endpoints, FN) has both",
              "directed significant hits and connected graph support, but its observed",
              "ancestor is called duplication and excludes the native pair. OF Case0030",
              "(unflagged FN) is co-clustered at MCL but absent from native pairs; a",
              "matching event adapter is unavailable, so its tree-level cause is unknown.", "",
              "Search-hit absence does not distinguish prefilter rejection from scoring",
              "or significance thresholds. Graph/group support is not orthology. A",
              "recorded duplication establishes the pipeline's decision, not correctness",
              "of an inferred tree or the underlying evolutionary history.", "",
              "[Selection and empty bins](controlled_fragment_trace_selection_20261009_v1/selection.json),",
              "[all retained observations](controlled_fragment_stage_trace_20261010_v1/report.json),",
              "[stage table](controlled_fragment_stage_trace_20261010_v1/stages.tsv),",
              "[independent readback](controlled_fragment_stage_readback_20261010_v2.json),",
              "[actual execution and preserved refusal](controlled_fragment_trace_execution_20261010_v1.json).", "", ""]
    limitations = """### Representative Trace Limitations

Minimum-hash retrospective examples are not a random family sample or new
independent generalization test. Category overlap is allowed; realized unique
cases do not remove dependence. Selected counts cannot quantify the prevalence
or causal contribution of search, clustering or reconciliation errors. Native
node calls do not validate inferred trees against biology. Unavailable OF
search/event evidence remains null, not a negative result. This diagnostic
does not establish natural-fragment recovery, initial-HMM benefit, a default
parameter improvement, isolated timing performance or publication readiness.
The full benchmark uncertainty, failed-cell and shared-host limitations remain.

"""
    return (methods, "\n".join(lines), limitations), rows


def manuscript(parent, docs):
    require(all(parent.count(a) == 1 for a in ANCHORS) and parent.count(OLD_STATUS) == 1, "Ambiguous parent anchors/status")
    inserts, rows = sections(docs)
    revised = parent.replace(OLD_STATUS, NEW_STATUS, 1)
    for anchor, section in zip(ANCHORS, inserts):
        require(section not in parent, "Already integrated")
        revised = revised.replace(anchor, section + anchor, 1)
    restored = revised.replace(NEW_STATUS, OLD_STATUS, 1)
    for section in inserts:
        restored = restored.replace(section, "", 1)
    require(restored == parent, "Unscoped parent alteration")
    return revised, inserts, rows


def run(root, output, revised_path):
    root = Path(root).resolve()
    output, revised_path = Path(output).absolute(), Path(revised_path).absolute()
    for path in (output, revised_path):
        if path.exists() or path.is_symlink():
            raise FileExistsError(path)
    require(output.parent.resolve() == revised_path.parent.resolve() == root / "benchmark_tools/results", "Preserve relative links")
    docs, refs = {}, {}
    for key, (name, sha) in PINS.items():
        path = root / name
        refs[key] = record(path)
        require(refs[key]["sha256"] == sha, "Changed input " + key)
        docs[key] = path.read_text() if key == "parent" else json.loads(path.read_text())
    revised, inserts, rows = manuscript(docs["parent"], docs)
    source = record(__file__)
    output.mkdir()
    for name, section in zip(("methods", "results", "limitations"), inserts):
        (output / (name + "_section.md")).write_text(section)
    write_tsv(output / "representative_locations.tsv", rows, FIELDS)
    (output / "claims.md").write_text("# Retained Fragment Trace Claims\n\nSupported: 53 fixed representatives/106 observations independently read back; native stage decisions localized where observed.\n\nUnestablished: error prevalence, causal search or phylogeny defects, natural biological correctness, general superiority, initial-HMM benefit, independent confirmation or publication readiness.\n")
    revised_path.write_text(revised)
    for ref in [*refs.values(), source]:
        require(record(ref["path"]) == ref, "Input changed while writing")
    result = dict(schema="controlled_fragment_trace_integration_v1", inputs=refs, source=source,
                  outputs=[record(p) for p in sorted(output.iterdir())], manuscript=record(revised_path),
                  location_rows=len(rows), selected_cases=53, stage_rows=106,
                  bounded_status_replacement=dict(before=OLD_STATUS, after=NEW_STATUS),
                  parent_body_unchanged_except_insertions_and_recorded_status=True,
                  publication_ready=False, new_inference_or_scoring=False, new_bootstrap_draws=0,
                  manuscript_rendered=False, visual_reviewed=False)
    with (output / "manifest.json").open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return result
