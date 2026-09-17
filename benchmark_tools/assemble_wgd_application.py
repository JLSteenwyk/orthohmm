"""Score checksum-pinned, admitted WGD outputs and export the complete cohort."""

import argparse
import csv
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.read_wgd_native_groups import read_orthohmm, read_orthofinder_root_ids, read_species_table
from benchmark_tools.report_wgd_application import DIAGNOSTIC, METHODS, assemble
from benchmark_tools.report_ygob_validation import read_checkpoint
from benchmark_tools.run_wgd_application import pinned
from benchmark_tools.snapshot_orthohmm_input_order import record
from benchmark_tools.validate_scaling_outputs import input_universe

ADMISSIONS = (
    ("biological_wgd_high_admission_20260917.json", "10b5d3feb7cd5d208856dc59f0fa71c614afb7d4e75025218ccca912fd2ad1b4"),
    ("biological_wgd_satellite_admission_20260917.json", "7c95f62911ffc906a54e5d9f2857883ed05920deec2d9508b459e15be70db258"),
    ("biological_wgd_orthofinder_admission_v2_20260917.json", "b2b7769e41310f17a9d21d3cedf593a47fa517ee39195ec4ff00f4d27722df5f"),
    ("biological_wgd_sonic_admission_20260917.json", "a1bf34304dd816024dabe009f14093f4bbc9aa733869bc5bc9652a2888e55432"),
)
LABELS = dict(zip((*METHODS, DIAGNOSTIC), ("OrthoHMM high sensitivity", "OrthoHMM phylogeny",
                                        "OrthoFinder full (root HOGs)", "SonicParanoid", "OrthoFinder MCL checkpoint (diagnostic)")))


def unchanged(item):
    if record(item["path"]) != item:
        raise ValueError("Admitted artifact changed: " + item["path"])


def artifact(admission, basename):
    records = []
    for row in admission["native"]["checked_files"]:
        if Path(row["path"]).name == basename and row not in records:
            records.append(row)
    if len(records) != 1:
        raise ValueError("Missing or ambiguous admitted artifact: " + basename)
    return Path(records[0]["path"])


def run(repo):
    results = repo / "benchmark_tools/results"
    spec = pinned({"path": str(results / "biological_wgd_execution_20260917.json"),
                   "sha256": "c43704020c56ead316678461f5be3e8d4efc43a56ff5ced4f0e3cfd3c189025c"})
    plan = pinned(spec["command_plan"])
    prepared = pinned(plan["inputs"])
    if len(prepared["cohort_pairs"]) != 240 or prepared["eligible_counts"] != {"all": 240, "split": 239, "homolog_support": 231}:
        raise ValueError("Changed fixed cohort population")
    for item in (prepared["reference"], prepared["cohort"], prepared["protocol"], plan["contrast_protocol"]):
        unchanged(item)
    owners, species = input_universe({**prepared, "proteomes": 4})
    reference = json.loads(Path(prepared["reference"]["path"]).read_text())
    admitted, provenance = {}, []
    for method, (name, digest) in zip(METHODS, ADMISSIONS):
        path = results / name
        admission = pinned({"path": str(path), "sha256": digest})
        if admission["method"] != method or admission["status"] not in {
                "orthohmm_application_native_outputs_admitted", "comparator_application_native_outputs_admitted"}:
            raise ValueError("Wrong method admission")
        for item in admission["native"]["checked_files"] + [admission["receipt"]]:
            unchanged(item)
        if method == "orthohmm_high_sensitivity":
            groups = read_orthohmm(artifact(admission, "orthohmm_orthogroups.txt"), "named_groups", owners)
        elif method == "orthohmm_satellite_v2":
            groups = read_orthohmm(artifact(admission, "orthohmm_root_hogs.tsv"), "root_hogs", owners)
        elif method == "orthofinder_full":
            ids = artifact(admission, "SequenceIDs.txt")
            groups = read_orthofinder_root_ids(artifact(admission, "N0.ids.tsv"), ids, owners, {s: s for s in species})
            clusters = [Path(r["path"]) for r in admission["native"]["checked_files"]
                        if Path(r["path"]).name.startswith("clusters_OrthoFinder_I") and r["path"].endswith("_id_pairs.txt")]
            if len(clusters) != 1:
                raise ValueError("Ambiguous admitted MCL checkpoint")
            admitted[DIAGNOSTIC] = {"status": "admitted", "groups": read_checkpoint(clusters[0], ids, owners)}
        else:
            columns = {Path(r["path"]).name: Path(r["path"]).stem for r in prepared["inputs"]}
            groups = read_species_table(artifact(admission, "ortholog_groups.tsv"), "sonicparanoid", owners, columns)
        admitted[method] = {"status": "admitted", "groups": groups}
        provenance.append(record(path))
    report = assemble(prepared, reference, owners, admitted)
    report.update(admissions=provenance, input_manifest=plan["inputs"], protocol=prepared["protocol"],
                  contrast_protocol=plan["contrast_protocol"],
                  sources=[record(Path(__file__).with_name(name)) for name in
                           ("assemble_wgd_application.py", "report_wgd_application.py", "score_wgd_application.py",
                            "bootstrap_wgd_application.py", "read_wgd_native_groups.py")])
    return report


def write_table(report, path):
    methods = (*METHODS, DIAGNOSTIC)
    lookup = {m: {tuple(r["orf_pair"]): r for r in report["methods"][m]["rows"]} for m in methods}
    common = ("orf_pair", "experimental_class", "split_eligible", "reference_eligible", "input_reasons", "reference_reasons")
    measures = ("assignment_state", "anchor_groups", "anchor_group_sizes", "separation_rate",
                "supported_separation_rate", "mean_non_scer_coverage", "homolog_support_by_anchor",
                "coverage_numerator", "coverage_denominator", "foreign_pillar_members", "unmapped_members",
                "union_size", "pillar_native_group_count", "unassigned_pillar_members")
    fields = list(common) + [f"{m}.{key}" for m in methods for key in measures]
    with path.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for pair in sorted(lookup[METHODS[0]]):
            row = {key: lookup[METHODS[0]][pair].get(key) for key in common}
            row.update({f"{m}.{key}": lookup[m][pair][key] for m in methods for key in measures})
            writer.writerow({key: json.dumps(value, sort_keys=True) if isinstance(value, (list, dict, bool)) else value
                             for key, value in row.items()})


def markdown(report):
    lines = ["# WGD Biological Application", "", "Development-exposed application; not orthology F1 or independent generalization.",
             "All 240 experimental pairs retained; 239 input-eligible, 231 shared-pillar eligible.", "",
             "| Method | Separated /239 | Supported /231 | Mean homolog coverage (%) | Unsupported splits /231 | Foreign-pillar cases /231 |",
             "| --- | ---: | ---: | ---: | ---: | ---: |"]
    for name in (*METHODS, DIAGNOSTIC):
        entry = report["methods"][name]
        if entry["status"] != "evaluated":
            lines.append(f"| {LABELS[name]} | unavailable | unavailable | unavailable | unavailable | unavailable |")
            continue
        summary = entry["summary"]
        rows = [r for r in entry["rows"] if r["reference_eligible"]]
        separated = summary["assignment_states"].get("separated", 0)
        supported = sum(r["supported_separation_rate"] for r in rows)
        coverage = summary["endpoints"]["mean_non_scer_coverage"]["mean"] * 100
        unsupported = sum(r["separation_rate"] == 1 and r["supported_separation_rate"] == 0 for r in rows)
        foreign = sum(bool(r["foreign_pillar_members"]) for r in rows)
        lines.append(f"| {LABELS[name]} | {separated} | {supported} | {coverage:.3f} | {unsupported} | {foreign} |")
    lines += ["", "Coverage can be high for merged groups. Foreign-pillar counts exclude unmapped members; neither is orthology precision.",
              "", "## Paired Differences", "", "All contrasts use 231 pairs. 20,000 paired pillar replicates; seed 20260920; differences in percentage points.",
              "", "| First minus second | Endpoint | Difference | Bonferroni12 interval |", "| --- | --- | ---: | --- |"]
    for row in report["uncertainty"]["comparisons"]:
        bounds = row.get("bonferroni12_pp")
        interval = f"[{bounds[0]:.3f}, {bounds[1]:.3f}]" if bounds is not None else "unavailable"
        delta = f"{row['difference_pp']:.3f}" if row["difference_pp"] is not None else "unavailable"
        lines.append(f"| {LABELS[row['first']]} minus {LABELS[row['second']]} | {row['endpoint']} | {delta} | {interval} |")
    lines += ["", "## Limitations", "", *["- " + text for text in report["limitations"]]]
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("repo", "output", "table", "markdown"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    if any(p.exists() for p in (args.output, args.table, args.markdown)):
        raise FileExistsError("Refusing to overwrite application results")
    report = run(args.repo.resolve())
    with args.output.open("x") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    write_table(report, args.table)
    with args.markdown.open("x") as handle:
        handle.write(markdown(report))
