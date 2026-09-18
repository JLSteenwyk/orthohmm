"""Re-audit corrected sequence-control counts and run the frozen SwissTrees protocol."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_qfo_sequence_swiss import audit
from benchmark_tools.bootstrap_qfo_sequence import bootstrap, METRICS
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen

PROTOCOL_SHA = "e853c2f9fbfe469e34fc6ebd0d00692c9fd9c91327750931e01131c0adb5e277"
SOURCES = {
    "audit_qfo_sequence_swiss.py": "e047fb1df9b7c7a860aa20ec4bc486191636b81a3f76517b48fa1ffb34e00631",
    "bootstrap_qfo_sequence.py": "d329964729862ff0d203b9d91542a5472121543d82f66d091d75ca0dd4e34202",
    "bootstrap_qfo_swiss_stages.py": "0598773a8bd04846f24e4cb27ab7b3529410c11451cbb399b9148a320fd06435",
    "audit_qfo_swiss_counts.py": "ff8e5314c98bdffce54732df94611e9d340d854e70392358c495b636b67452ae",
    "audit_qfo_factorial_swiss.py": "50a2e722c0bd57fbf8fb96412ce3e4f788132491983563dcc75b26032ee13f28",
    "audit_qfo_corrected_factorial_swiss.py": "375f2b4fd511c04e562827aebacb3b7522ca979f76b665652dd209940542f912",
    "export_qfo_corrected_factorial.py": "8a508b229d6aee560a65c7700d4281f4b8fd550b748ef9062e18bc9f84717ecd",
    "run_qfo_sequence_assessment.py": "8375b675c4baf6ab5235632b0dddae62975bf91879e6637b6adbecc0ca5f3ad6",
    "run_qfo_corrected_factorial_assessment.py": "ada6c4b62aef45940c9d8b3d7d670cd5ef224c356d48f8088f9b4892d4c7f164",
    "validate_qfo_native_assessment.py": "4786081a21ce04b8ae25f1f3e83ed8a1766472ff7e4476ae70752fe7f0656da5",
    "prepare_ob_candidate_neighborhood.py": "03843ed17fa9c44ea1ce5249782cdf009d8c7db9ceb85c62aa85bc2b4156ec29",
    "run_simulation_methods.py": "cccf579c221c51e8d03d512264efc3d0a61f4c0a598863d2423bdbd90f7921be",
}


def markdown(result):
    lines = ["# Corrected QfO Sequence-Control SwissTrees Intervals", "",
        f"{result['replicates']} shared family draws; seed {result['seed']}; six adjusted endpoints.",
        "Differences are candidate minus initial HMM, in raw 0-to-1 units.", "",
        "| Contrast | Metric | Difference | Nominal 95% CI | Adjusted CI | Wins/ties/losses |",
        "|---|---|---:|---|---|---|"]
    for contrast in result["comparisons"]:
        for metric in METRICS:
            row = contrast["metrics"][metric]
            intervals = ["[%.6f, %.6f]" % tuple(row[key])
                         for key in ("paired_percentile_ci", "bonferroni_percentile_ci")]
            lines.append(f"| {contrast['candidate']} - {contrast['reference']} | {metric} | {row['difference']:+.6f} | "
                + " | ".join(intervals) + f" | {row['family_wins']}/{row['family_ties']}/{row['family_losses']} |")
    return "\n".join([*lines, "", *["- " + s for s in result["limitations"]], ""])


def run(counts_path, counts_sha, protocol_path, output, markdown_path):
    if output.resolve() == markdown_path.resolve() or any(p.exists() or p.is_symlink() for p in (output, markdown_path)):
        raise FileExistsError("Require distinct fresh output paths")
    protocol = record(protocol_path)
    if protocol["sha256"] != PROTOCOL_SHA:
        raise ValueError("Changed frozen protocol")
    helpers = [record(Path(__file__).with_name(name)) for name in SOURCES]
    for name, item in zip(SOURCES, helpers):
        if item["sha256"] != SOURCES[name]:
            raise ValueError("Changed frozen implementation: " + name)
    data = read_frozen(counts_path, counts_sha)
    source = helpers[0]
    if data["source"] != source:
        raise ValueError("Counts do not identify the frozen auditor")
    inventory, baseline = data["admission_inventory"], data["baseline_audit"]
    checked = [record(counts_path), protocol, source, *helpers, inventory, baseline,
               *data["checked_inputs"], *data["helpers"]]
    for item in checked:
        check(item)
    # Reconstruct from raw evidence, not just the self-reported count-audit tag.
    rebuilt = audit(Path(inventory["path"]), inventory["sha256"], Path(baseline["path"]))
    if data != rebuilt:
        raise ValueError("Saved counts differ from fresh source-bound reconstruction")
    result = bootstrap(rebuilt)
    for item in checked:
        check(item)
    result.update(status="paired_corrected_sequence_swiss_intervals", uncertainty_admitted=True,
        counts=record(counts_path), protocol=protocol, source=record(__file__), helpers=helpers,
        checked_inputs=checked, admission_inventory=inventory, baseline_audit=baseline)
    result["limitations"] = [s for s in result["limitations"] if s != "Numerical output alone does not verify source provenance."]
    result["limitations"].append("Source-bound counts were reconstructed; upstream inference/scoring admissions were not rerun.")
    rendered = markdown(result)
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    with markdown_path.open("x") as stream:
        stream.write(rendered)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("counts", "protocol", "output", "markdown"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--counts-sha256", required=True)
    args = parser.parse_args()
    run(args.counts.resolve(), args.counts_sha256, args.protocol.resolve(), args.output.absolute(), args.markdown.absolute())
