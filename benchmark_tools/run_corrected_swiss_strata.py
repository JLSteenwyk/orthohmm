"""Reconstruct corrected primary-comparison counts before stratified resampling."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_qfo_corrected_factorial_swiss import audit as audit_factorial
from benchmark_tools.audit_qfo_corrected_swiss import audit as audit_comparator
from benchmark_tools.audit_qfo_swiss_counts import REFERENCE_SHA
from benchmark_tools.bootstrap_corrected_swiss_strata import bootstrap, METHODS
from benchmark_tools.bootstrap_qfo_factorial import CELLS
from benchmark_tools.prepare_corrected_swiss_sequence_strata import define_strata, PROTOCOL_SHA
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen

STRATA_SHA = "912363276fa5456f11aeff9f398fce71aec8d377c8c8eda7b144105945e554f1"
SOURCES = {
    "bootstrap_qfo_factorial.py": "be09876ae4de31b923818385bd3d40d8d5e7df177215fc45a2101c7d30d1c1fa",
    "prepare_qfo_corrected_orthofinder_pairs.py": "c804a4070a8835c98f94b5079184de9123d7c2d179700ecfcc0bb6d4da602cd5",
    "prepare_qfo_corrected_fastoma_pairs.py": "c4437cd13a61573100d5a3fabdbb2fb3d15c6e2f2570f13b4b4fc60631953069",
    "prepare_qfo_corrected_orthomcl_pairs.py": "87b229469466cf330259db7f1ec5cd03609001bbb5d594ef49583bd4f832b23e",
    "run_qfo_corrected_factorial_assessment.py": "ada6c4b62aef45940c9d8b3d7d670cd5ef224c356d48f8088f9b4892d4c7f164",
    "bootstrap_corrected_swiss_strata.py": "0ba971cf15ba61eb4c03d10c252fae25a3f5d3cf4f621eff3939ac81bb71db7e",
    "prepare_corrected_swiss_sequence_strata.py": "41c5cf70ca40debeb19476b636ea60ce115579e8015c29ada59738368383e26a",
    "audit_qfo_corrected_swiss.py": "f1ab13b089ca1cba095ce5e296456cf121067fbcda7b94aac7f6417da8654cf5",
    "audit_qfo_corrected_factorial_swiss.py": "375f2b4fd511c04e562827aebacb3b7522ca979f76b665652dd209940542f912",
    "audit_qfo_factorial_swiss.py": "50a2e722c0bd57fbf8fb96412ce3e4f788132491983563dcc75b26032ee13f28",
    "export_qfo_corrected_comparison.py": "4239f8d263295b313c5cc27b866bf7c6494feac44a8fbfe18e160bd92ff69bfb",
    "export_qfo_corrected_factorial.py": "8a508b229d6aee560a65c7700d4281f4b8fd550b748ef9062e18bc9f84717ecd",
    "audit_qfo_swiss_counts.py": "ff8e5314c98bdffce54732df94611e9d340d854e70392358c495b636b67452ae",
    "bootstrap_qfo_swiss_stages.py": "0598773a8bd04846f24e4cb27ab7b3529410c11451cbb399b9148a320fd06435",
    "validate_qfo_native_assessment.py": "4786081a21ce04b8ae25f1f3e83ed8a1766472ff7e4476ae70752fe7f0656da5",
    "prepare_ob_candidate_neighborhood.py": "03843ed17fa9c44ea1ce5249782cdf009d8c7db9ceb85c62aa85bc2b4156ec29",
    "run_simulation_methods.py": "cccf579c221c51e8d03d512264efc3d0a61f4c0a598863d2423bdbd90f7921be",
}


def assemble(factorial, comparator, strata):
    if (factorial["status"] != "corrected_qfo_factorial_swiss_counts_verified"
            or comparator["status"] != "corrected_comparator_swiss_counts_verified"
            or comparator["method"] != "orthofinder_full"
            or comparator["method_key"] != "orthofinder_3_1_5_full"):
        raise ValueError("Require corrected factorial and full OrthoFinder counts")
    for report in (factorial, comparator):
        if (report["reference"]["sha256"] != REFERENCE_SHA or report["shared_represented_genes"]
                or report["reference_relation_count"] != 10765):
            raise ValueError("Changed or overlapping reference universe")
    if (strata["status"] != "corrected_swiss_sequence_strata_prepared_unscored"
            or strata["prediction_statistics_evaluated"] is not False):
        raise ValueError("Require frozen input-only strata")
    families = strata["family_memberships"]
    rebuilt = define_strata(families, strata["genes"])
    if any(strata[k] != v for k, v in rebuilt.items()):
        raise ValueError("Strata disagree with frozen descriptors")
    if len(families) != 18 or sum(map(len, families.values())) != 563:
        raise ValueError("Changed sequence-stratum reference inventory")
    cells = factorial["cells"]
    if [r["cell"] for r in cells] != list(CELLS):
        raise ValueError("Changed corrected factorial inventory/order")
    selected = (cells[4], cells[7], comparator)
    counts = {}
    for method, report in zip(METHODS, selected):
        rows = report["families"]
        if len(rows) != 18 or {r["family"] for r in rows} != set(families):
            raise ValueError("Changed count family inventory")
        counts[method] = {}
        for row in rows:
            if sorted(row["represented_genes"]) != sorted(families[row["family"]]):
                raise ValueError("Counts do not match sequence-stratum members")
            counts[method][row["family"]] = row
    membership = {f: b for b, key in (("lower", "lower_entropy"), ("higher", "higher_entropy"),
                                     ("missing", "missing_entropy")) for f in strata["primary_strata"][key]}
    return counts, {f: membership[f] for f in sorted(families)}


def run(inventory_path, inventory_sha, of_path, of_sha, baseline_path, strata_path, protocol_path, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError("Require a fresh result path")
    helpers = [record(Path(__file__).with_name(name)) for name in SOURCES]
    if any(item["sha256"] != SOURCES[name] for name, item in zip(SOURCES, helpers)):
        raise ValueError("Changed frozen analysis implementation")
    protocol = record(protocol_path)
    if protocol["sha256"] != PROTOCOL_SHA:
        raise ValueError("Changed frozen protocol")
    strata = read_frozen(strata_path, STRATA_SHA)
    if strata["protocol"]["sha256"] != PROTOCOL_SHA:
        raise ValueError("Wrong strata protocol binding")
    checked = [*helpers, protocol, record(strata_path), record(inventory_path), record(of_path),
               record(baseline_path), strata["source"], strata["helper"], *strata["inputs"],
               strata["descriptor_update_protocol"], strata["native_sequence_audit"]]
    for item in checked:
        check(item)
    factorial = audit_factorial(inventory_path, inventory_sha, baseline_path)
    comparator = audit_comparator(of_path, of_sha, baseline_path)
    counts, membership = assemble(factorial, comparator, strata)
    result = bootstrap(counts, membership)
    checked.extend(factorial["checked_inputs"])
    checked.extend(comparator["checked_records"])
    for item in checked:
        check(item)
    result.update(status="corrected_swiss_primary_stratified_intervals", scientific_inputs_admitted=True,
                  uncertainty_admitted=True, source=record(__file__), helpers=helpers, checked_inputs=checked,
                  strata=record(strata_path), protocol=protocol,
                  reconstructed_counts=dict(factorial=factorial, orthofinder=comparator),
                  method_bindings=dict(high_sensitivity="corrected replay p1c0r0",
                                       phylogenetic="corrected p1c1r1 native pairs",
                                       orthofinder_full="corrected OrthoFinder 3.1.5 native pairs"))
    result["limitations"] = [s for s in result["limitations"] if not s.startswith("No prediction or input provenance")]
    result["limitations"].extend([
        "Raw scoring counts were reconstructed; upstream inference/scoring admissions were not rerun.",
        "Frozen sequence descriptors were reused and their bins recomputed; source FASTAs were not reread.",
        "Only the three primary configurations are included; all-method and secondary-stratum displays remain separate."])
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("inventory", "orthofinder-admission", "baseline", "strata", "protocol", "output"):
        parser.add_argument("--" + flag, type=Path, required=True)
    parser.add_argument("--inventory-sha256", required=True)
    parser.add_argument("--orthofinder-admission-sha256", required=True)
    args = parser.parse_args()
    run(args.inventory.resolve(), args.inventory_sha256, args.orthofinder_admission.resolve(),
        args.orthofinder_admission_sha256, args.baseline.resolve(), args.strata.resolve(),
        args.protocol.resolve(), args.output.absolute())
