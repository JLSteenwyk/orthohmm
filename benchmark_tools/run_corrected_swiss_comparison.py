"""Audit corrected raw comparator counts before the frozen paired analysis."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_corrected_swiss_comparison import audit
from benchmark_tools.bootstrap_corrected_swiss_comparators import bootstrap
from benchmark_tools.bootstrap_qfo_swiss_comparators import PROTOCOL_SHA
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

RELEASE_PROTOCOL_SHA = "b3603bc9b51ce1708f49a0ee816eb02d3cb52c66fb00f07eeb0d66f840a03e1e"
SOURCES = {
    "audit_corrected_swiss_comparison.py": "2a3e8681febff20aae43b101c55747d616b1bd45d2107fa2c1d04b06aeb67616",
    "bootstrap_corrected_swiss_comparators.py": "a5c339960936d0a41a66023059e4390db46e0cdc505194d209f56af060b97050",
    "audit_qfo_corrected_factorial_swiss.py": "375f2b4fd511c04e562827aebacb3b7522ca979f76b665652dd209940542f912",
    "audit_qfo_corrected_swiss.py": "f1ab13b089ca1cba095ce5e296456cf121067fbcda7b94aac7f6417da8654cf5",
    "audit_qfo_factorial_swiss.py": "50a2e722c0bd57fbf8fb96412ce3e4f788132491983563dcc75b26032ee13f28",
    "audit_qfo_swiss_counts.py": "ff8e5314c98bdffce54732df94611e9d340d854e70392358c495b636b67452ae",
    "audit_qfo_swiss_comparators.py": "be5b5b17fcd5591145165da2cac24d7aa424807564369eca691ba33916e15efc",
    "bootstrap_qfo_swiss_comparators.py": "c7b418c92a8b12dc06bb1a51880f7ef82c3451a8af4a876cd3c2cbd1802f9122",
    "bootstrap_qfo_swiss_stages.py": "0598773a8bd04846f24e4cb27ab7b3529410c11451cbb399b9148a320fd06435",
    "export_qfo_complete_comparison.py": "b2b0740a9b99eca52cffc59f211b4f94077d92b112d03b304a2b7f8c6849dfee",
    "export_qfo_corrected_comparison.py": "4239f8d263295b313c5cc27b866bf7c6494feac44a8fbfe18e160bd92ff69bfb",
    "export_qfo_corrected_factorial.py": "8a508b229d6aee560a65c7700d4281f4b8fd550b748ef9062e18bc9f84717ecd",
    "run_qfo_corrected_factorial_assessment.py": "ada6c4b62aef45940c9d8b3d7d670cd5ef224c356d48f8088f9b4892d4c7f164",
    "prepare_qfo_corrected_orthofinder_pairs.py": "c804a4070a8835c98f94b5079184de9123d7c2d179700ecfcc0bb6d4da602cd5",
    "prepare_qfo_corrected_fastoma_pairs.py": "c4437cd13a61573100d5a3fabdbb2fb3d15c6e2f2570f13b4b4fc60631953069",
    "prepare_qfo_corrected_orthomcl_pairs.py": "87b229469466cf330259db7f1ec5cd03609001bbb5d594ef49583bd4f832b23e",
    "publication_comparison.py": "99f55e93ab60b2220663b631980076c79fe99cb9569d8a6596f7afe12fcabb5b",
    "validate_qfo_native_assessment.py": "4786081a21ce04b8ae25f1f3e83ed8a1766472ff7e4476ae70752fe7f0656da5",
    "prepare_ob_candidate_neighborhood.py": "03843ed17fa9c44ea1ce5249782cdf009d8c7db9ceb85c62aa85bc2b4156ec29",
    "run_simulation_methods.py": "cccf579c221c51e8d03d512264efc3d0a61f4c0a598863d2423bdbd90f7921be",
}


def run(comparison, comparison_sha, baseline, protocol_path, release_protocol_path, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError("Require fresh corrected comparison result")
    source = record(__file__)
    helpers = [record(Path(__file__).with_name(name)) for name in SOURCES]
    if any(item["sha256"] != SOURCES[name] for name, item in zip(SOURCES, helpers)):
        raise ValueError("Changed frozen corrected comparison implementation")
    protocol, release_protocol = record(protocol_path), record(release_protocol_path)
    if protocol["sha256"] != PROTOCOL_SHA or release_protocol["sha256"] != RELEASE_PROTOCOL_SHA:
        raise ValueError("Changed corrected comparator protocol")
    checked = [source, *helpers, protocol, release_protocol, record(comparison), record(baseline)]
    for item in checked:
        check(item)
    counts = audit(comparison, comparison_sha, baseline)
    result = bootstrap(counts)
    if not result["protocol_controls_match"]:
        raise ValueError("Changed corrected comparator bootstrap controls")
    checked.extend(counts["checked_inputs"])
    for item in checked:
        check(item)
    estimated = sum(row["status"] == "estimated" for row in result["comparisons"])
    result.update(status="corrected_swiss_comparison_intervals_audited", source=source, helpers=helpers,
        protocol=protocol, release_protocol=release_protocol, comparison=record(comparison),
        reconstructed_counts=counts, checked_inputs=checked, scientific_inputs_admitted=True,
        estimated_contrasts=estimated, complete_panel=estimated == 8, uncertainty_admitted=estimated > 0)
    result["limitations"].append("Raw prediction counts and score arithmetic were reconstructed; upstream native inference/scoring admissions were not rerun.")
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("comparison", "baseline", "protocol", "release-protocol", "output"):
        parser.add_argument("--" + flag, type=Path, required=True)
    parser.add_argument("--comparison-sha256", required=True)
    args = parser.parse_args()
    run(args.comparison.resolve(), args.comparison_sha256, args.baseline.resolve(), args.protocol.resolve(),
        args.release_protocol.resolve(), args.output.absolute())
