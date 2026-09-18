"""Apply the frozen paired-family protocol to corrected factorial counts."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.bootstrap_qfo_factorial import bootstrap, markdown, PROTOCOL_SHA
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen

CORRECTED_PROTOCOL_SHA = "b3603bc9b51ce1708f49a0ee816eb02d3cb52c66fb00f07eeb0d66f840a03e1e"
ENGINE_SHA = "be09876ae4de31b923818385bd3d40d8d5e7df177215fc45a2101c7d30d1c1fa"
AUDITOR_SHA = "375f2b4fd511c04e562827aebacb3b7522ca979f76b665652dd209940542f912"


def calculate(data, replicates=100000, seed=20260922):
    if (data.get("status") != "corrected_qfo_factorial_swiss_counts_verified"
            or data.get("publication_ready") is not False or data.get("uncertainty_admitted") is not False):
        raise ValueError("Require verified unresampled corrected factorial counts")
    # Only the schema tag changes for the shared numerical engine; counts are not replaced.
    adapted = {**data, "status": "qfo_factorial_swiss_counts_verified"}
    result = bootstrap(adapted, replicates=replicates, seed=seed)
    result["status"] = "paired_corrected_qfo_factorial_swiss_intervals"
    result["input_release"] = "QfO 2020_04 corrected UP000008143"
    result["limitations"].append("Corrected-release counts only; no substitution or pooling of historical prediction counts.")
    return result


def run(counts_path, counts_sha, protocol_path, corrected_protocol_path, output, markdown_path):
    if output == markdown_path or any(p.exists() or p.is_symlink() for p in (output, markdown_path)):
        raise FileExistsError("Require distinct fresh result paths")
    counts = record(counts_path)
    protocol, corrected_protocol = record(protocol_path), record(corrected_protocol_path)
    engine = record(Path(__file__).with_name("bootstrap_qfo_factorial.py"))
    auditor = record(Path(__file__).with_name("audit_qfo_corrected_factorial_swiss.py"))
    if (protocol["sha256"] != PROTOCOL_SHA or corrected_protocol["sha256"] != CORRECTED_PROTOCOL_SHA
            or engine["sha256"] != ENGINE_SHA or auditor["sha256"] != AUDITOR_SHA):
        raise ValueError("Frozen protocol or implementation changed")
    data = read_frozen(counts_path, counts_sha)
    if data["source"] != auditor:
        raise ValueError("Counts do not identify the corrected count auditor")
    checked = [counts, protocol, corrected_protocol, engine, auditor,
               data["admission_inventory"], data["baseline_audit"], *data["checked_inputs"], *data["helpers"]]
    for item in checked:
        check(item)
    result = calculate(data)
    for item in checked:
        check(item)
    result.update(counts=counts, protocol=protocol, corrected_protocol=corrected_protocol,
                  checked_inputs=checked, source=record(__file__),
                  helpers=[engine, *[record(Path(__file__).with_name(n)) for n in (
                      "audit_qfo_swiss_counts.py", "bootstrap_qfo_swiss_stages.py")]])
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    with markdown_path.open("x") as stream:
        stream.write(markdown(result).replace("# QfO Factorial SwissTrees Intervals",
                                             "# Corrected-Release QfO Factorial SwissTrees Intervals", 1))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("counts", "protocol", "corrected-protocol", "output", "markdown"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--counts-sha256", required=True)
    args = parser.parse_args()
    run(args.counts.resolve(), args.counts_sha256, args.protocol.resolve(), args.corrected_protocol.resolve(),
        args.output.absolute(), args.markdown.absolute())
