"""Check strict FastOMA row conversion against retained original-release pairs."""

import argparse
import json
from itertools import zip_longest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.inventory_swiss_sequences import PREPARED_SHA
from benchmark_tools.fastoma_to_pairwise import input_owners, iter_pairs


def compare(native, retained, owners):
    count = 0
    with retained.open("rt", encoding="ascii", newline="") as stream:
        for pair, line in zip_longest(iter_pairs(native, owners), stream):
            if pair is None or line is None or line != "\t".join(pair) + "\n":
                raise ValueError(f"Retained FastOMA conversion differs at row {count + 1}")
            count += 1
    if not count:
        raise ValueError("Empty FastOMA pair output")
    return count


def audit(root, output):
    if output.exists():
        raise FileExistsError(output)
    prepared_path = root / "benchmark_tools/results/qfo_factorial_prepared_20260917.json"
    prepared = read_frozen(prepared_path, PREPARED_SHA)
    native = root / "qfo_benchmark/results/fastoma/output/orthologs.tsv.gz"
    retained = root / "qfo_benchmark/results/fastoma/pairs.tsv"
    records = [record(prepared_path), record(native), record(retained), *prepared["input_fastas"],
               record(Path(__file__).with_name("fastoma_to_pairwise.py"))]
    for item in records:
        check(item)
    owners = input_owners([Path(row["path"]) for row in prepared["input_fastas"]])
    count = compare(native, retained, owners)
    for item in records:
        check(item)
    result = {"status": "strict_fastoma_rows_match_retained_conversion", "source": record(__file__),
              "inputs": records, "input_accessions": len(owners), "validated_pair_rows": count,
              "corrected_release_result": False, "accuracy_evaluated": False,
              "limitations": ["Original-release conversion audit, not a corrected-release result.",
                              "Order and multiplicity preserved; global pair uniqueness is not asserted.",
                              "Native workflow completeness and biological orthology correctness are not established.",
                              "No retained predictions or scores were modified."]}
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    audit(args.root.resolve(), args.output.resolve())
