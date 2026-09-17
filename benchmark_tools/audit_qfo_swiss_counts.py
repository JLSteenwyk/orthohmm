"""Recover and verify SwissTrees family counts before uncertainty analysis."""

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.report_qfo_recovered_stages import ADMISSION_SHA, STAGES
from benchmark_tools.validate_qfo_native_assessment import validate_records

SCORER_SHA = "262c8d1f06527461e11391682ccaccc637675f84a5064c1658b80bf1b956bd51"
IMAGE_SHA = "abed2b0ff8bb033a2372d0cd5282fd1d6c8ab5fa07496afc48c82f54ed65ea21"
REFERENCE_SHA = "40528e3537ba57e345c443c8024f99dc3940a66e2d102703754b41bb29858ecd"
HEADER = "# Dataset<tab>Protein ID 1<tab>Protein ID 2<tab>Correctness (TP:True positive, FP: False positive, TN: True Negative. FN: False Negative)"
LABELS = ("TP", "FP", "FN", "TN")


def read_raw(path, families):
    counts = {family: Counter({label: 0 for label in LABELS}) for family in families}
    truth, members = {}, {family: set() for family in families}
    with gzip.open(path, "rt", newline="") as handle:
        if handle.readline().rstrip("\r\n") != HEADER:
            raise ValueError("Unexpected native raw header")
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) != 4:
                raise ValueError("Invalid raw row width")
            family, a, b, label = row
            if family not in counts or not a or not b or a == b or label not in LABELS:
                raise ValueError("Invalid raw family/pair/label")
            key = (family, *sorted((a, b)))
            if key in truth:
                raise ValueError("Duplicate family pair")
            truth[key] = label in ("TP", "FN")
            members[family].update((a, b))
            counts[family][label] += 1
    if any(not sum(c.values()) for c in counts.values()):
        raise ValueError("Empty reference family")
    return counts, truth, members


def statistics(counts):
    # Each stored relation contributes 1/2, then native code adds a prior of 1.
    tp, fp, fn = (counts[name] / 2 + 1 for name in ("TP", "FP", "FN"))
    precision, recall = tp / (tp + fp), tp / (tp + fn)
    return {"PPV": precision, "TPR": recall,
            "F1": 2 * precision * recall / (precision + recall)}


def verify_native(counts, native):
    values = {family: statistics(row) for family, row in counts.items()}
    for family, scores in values.items():
        for metric in ("TPR", "PPV"):
            expected = native["SwissTrees-" + family, metric]["metrics"]["value"]
            if not math.isclose(scores[metric], expected, rel_tol=0, abs_tol=5e-8):
                raise ValueError("Family counts do not reproduce native metric")
    means = {metric: sum(v[metric] for v in values.values()) / len(values) for metric in ("TPR", "PPV")}
    for metric, value in means.items():
        if not math.isclose(value, native["SwissTrees", metric]["metrics"]["value"], rel_tol=0, abs_tol=5e-8):
            raise ValueError("Family macro mean differs from native aggregate")
    means["F1"] = 2 * means["TPR"] * means["PPV"] / (means["TPR"] + means["PPV"])
    return values, means


def reference_orientation(image, reference):
    if "'" in str(reference):
        raise ValueError("Unsupported Darwin path quoting")
    script = Path(__file__).with_name("audit_qfo_swiss_reference.drw")
    program = f"reference := '{reference}':\n" + script.read_text()
    result = subprocess.run(["singularity", "exec", str(image), "darwin", "-E"],
                            input=program, capture_output=True, text=True, timeout=60, check=True)
    rows = re.findall(r"^SWISS_REFERENCE\t([^\t\n]+)\t(\d+)\t(\d+)\t(\d+)$", result.stdout, re.MULTILINE)
    if len(rows) != 18 or len({r[0] for r in rows}) != 18:
        raise ValueError("Darwin reference inventory incomplete")
    inventory = {name: {"mapped_proteins": int(n), "forward_relations": int(f), "other_relations": int(r)}
                 for name, n, f, r in rows}
    if any(r["other_relations"] != 0 or r["mapped_proteins"] <= 5 for r in inventory.values()):
        raise ValueError("Native one-direction counting assumption failed")
    container_source = subprocess.run(["singularity", "exec", str(image), "cat", "/benchmark/RefPhyloTest.drw"],
                                      capture_output=True, timeout=30, check=True).stdout
    if hashlib.sha256(container_source).hexdigest() != SCORER_SHA:
        raise ValueError("Container scorer differs from pinned inspected scorer")
    return inventory, {"source": record(script), "stdout": result.stdout, "stderr": result.stderr,
                       "container_scorer_sha256": SCORER_SHA}


def audit(admission_path, raw_root, scorer, image, reference):
    admission_record, scorer_record = record(admission_path), record(scorer)
    if admission_record["sha256"] != ADMISSION_SHA or scorer_record["sha256"] != SCORER_SHA:
        raise ValueError("Changed frozen admission or native scorer")
    image_record, reference_record = record(image), record(reference)
    if image_record["sha256"] != IMAGE_SHA or reference_record["sha256"] != REFERENCE_SHA:
        raise ValueError("Changed native image or reference")
    orientation, orientation_evidence = reference_orientation(image, reference)
    admission = json.loads(admission_path.read_text())
    if admission["status"] != "four_stage_assessments_checked" or [r["stage"] for r in admission["records"]] != list(STAGES):
        raise ValueError("Changed stage admission")
    records, expected_truth, expected_members, expected_families = [], None, None, None
    for index, row in enumerate(admission["records"]):
        if row["status"] != "admitted" or row["index"] != index:
            raise ValueError("Unadmitted stage")
        families = row["assessment"]["swiss_reference_families"]
        if len(families) != 18 or len(set(families)) != 18:
            raise ValueError("Wrong reference inventory")
        if expected_families is None:
            expected_families = families
        if families != expected_families:
            raise ValueError("Stage family inventories differ")
        native = validate_records(row["assessment"]["native_assessments"], row["participant"], set(families))
        paths = list((raw_root / f"checked_v2_{index}/results/SwissTrees").glob("*raw.txt.gz"))
        if len(paths) != 1:
            raise ValueError("Ambiguous raw evidence")
        raw_record = record(paths[0])
        counts, truth, members = read_raw(paths[0], families)
        if set(orientation) != set(families) or any(
                sum(counts[f].values()) != orientation[f]["forward_relations"] or
                len(members[f]) != orientation[f]["mapped_proteins"] for f in families):
            raise ValueError("Raw evidence does not cover the native reference inventory")
        values, aggregate = verify_native(counts, native)
        if expected_truth is None:
            expected_truth, expected_members = truth, members
        if truth != expected_truth or members != expected_members:
            raise ValueError("Stage reference relations or members differ")
        check(raw_record)
        records.append({"stage": row["stage"], "raw_file": raw_record, "aggregate": aggregate,
                        "families": [{"family": family, "counts_without_prior": dict(counts[family]),
                                      "statistics_with_prior": values[family],
                                      "represented_genes": sorted(members[family])} for family in families]})
    owners = defaultdict(list)
    for family, genes in expected_members.items():
        for gene in genes:
            owners[gene].append(family)
    overlaps = {gene: sorted(names) for gene, names in sorted(owners.items()) if len(names) > 1}
    for row in records:
        check(row["raw_file"])
    check(admission_record)
    check(scorer_record)
    check(image_record)
    check(reference_record)
    return {"status": "raw_swiss_family_counts_verified", "admission": admission_record,
            "native_scorer": scorer_record, "source": record(Path(__file__)), "stages": records,
            "families": expected_families, "reference_relation_count": len(expected_truth),
            "reference_orientation": orientation, "orientation_evidence": orientation_evidence,
            "container": image_record, "reference": reference_record,
            "count_conversion": "native confusion count = raw one-direction relation count / 2 + 1",
            "shared_represented_genes": overlaps, "publication_ready": False,
            "limitations": ["No bootstrap intervals or significance established by count validation.",
                            "Finite decimal native output compared with absolute tolerance 5e-8.",
                            "Shared evolutionary history and merged predictions can correlate distinct families.",
                            "Raw artifact hashes are current snapshots, not retroactive execution traces.",
                            "Only SwissTrees is covered; other QfO challenges need separate uncertainty methods."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("admission", "raw-root", "scorer", "image", "reference", "output"):
        parser.add_argument("--" + flag, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = audit(args.admission.resolve(), args.raw_root.resolve(), args.scorer.resolve(),
                   args.image.resolve(), args.reference.resolve())
    with args.output.open("x") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
