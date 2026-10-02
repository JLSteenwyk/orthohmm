"""Check retained SwissTrees table arithmetic, not raw annotation admission."""

import argparse
import csv
from fractions import Fraction
import hashlib
import io
import json
from pathlib import Path

METRICS = ("F1", "PPV", "TPR")
REFERENCE = "orthofinder_3_1_5_full"
COUNTS = "qfo_recovered_swiss_uncertainty_22178.json"
METHODS = {"orthohmm_high_sensitivity", "orthohmm_phylogeny_satellite_v2", REFERENCE,
           "orthofinder_3_1_5_sequence_only", "sonicparanoid_2_0_9",
           "proteinortho_6_3_6", "fastoma_0_3_5", "orthomcl_1_4"}
PANELS = {
    "descriptive": ("corrected_swiss_sequence_strata_20260918.json",
                    "e8cc9ae62434192f33e56c70f6d54d5ccc49b51ef21c9b23789940ea04e51fc9"),
    "identity": ("corrected_swiss_identity_admission_22102.json",
                 "8d1684de6cd1aa98fdefdd47de78d1b9eebde47d8a5e205460891b982a54195e"),
    "fragment": ("swiss_historical_fragment_admission_22117.json",
                 "c93c436ea899cae3a5b95857ecef5c0e746ada9d8833e7283f467227c1c93862"),
    "duplication": ("swiss_duplication_features_v2_20260923.json",
                    "2269e83e9a8bdd80c5286be773ca9eb83234908dfa4e223ae6769a3d7e6fa991"),
}


def record(path, data):
    return dict(path=str(path.resolve()), bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def read_pinned(path, expected):
    data = path.read_bytes()
    actual = record(path, data)
    if (actual["sha256"] != expected["sha256"]
            or ("bytes" in expected and actual["bytes"] != expected["bytes"])):
        raise ValueError("Changed retained bytes: " + path.name)
    return data, actual


def number_matches(value, expected):
    if expected is None:
        if value is not None:
            raise ValueError("Missing score was imputed")
        return
    if type(value) not in (int, float, str):
        raise ValueError("Invalid numeric score")
    try:
        difference = abs(Fraction(str(value)) - expected)
    except (ValueError, ZeroDivisionError):
        raise ValueError("Nonfinite or invalid numeric score") from None
    if difference > Fraction("1e-12"):
        raise ValueError("Rational score differs")


def aggregate(values):
    p = sum((v[0] for v in values), Fraction()) / len(values)
    r = sum((v[1] for v in values), Fraction()) / len(values)
    return dict(F1=2*p*r/(p+r), PPV=p, TPR=r)


def rational_rows(counts, bins, differences):
    if (counts["status"] != "corrected_swiss_comparison_intervals_audited"
            or counts["scientific_inputs_admitted"] is not True):
        raise ValueError("Require admitted sufficient counts")
    families = counts["families"]
    if (not families or len(set(families)) != len(families)
            or bins.get("all") != families):
        raise ValueError("Invalid family universe or all-family bin")
    for members in bins.values():
        if len(set(members)) != len(members) or not set(members) <= set(families):
            raise ValueError("Invalid stratum membership")
    methods = counts["reconstructed_counts"]["methods"]
    names = [m["method"] for m in methods]
    if not methods or len(set(names)) != len(names) or set(names) != set(counts["point_estimates"]):
        raise ValueError("Invalid method inventory")
    if differences and REFERENCE not in names:
        raise ValueError("Missing reference method")
    rows = []
    for method in methods:
        values = {}
        if method["status"] == "counts_verified":
            if [r["family"] for r in method["families"]] != families:
                raise ValueError("Incomplete or reordered method family counts")
            for family in method["families"]:
                raw = family["counts_without_prior"]
                if (set(raw) != {"TP", "FP", "FN", "TN"}
                        or any(type(v) is not int or v < 0 for v in raw.values())):
                    raise ValueError("Invalid family sufficient counts")
                # Native directed counts are halved, then one is added to each
                # confusion category: these equivalent ratios stay exact.
                tp, fp, fn = (raw[k] for k in ("TP", "FP", "FN"))
                pair = (Fraction(tp+2, tp+fp+4), Fraction(tp+2, tp+fn+4))
                values[family["family"]] = pair
                for key, value in aggregate([pair]).items():
                    number_matches(family["statistics_with_prior"][key], value)
            for key, value in aggregate(list(values.values())).items():
                number_matches(counts["point_estimates"][method["method"]][key], value)
        elif method["status"] != "not_admitted" or counts["point_estimates"][method["method"]] is not None:
            raise ValueError("Invalid method admission or imputed estimates")
        for name, members in bins.items():
            available = bool(members) and bool(values)
            score = aggregate([values[f] for f in members]) if available else dict.fromkeys(METRICS)
            rows.append(dict(method=method["method"], stratum=name, families=len(members),
                family_members=members, prediction_semantics=method.get("prediction_semantics"),
                status="descriptive" if available else "empty_bin" if values else "method_not_admitted",
                **score))
    if differences:
        reference = {r["stratum"]: r for r in rows if r["method"] == REFERENCE}
        for row in rows:
            for metric in METRICS:
                left, right = row[metric], reference[row["stratum"]][metric]
                row["delta_" + metric] = None if left is None or right is None else left-right
    return rows


def panel_bins(kind, feature, families):
    if feature["prediction_statistics_evaluated"] is not False or feature["publication_ready"] is not False:
        raise ValueError("Require unscored descriptive features")
    members = feature["families"] if kind in ("fragment", "duplication") else feature["family_memberships"]
    if set(members) != set(families):
        raise ValueError("Feature family universe differs")
    bins = {"all": families}
    if kind == "fragment":
        sources = [{prefix + key: value for key, value in feature[field].items()}
                   for prefix, field in (("historical_", "strata"), ("baseline_only_", "baseline_only_strata"))]
    elif kind == "descriptive":
        sources = [feature["primary_strata"], feature["secondary_strata"]]
    else:
        sources = [feature["strata" if kind == "identity" else "primary_strata"]]
    for source in sources:
        if set(bins) & set(source):
            raise ValueError("Duplicate stratum names")
        bins.update(source)
    return bins


def compare_rows(retained, expected, *, tsv=False):
    indexed = {(r["method"], r["stratum"]): r for r in retained}
    keys = {(r["method"], r["stratum"]) for r in expected}
    if len(indexed) != len(retained) or set(indexed) != keys:
        raise ValueError("Table row inventory differs")
    for row in expected:
        actual = indexed[row["method"], row["stratum"]]
        fields = set(row) - ({"family_members"} if tsv else set())
        if set(actual) != fields:
            raise ValueError("Table field inventory differs")
        for key in fields:
            value = actual[key]
            if key in METRICS or key.startswith("delta_"):
                number_matches(None if tsv and value == "NA" else value, row[key])
            elif not tsv and type(value) is not type(row[key]):
                raise ValueError("Table metadata type differs: " + key)
            elif value != (("NA" if row[key] is None else str(row[key])) if tsv else row[key]):
                raise ValueError("Table metadata differs: " + key)


def verify(base):
    checked, panels = [], []
    for kind, (feature_name, manifest_sha) in PANELS.items():
        directory = base / f"swiss_{kind}_strata_20260926"
        raw, pin = read_pinned(directory / "manifest.json", {"sha256": manifest_sha})
        checked.append(pin)
        manifest = json.loads(raw)
        if manifest["new_inferential_claims"] is not False or manifest["publication_ready"] is not False:
            raise ValueError("Unexpected retained claim flags")
        inputs = []
        for name, source in zip((COUNTS, feature_name), manifest["inputs"][:2]):
            if Path(source["path"]).name != name:
                raise ValueError("Unexpected retained input identity")
            raw, pin = read_pinned(base / name, source)
            checked.append(pin)
            inputs.append(json.loads(raw))
        counts, features = inputs
        if len(counts["families"]) != 18 or {r["method"] for r in counts["reconstructed_counts"]["methods"]} != METHODS:
            raise ValueError("Frozen eighteen-family/eight-method panel differs")
        expected = rational_rows(counts, panel_bins(kind, features, counts["families"]), kind != "descriptive")
        compare_rows(manifest["rows"], expected)
        outputs = {Path(r["path"]).name: r for r in manifest["outputs"]}
        if len(outputs) != len(manifest["outputs"]) or set(outputs) != {"scores.tsv", "scores.md"}:
            raise ValueError("Retained output inventory differs")
        for name, source in outputs.items():
            raw, pin = read_pinned(directory / name, source)
            checked.append(pin)
            if name == "scores.tsv":
                reader = csv.DictReader(io.StringIO(raw.decode("utf-8")), delimiter="\t")
                if len(set(reader.fieldnames)) != len(reader.fieldnames):
                    raise ValueError("Duplicate TSV columns")
                table = list(reader)
                compare_rows(table, expected, tsv=True)
        metrics = [k for k in expected[0] if k in METRICS or k.startswith("delta_")]
        panels.append(dict(kind=kind, rows=len(expected), methods=8, families=18,
            score_and_difference_cells=len(expected)*len(metrics),
            missing_cells=sum(row[k] is None for row in expected for k in metrics),
            strata=list(dict.fromkeys(r["stratum"] for r in expected))))
    for pin in checked:
        read_pinned(Path(pin["path"]), pin)
    return dict(status="swiss_descriptive_table_arithmetic_verified", panels=panels,
        rows=sum(p["rows"] for p in panels),
        score_and_difference_cells=sum(p["score_and_difference_cells"] for p in panels),
        absolute_tolerance="1e-12", checked_records=checked,
        checker=record(Path(__file__), Path(__file__).read_bytes()), publication_ready=False,
        raw_inputs_revalidated=False, bootstrap_intervals_recomputed=False,
        limitations=["Retained derived counts and bin assignments only; not independent biological validation.",
            "No annotation, alignment, tree, source reference or raw-QfO scoring regeneration/admission.",
            "Exact rational arithmetic checks JSON/TSV values and metadata; Markdown bytes are pinned, not regenerated.",
            "Historical absolute paths are provenance and are never opened; no source/data rights clearance.",
            "Development-exposed descriptive associations are not causal or confirmatory claims."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    result = verify(args.results.resolve())
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
