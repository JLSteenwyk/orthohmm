"""Corrected-release SwissTrees input descriptors and unscored composition bins."""

import argparse
import json
import math
from pathlib import Path
import statistics
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.inventory_swiss_sequences import collect
from benchmark_tools.audit_qfo_input_sequences import sequence_identity
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from Bio import SeqIO

STAGING_SHA = "07a890eb816944f946d46039559a6046c2a3664f033eaa0f30f35bb25b9c9ab8"
COUNTS_SHA = "2995868b0407ceda2e7422db6c3fc3b99523716e4bec8b1b43e15f296e40769b"
ORIGINAL_SHA = "f4a147ff1ecdf0c8df046ab1a5c63069fef106fa627d52eaec008aae250527a5"
HELPER_SHA = "7f6dd0d3ada54ea2d34c091bc33a667e688518acf9405c8fbffaef7d949d5d0d"
PROTOCOL_SHA = "2c87ba1df8ee39dfc325eefcdf374e713ba8fab1279a0330c3ea600ab637da54"
UPDATE_SHA = "69b20dd4b656e9c6f7dbbb84bdc228d1919389a2554d389f790c59df09d8b4cd"
NATIVE_AUDIT_SHA = "480a09d0c60c95274d93fb49e15ae95b45a99d0cdcd17c72bf9b3af066cc12c4"
CHANGED_SEQUENCES = {"F6PXR2", "F7BEF0", "F7C949"}


def descriptor_updates(old, new):
    if len(old) != 549 or not set(old) <= set(new):
        raise ValueError("Old descriptor coverage changed")
    if any(set(row) != set(new[g]) for g, row in old.items()):
        raise ValueError("Descriptor schema changed")
    changed = {g: {k: dict(old=row[k], new=new[g][k]) for k in row if row[k] != new[g][k]}
               for g, row in old.items() if row != new[g]}
    if set(changed) != CHANGED_SEQUENCES | {"A0A6I8Q293"}:
        raise ValueError("Unexpected cross-release descriptor changes")
    expected = {"canonical_entropy_bits", "canonical_residues", "description", "length", "maximum_canonical_frequency"}
    if any(set(changed[g]) != expected for g in CHANGED_SEQUENCES):
        raise ValueError("Unexpected changed descriptor fields")
    header = changed["A0A6I8Q293"]
    if set(header) != {"description"} or header["description"]["old"].replace(" PE=4 ", " PE=3 ") != header["description"]["new"]:
        raise ValueError("Unexpected header-only update")
    return changed


def define_strata(families, genes):
    if set(genes) != {gene for members in families.values() for gene in members}:
        raise ValueError("Incomplete or extra feature coverage")
    features = {}
    for family, members in sorted(families.items()):
        if not members or len(set(members)) != len(members):
            raise ValueError("Empty or duplicate family members")
        rows = [genes[g] for g in members]
        for row in rows:
            n, length, entropy = row["canonical_residues"], row["length"], row["canonical_entropy_bits"]
            if type(length) is not int or type(n) is not int or not 0 <= n <= length or length <= 0:
                raise ValueError("Invalid protein lengths")
            if (n == 0) != (entropy is None):
                raise ValueError("Entropy missingness differs from canonical count")
            if entropy is not None and (isinstance(entropy, bool) or not math.isfinite(entropy)
                                       or not 0 <= entropy <= math.log2(20) + 1e-12):
                raise ValueError("Invalid canonical entropy")
        eligible = [r["canonical_residues"] >= 20 and r["canonical_residues"] / r["length"] >= .9 for r in rows]
        median = statistics.median(r["canonical_entropy_bits"] / math.log2(20) for r in rows) if all(eligible) else None
        concentrated = [r["canonical_entropy_bits"] / math.log2(20) < .8 if ok else None
                        for r, ok in zip(rows, eligible)]
        composition = ("concentrated" if any(v is True for v in concentrated) else
                       "missing" if any(v is None for v in concentrated) else "not_concentrated")
        lengths = [r["length"] for r in rows]
        features[family] = dict(genes=len(rows), median_normalized_entropy=median,
            composition=composition, median_length=statistics.median(lengths), minimum_length=min(lengths),
            short_relative=min(lengths) < .5 * statistics.median(lengths),
            explicit_fragment_descriptions=sum(r["explicit_fragment_description"] for r in rows))
    values = [r["median_normalized_entropy"] for r in features.values() if r["median_normalized_entropy"] is not None]
    cutoff = statistics.median(values) if values else None
    primary = {key: [] for key in ("lower_entropy", "higher_entropy", "missing_entropy")}
    secondary = {key: [] for key in ("concentrated", "not_concentrated", "missing", "short_relative",
                                     "not_short_relative", "explicit_fragment", "no_explicit_fragment")}
    for family, row in features.items():
        entropy = row["median_normalized_entropy"]
        primary["missing_entropy" if entropy is None else "lower_entropy" if entropy <= cutoff else "higher_entropy"].append(family)
        secondary[row["composition"]].append(family)
        secondary["short_relative" if row["short_relative"] else "not_short_relative"].append(family)
        secondary["explicit_fragment" if row["explicit_fragment_descriptions"] else "no_explicit_fragment"].append(family)
    return dict(family_features=features, primary_strata=primary, secondary_strata=secondary,
                median_family_entropy_cutoff=cutoff)


def prepare(root, output):
    if output.exists():
        raise FileExistsError(output)
    results = root / "benchmark_tools/results"
    paths = [results / name for name in ("qfo_corrected_staging_manifest_20260918.json",
             "qfo_swiss_comparator_counts_20260917.json", "swiss_sequence_inventory_20260917.json")]
    identities = [record(p) for p in paths]
    if [r["sha256"] for r in identities] != [STAGING_SHA, COUNTS_SHA, ORIGINAL_SHA]:
        raise ValueError("Frozen input identities differ")
    helper = record(Path(__file__).with_name("inventory_swiss_sequences.py"))
    if helper["sha256"] != HELPER_SHA:
        raise ValueError("Feature implementation changed")
    protocol = record(results / "CORRECTED_SWISS_SEQUENCE_STRATA_PROTOCOL_20260918.md")
    if protocol["sha256"] != PROTOCOL_SHA:
        raise ValueError("Feature and outcome protocol changed")
    update = record(results / "CORRECTED_SWISS_DESCRIPTOR_UPDATE_20260918.md")
    native_audit = record(results / "qfo_original_input_sequence_audit_20260917.json")
    if [update["sha256"], native_audit["sha256"]] != [UPDATE_SHA, NATIVE_AUDIT_SHA]:
        raise ValueError("Descriptor update evidence changed")
    stage, counts, old = [json.loads(p.read_text()) for p in paths]
    if stage["total_sequences"] != 984137 or len(stage["input_fastas"]) != 78 or stage["inputs_normalized"]:
        raise ValueError("Unexpected corrected input inventory")
    # Only membership is used from the historical count file, never its scores.
    families = {r["family"]: r["represented_genes"] for r in counts["methods"][0]["families"]}
    if len(families) != 18 or sum(map(len, families.values())) != 563:
        raise ValueError("Changed reference universe")
    inventory = collect(families, stage["input_fastas"])
    if inventory["summary"]["matched_genes"] != 563 or inventory["summary"]["missing_genes"]:
        raise ValueError("Corrected release still lacks reference sequences")
    changes = descriptor_updates(old["genes"], inventory["genes"])
    expected = {r["accession"]: r["native_sequence"] for r in json.loads(Path(native_audit["path"]).read_text())["differences"]
                if r["accession"] in CHANGED_SEQUENCES}
    xenopus = next(r for r in stage["input_fastas"] if Path(r["path"]).name == "UP000008143_8364.fasta")
    check(xenopus)
    observed = {entry.id.split("|")[1]: sequence_identity(str(entry.seq)) for entry in SeqIO.parse(xenopus["path"], "fasta")
                if entry.id.split("|")[1] in CHANGED_SEQUENCES}
    check(xenopus)
    if set(expected) != CHANGED_SEQUENCES or observed != expected:
        raise ValueError("Corrected changed sequences differ from native database")
    strata = define_strata(families, inventory["genes"])
    for item in [*identities, helper, protocol, update, native_audit]:
        check(item)
    result = dict(status="corrected_swiss_sequence_strata_prepared_unscored", **inventory, **strata,
        recovered_accessions=sorted(set(inventory["genes"]) - set(old["genes"])), family_memberships=families,
        source=record(__file__), helper=helper, inputs=identities, fasta_inputs=stage["input_fastas"], protocol=protocol,
        descriptor_update_protocol=update, native_sequence_audit=native_audit,
        unchanged_original_descriptors=545, changed_original_descriptors=changes,
        changed_sequence_native_identities=observed,
        prediction_statistics_evaluated=False, publication_ready=False,
        limitations=["Input-only features on development-exposed curated families, not independent confirmation.",
            "No historical score is relabeled as a corrected-release result.",
            "Global entropy is not local low complexity; family entropy is not evolutionary divergence.",
            "Short-relative length is not a fragment diagnosis; absent fragment text does not prove completeness.",
            "Strata overlap and can differ in size, taxa, domain architecture and evolutionary history.",
            "No corrected accuracy contrast or uncertainty result is evaluated here."])
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.root.resolve(), args.output.absolute())
