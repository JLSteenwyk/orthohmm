"""Read retained OrthoFinder OrthoBench input, checkpoint and command evidence."""

import argparse
import json
from pathlib import Path
import shlex

from Bio import SeqIO

from benchmark_tools.orthofinder_mcl_to_orthogroups import load_sequence_ids, iter_mcl_clusters
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.summarize_matched_resources import verbose_time


def fasta(path):
    result = {}
    for seq in SeqIO.parse(path, "fasta"):
        if seq.id in result:
            raise ValueError("Duplicate FASTA ID")
        result[seq.id] = str(seq.seq)
    if not result:
        raise ValueError("Empty FASTA")
    return result


def species_ids(path):
    result = {}
    for line in path.read_text().splitlines():
        key, separator, value = line.partition(": ")
        if not separator or not key.isdigit() or int(key) in result or Path(value).name != value:
            raise ValueError("Malformed/duplicate species mapping")
        result[int(key)] = value
    if not result or set(result) != set(range(len(result))) or len(set(result.values())) != len(result):
        raise ValueError("Incomplete species mapping")
    return result


def compare_species(original, processed, mapping, index):
    internal = {key: gene for key, gene in mapping.items() if key.split("_")[0] == str(index)}
    if len(set(internal.values())) != len(internal):
        raise ValueError("Duplicate original IDs in mapping")
    mapped = {gene: processed[key] for key, gene in internal.items() if key in processed}
    shared = original.keys() & mapped.keys()
    return dict(input_genes=len(original), processed_genes=len(processed),
                mapping_ids_match_processed=set(internal) == set(processed),
                original_ids_match=set(mapped) == set(original),
                sequence_mismatch_ids=sorted(g for g in shared if original[g] != mapped[g]))


def canonical(groups, universe):
    groups = [tuple(sorted(group)) for group in groups]
    flat = [g for group in groups for g in group]
    if any(not g for g in groups) or len(flat) != len(set(flat)) or set(flat) != set(universe):
        raise ValueError("Partition does not cover the exact gene universe")
    return sorted(groups)


def read_groups(path, labeled):
    groups, labels = [], set()
    for line in path.read_text().splitlines():
        fields = line.split()
        if labeled:
            if not fields or not fields[0].endswith(":") or fields[0] in labels:
                raise ValueError("Invalid/duplicate native group label")
            labels.add(fields[0])
            fields = fields[1:]
        if not fields:
            raise ValueError("Empty group row")
        groups.append(fields)
    return groups


def log_command(text, prefix):
    rows = [line.strip()[len(prefix):] for line in text.splitlines() if line.strip().startswith(prefix)]
    if len(rows) != 1:
        raise ValueError("Expected one command record")
    value = rows[0]
    if prefix == "Command being timed: ":
        if not value.startswith('"') or not value.endswith('"'):
            raise ValueError("Malformed GNU-time command")
        value = value[1:-1]
    return shlex.split(value)


def audit(root):
    results = root / "benchmark_tools/results"
    inventory_path = results / "installed_ob_input_inventory_20260926.json"
    comparison_path = results / "publication_comparison_orthomcl_complete_20260916.json"
    readback_path = results / "retained_ob_comparator_readback_20260926.json"
    pinned = [(inventory_path, "8b429398e6c381fcdec8289c641a7f20969edcad5595e3d3f6930335601571bd"),
              (comparison_path, "094f842ad8211450519d4238c0b3a5020a7969049d4b7306f5465d5acbf914a3"),
              (readback_path, "5dcb4e65f277eb9debf2bbacb1d4e298c678ac13fc981e75c4ba78e4851eeb55")]
    checked = []
    for path, sha in pinned:
        item = record(path)
        if item["sha256"] != sha:
            raise ValueError("Changed pinned evidence")
        checked.append(item)
    inventory, comparison, readback = [json.loads(p.read_text()) for p, _ in pinned]
    base = root / "benchmarks/results/orthofinder_v3_diamond"
    native = base / "input/OrthoFinder/Results_Jul15"
    working = native / "WorkingDirectory"
    sequence = root / "benchmarks/results/orthofinder_v3_sequence_only_20260830"
    paths = [working / name for name in ("SpeciesIDs.txt", "SequenceIDs.txt", "clusters_OrthoFinder_I1.2.txt_id_pairs.txt")]
    paths += [native / "Log.txt", base / "run.log", base / "time.log", base / "tool_version.txt",
              base / "orthogroups_source.txt", sequence / "conversion_time.log"]
    checked.extend(record(p) for p in paths)
    species = species_ids(working / "SpeciesIDs.txt")
    mapping = load_sequence_ids(working / "SequenceIDs.txt")
    if len(set(mapping.values())) != len(mapping):
        raise ValueError("Duplicate mapped gene IDs")
    expected_inputs = {Path(r["path"]).name: r for r in inventory["inputs"]}
    if set(species.values()) != set(expected_inputs) or len(mapping) != inventory["genes"]:
        raise ValueError("Species/gene universe differs from frozen input inventory")
    rows = []
    seen = set()
    for index, name in species.items():
        staged, processed = base / "input" / name, working / f"Species{index}.fa"
        source = expected_inputs[name]
        check(source)
        a, b = record(staged), record(processed)
        checked.extend([source, a, b])
        original = fasta(staged)
        if seen & original.keys():
            raise ValueError("Gene ID reused across species")
        seen.update(original)
        row = compare_species(original, fasta(processed), mapping, index)
        row.update(species=name, staged=a, processed=b, reference_input=source,
                   staged_bytes_match=all(a[k] == source[k] for k in ("bytes", "sha256")))
        rows.append(row)
    full_record = next(r for r in comparison["methods"] if r["key"] == "orthofinder_3_1_5_full")["orthobench"]["prediction_provenance"]
    sequence_record = next(r for r in readback["rows"] if r["key"] == "orthofinder_3_1_5_sequence_only")["prediction"]
    for item in (full_record, sequence_record):
        check(item)
        checked.append(item)
    expected_full = native / "Orthogroups/Orthogroups.txt"
    if full_record != record(expected_full) or sequence_record != record(sequence / "orthogroups.txt"):
        raise ValueError("Wrong retained prediction binding")
    full_groups = canonical(read_groups(expected_full, True), seen)
    checkpoint_groups = canonical(([mapping[g] for g in group] for group in iter_mcl_clusters(paths[2])), seen)
    sequence_groups = canonical(read_groups(sequence / "orthogroups.txt", False), seen)
    native_log, time_log = (native / "Log.txt").read_text(), (base / "time.log").read_text()
    command = log_command(native_log, "Command Line: ")
    timed_command = log_command(time_log, "Command being timed: ")
    conversion_log = (sequence / "conversion_time.log").read_text()
    conversion = log_command(conversion_log, "Command being timed: ")
    expected_conversion = ["python", "benchmark_tools/orthofinder_mcl_to_orthogroups.py",
                           str(paths[2].relative_to(root)), str(paths[1].relative_to(root)),
                           str((sequence / "orthogroups.txt").relative_to(root))]
    gates = dict(staged_inputs_match=all(r["staged_bytes_match"] for r in rows),
                 processed_sequences_match=all(r["mapping_ids_match_processed"] and r["original_ids_match"] and not r["sequence_mismatch_ids"] for r in rows),
                 checkpoint_conversion_matches=checkpoint_groups == sequence_groups,
                 full_command_matches=command == timed_command,
                 expected_settings=command[1:] == ["-f", str(base / "input"), "-t", "32", "-a", "8", "-S", "diamond"],
                 conversion_command_matches=conversion == expected_conversion,
                 logged_version_and_completion="Started OrthoFinder version 3.1.5" in native_log and "OrthoFinder run completed" in native_log
                     and (base / "tool_version.txt").read_text().strip() == "OrthoFinder:v3.1.5")
    checked += [record(Path(__file__).with_name(name)) for name in
                ("audit_ob_orthofinder_provenance.py", "orthofinder_mcl_to_orthogroups.py", "summarize_matched_resources.py", "gnu_time_companion.py")]
    for item in checked:
        check(item)
    return dict(status="retained_orthofinder_ob_provenance_readback", all_checks_passed=all(gates.values()), checks=gates,
                species=rows, genes=len(seen), full_groups=len(full_groups), sequence_groups=len(sequence_groups),
                command=command, checkpoint_conversion_command=conversion,
                full_inference_resources=verbose_time(time_log), checkpoint_conversion_resources=verbose_time(conversion_log),
                sequence_only_inference_resources=None, checked_records=checked, publication_ready=False,
                limitations=["Retrospective retained-file consistency, not independently attested historical execution or binary identity.",
                             "Checkpoint readback shares the existing converter parser; not an independent MCL parser validation.",
                             "Sequence-only is a checkpoint from the full run; its conversion time is not inference runtime.",
                             "GNU-time resources are descriptive historical process measurements, not controlled comparative or process-tree memory evidence.",
                             "Full output semantics rely on the separately documented version-specific source audit, not filename inference.",
                             "No native inference or benchmark scoring rerun; no changes to retained scores."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.root.resolve())
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    raise SystemExit(0 if result["all_checks_passed"] else 1)
