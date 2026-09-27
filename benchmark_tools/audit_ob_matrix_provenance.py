"""Bind retained SonicParanoid/Proteinortho OrthoBench groups to native tables."""

import argparse
import csv
import json
from pathlib import Path
import shlex

from benchmark_tools.audit_ob_orthofinder_provenance import canonical, fasta, read_groups
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


SCHEMAS = {"sonicparanoid": ("group_id", "group_size", "sp_in_grp", "seed_ortholog_cnt"),
           "proteinortho": ("# Species", "Genes", "Alg.-Conn.")}


def native_groups(path, method, owners):
    metadata = SCHEMAS[method]
    groups, seen, labels = [], set(), set()
    with path.open(newline="") as stream:
        reader = csv.reader(stream, delimiter="\t")
        header = next(reader)
        species = header[len(metadata):]
        if tuple(header[:len(metadata)]) != metadata or len(set(species)) != len(species) or set(species) != set(owners.values()):
            raise ValueError("Native species-column inventory differs")
        for row in reader:
            if len(row) != len(header):
                raise ValueError("Truncated/extra native columns")
            members, occupied = [], 0
            for sp, cell in zip(species, row[len(metadata):]):
                if cell == "*" or not cell:
                    continue
                genes = cell.split(",")
                if any(not gene or owners.get(gene) != sp for gene in genes):
                    raise ValueError("Unknown gene or incorrect species column")
                occupied += 1
                members.extend(genes)
            count, species_count = (int(row[1]), int(row[2] if method == "sonicparanoid" else row[0]))
            if not members or len(members) != count or occupied != species_count:
                raise ValueError("Native gene/species counts differ")
            if len(members) != len(set(members)) or seen.intersection(members):
                raise ValueError("Duplicate native membership")
            if method == "sonicparanoid":
                if not row[0] or row[0] in labels:
                    raise ValueError("Duplicate/empty group ID")
                labels.add(row[0])
            seen.update(members)
            groups.append(tuple(members))
    if not groups:
        raise ValueError("Empty native table")
    padded = groups + [(gene,) for gene in sorted(owners.keys() - seen)]
    return canonical(padded, owners), len(groups), len(seen), len(owners) - len(seen)


def sonic_metadata(text):
    fields = {}
    for line in text.splitlines():
        if "\t" in line:
            key, value = line.split("\t", 1)
            if key in fields:
                raise ValueError("Duplicate run-info field")
            fields[key] = value
    if "SonicParanoid 2.0.9" not in text or fields.get("Input proteomes:") != "12":
        raise ValueError("Unexpected SonicParanoid version/input count")
    return fields


def audit(root):
    reports = root / "benchmark_tools/results"
    inventory_path, readback_path = (reports / name for name in
        ("installed_ob_input_inventory_20260926.json", "retained_ob_comparator_readback_20260926.json"))
    checked = [record(inventory_path), record(readback_path)]
    if [r["sha256"] for r in checked] != ["8b429398e6c381fcdec8289c641a7f20969edcad5595e3d3f6930335601571bd", "5dcb4e65f277eb9debf2bbacb1d4e298c678ac13fc981e75c4ba78e4851eeb55"]:
        raise ValueError("Changed frozen input/prediction evidence")
    inventory, readback = json.loads(inventory_path.read_text()), json.loads(readback_path.read_text())
    rows = []
    for method, key in (("sonicparanoid", "sonicparanoid_2_0_9"), ("proteinortho", "proteinortho_6_3_6")):
        base = root / "benchmarks/results" / (method + "_run")
        staged, owners = [], {}
        for source in inventory["inputs"]:
            check(source)
            path = base / "input" / Path(source["path"]).name
            item = record(path)
            checked.extend([source, item])
            staged.append(dict(source=source, staged=item,
                               bytes_match=all(source[k] == item[k] for k in ("sha256", "bytes"))))
            genes = fasta(path)
            if owners.keys() & genes.keys():
                raise ValueError("Duplicate input gene IDs")
            owners.update({gene: path.name for gene in genes})
        if len(owners) != inventory["genes"]:
            raise ValueError("Input count mismatch")
        log = base.with_name(method + "_run.log")
        if method == "sonicparanoid":
            run = base / "output/runs/orthobench"
            table, info = run / "ortholog_groups/ortholog_groups.tsv", run / "run_info.txt"
        else:
            table, info = base / "orthobench.proteinortho.tsv", base / "orthobench.info"
        checked.extend(record(p) for p in (table, info, log))
        evidence = next(r for r in readback["rows"] if r["key"] == key)["prediction"]
        check(evidence)
        checked.append(evidence)
        padded, count, assigned, unassigned = native_groups(table, method, owners)
        normalized = canonical(read_groups(Path(evidence["path"]), False), owners)
        text, log_text = info.read_text(), log.read_text()
        if method == "sonicparanoid":
            settings = sonic_metadata(text)
            elapsed = [float(line.split("\t", 1)[1]) for line in log_text.splitlines() if line.startswith("Total elapsed time (seconds):\t")]
            if len(elapsed) != 1 or elapsed[0] <= 0:
                raise ValueError("Missing/ambiguous native elapsed time")
            metadata = dict(settings=settings, command_argv=None, reported_version="2.0.9",
                            tool_reported_elapsed_seconds=elapsed[0], cpu_seconds=None, peak_rss_kib=None)
        else:
            calls = [shlex.split(line) for line in text.splitlines() if line.startswith("proteinortho ")]
            if not calls or "version=6.3.6," not in text or "All finished." not in log_text:
                raise ValueError("Missing Proteinortho command/version/completion evidence")
            metadata = dict(system_calls=calls, reported_version="6.3.6", recorded_invocations=len(calls),
                            tool_reported_elapsed_seconds=None, cpu_seconds=None, peak_rss_kib=None,
                            invocation_note="Multiple retained calls do not establish a single uninterrupted run.")
        rows.append(dict(key=key, staged_inputs=staged, input_hashes_match=all(r["bytes_match"] for r in staged),
                         native_table=record(table), normalized_prediction=evidence,
                         native_groups=count, native_assigned_genes=assigned, singleton_padding=unassigned,
                         conversion_partition_matches=padded == normalized, native_metadata=metadata,
                         historical_consumption_proven=False))
    checked.extend(record(Path(__file__).with_name(name)) for name in
                   ("audit_ob_matrix_provenance.py", "audit_ob_orthofinder_provenance.py"))
    for item in checked:
        check(item)
    return dict(status="retained_ob_matrix_provenance_readback", rows=rows, checked_records=checked,
                all_checks_passed=all(r["input_hashes_match"] and r["conversion_partition_matches"] for r in rows),
                publication_ready=False, inference_rerun=False,
                limitations=["Staged hashes and native-table conversion consistency do not independently attest historical consumption or binary identity.",
                             "Tool-reported SonicParanoid elapsed time is not externally measured end-to-end time; no CPU/RSS record is claimed.",
                             "Proteinortho info retains multiple invocations, with no independently identified per-invocation timing or outcome.",
                             "Singleton padding preserves pair co-membership but is not native assigned-gene coverage.",
                             "No controlled timing, native rerun or changed benchmark score."])


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
