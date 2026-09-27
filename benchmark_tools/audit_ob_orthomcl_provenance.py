"""Read retained OrthoMCL input, MCL conversion, command and timestamp evidence."""

import argparse
from datetime import datetime, timedelta, timezone
import json
from pathlib import Path
import shlex
import sys

from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.audit_ob_orthofinder_provenance import fasta
from benchmark_tools.compare_ob_orthomcl_runs import partition
from benchmark_tools.normalize_three_kingdoms_orthogroups import iter_orthomcl
from benchmark_tools.orthofinder_mcl_to_orthogroups import iter_mcl_clusters
from benchmark_tools.explain_ob_input_differences import differences

OUTPUTS = {"Apr_20":"6234d4fd517826fdf4fd9c25cca14b8d0bf58aef089c407eae7a57324e53e32a",
           "Jul_25":"1add0db6640e0451184975b5cbdf085ff404bbf24dadb66402921cd4dfd8a227"}


def parse_timestamp(value):
    parsed = datetime.strptime(value, "%a %b %d %I:%M:%S %p EDT %Y")
    if parsed.strftime("%a") != value.split()[0]:
        raise ValueError("Timestamp weekday disagrees with date")
    return parsed.replace(tzinfo=timezone(timedelta(hours=-4)))


def duration(start, end):
    a, b = parse_timestamp(start), parse_timestamp(end)
    seconds = (b-a).total_seconds()
    if seconds <= 0:
        raise ValueError("Nonpositive logged duration")
    return dict(start=a.isoformat(), end=b.isoformat(), seconds=seconds,
                timezone="Explicit EDT (UTC-04:00) in native log")


def after_marker(lines, marker):
    positions = [i for i, line in enumerate(lines) if line.strip() == marker]
    if len(positions) != 1:
        raise ValueError(f"Missing or repeated marker: {marker}")
    return next(line.strip() for line in lines[positions[0]+1:] if line.strip())


def flag(tokens, name):
    if tokens.count(name) != 1 or tokens.index(name)+1 == len(tokens):
        raise ValueError(f"Missing/repeated command option: {name}")
    return tokens[tokens.index(name)+1]


def mapping(path, species_genes):
    expected = {frozenset(genes):name for name,genes in species_genes.items()}
    if len(expected) != len(species_genes):
        raise ValueError("Ambiguous expected species membership")
    result, seen = {}, set()
    for line in path.read_text().splitlines():
        name, sep, remainder = line.partition(":")
        genes = remainder.split()
        group = frozenset(genes)
        if (not sep or not name or name in result or not genes or len(genes) != len(group)
                or group not in expected or seen & group):
            raise ValueError("Invalid species-gene mapping")
        result[name] = dict(input_species=expected[group],genes=len(group))
        seen.update(group)
    if len(result) != len(species_genes):
        raise ValueError("Incomplete species mapping")
    return result


def read_index(path, universe):
    index, seen = {}, set()
    for line in path.read_text().splitlines():
        fields = line.split()
        if (len(fields) != 2 or not fields[0].isdigit() or fields[0] in index
                or fields[1] not in universe or fields[1] in seen):
            raise ValueError("Invalid MCL sequence index")
        index[fields[0]] = fields[1]
        seen.add(fields[1])
    if set(map(int,index)) != set(range(len(index))) or not index:
        raise ValueError("Noncontiguous or empty MCL index")
    return index


def commands(directory):
    native = (directory / "orthomcl.log").read_text().splitlines()
    parameters = (directory / "parameter.log").read_text().splitlines()
    result = {key:shlex.split(after_marker(native,marker)) for key,marker in (
        ("formatdb","Native FormatDB command:"),("blastall","Native BLASTALL command:"),
        ("mcl","Run MCL program"))}
    result["orthomcl"] = shlex.split(after_marker(parameters,"########################COMMAND######################"))
    for tokens,key,value in ((result["formatdb"],"-i",str(directory / "tmp/all.fa")),
        (result["blastall"],"-i",str(directory / "tmp/all.fa")),
        (result["blastall"],"-d",str(directory / "tmp/all.fa")),
        (result["blastall"],"-o",str(directory / "tmp/all.blast")),
        (result["mcl"],"-o",str(directory / "tmp/all_ortho.mcl"))):
        if flag(tokens,key) != value:
            raise ValueError("Native command points to another run")
    start = [line.split(":",1)[1].strip() for line in native if line.startswith("Start Time:")]
    end = [line.split(":",1)[1].strip() for line in native if line.startswith("End Time:")]
    if len(start) != 1 or len(end) != 1:
        raise ValueError("Ambiguous run timestamps")
    if (start[0] != after_marker(parameters,"######################START TIME#####################")
            or end[0] != after_marker(parameters,"########################END TIME#####################")):
        raise ValueError("Native and parameter log timestamps disagree")
    return dict(commands=result, logged_duration=duration(start[0],end[0]),
                blast_threads=int(flag(result["blastall"],"-a")),
                blast_evalue=float(flag(result["blastall"],"-e")),
                mcl_inflation=float(flag(result["mcl"],"-I")))


def audit(repo, base):
    inventory_path = repo / "benchmark_tools/results/installed_ob_input_inventory_20260926.json"
    if record(inventory_path)["sha256"] != "8b429398e6c381fcdec8289c641a7f20969edcad5595e3d3f6930335601571bd":
        raise ValueError("Changed input inventory")
    inventory = json.loads(inventory_path.read_text())
    reference_path = repo / "benchmark_tools/results/orthobench_paired_uncertainty_20260916.json"
    if record(reference_path)["sha256"] != "660ead29c5b6ac0b8278cd1e62cdcdb0a513db81dda317e634b805e661d70ba9":
        raise ValueError("Changed reference inventory")
    reference_records = json.loads(reference_path.read_text())["inputs"]["references"]
    records = [record(inventory_path),*inventory["inputs"],record(reference_path),*reference_records]
    for item in records:
        check(item)
    refs = {Path(r["path"]).name:set(Path(r["path"]).read_text().splitlines()) for r in reference_records}
    if len(refs) != 70:
        raise ValueError("Wrong reference-family inventory")
    expected, by_species = {}, {}
    for item in inventory["inputs"]:
        path = Path(item["path"])
        sequences = fasta(path)
        if expected.keys() & sequences.keys():
            raise ValueError("Duplicate input gene")
        expected.update(sequences)
        by_species[path.name] = set(sequences)
    if len(expected) != 251378 or len(by_species) != 12:
        raise ValueError("Wrong OrthoBench universe")
    rows = []
    for label, checksum in OUTPUTS.items():
        directory = base / label
        paths = [directory / p for p in ("all_orthomcl.out","orthomcl.log","parameter.log",
                 "tmp/all.fa","tmp/all.gg","tmp/all_ortho.idx","tmp/all_ortho.mcl")]
        inputs = [record(p) for p in paths]
        if inputs[0]["sha256"] != checksum:
            raise ValueError("Changed native partition")
        records.extend(inputs)
        aggregate = fasta(directory / "tmp/all.fa")
        delta = differences(expected,aggregate,refs)
        if delta["only_original"] or delta["only_staged"]:
            raise ValueError("Native aggregate gene universe differs from frozen inputs")
        species = mapping(directory / "tmp/all.gg",by_species)
        index = read_index(directory / "tmp/all_ortho.idx",set(expected))
        translated = ((str(i),tuple(index[g] for g in genes)) for i,genes in
                      enumerate(iter_mcl_clusters(directory / "tmp/all_ortho.mcl")))
        mcl, mcl_genes = partition(translated)
        native, native_genes = partition(iter_orthomcl(directory / "all_orthomcl.out"))
        if mcl != native or mcl_genes != native_genes or mcl_genes != set(index.values()):
            raise ValueError("MCL conversion differs from native partition")
        details = commands(directory)
        filenames = flag(details["commands"]["orthomcl"],"--fa_files").split(",")
        if set(filenames) != set(by_species) or len(filenames) != len(by_species):
            raise ValueError("Command input inventory differs")
        rows.append(dict(run=label,**details,species_mapping=species,input_genes=len(aggregate),
            sequences_exact=aggregate == expected,input_differences=delta,
            mcl_conversion_exact=True,groups=len(native),assigned_genes=len(native_genes),
            unassigned_genes=len(expected)-len(native_genes),input_records=inputs))
    for item in records:
        check(item)
    return dict(status="retained_orthomcl_provenance_readback",runs=rows,
        input_identity_all_runs=all(r["sequences_exact"] for r in rows),
        checked_records=records,source=record(__file__),scores_changed=False,
        limitations=["Retained log timestamps are descriptive wall intervals, not controlled timing or scheduler accounting.",
            "BLAST thread requests are not measured CPU use; peak RSS and CPU time remain unknown.",
            "Executable version-looking paths are not historical binary attestations.",
            "No complete BLAST/BPO/matrix or transitive execution audit; no inference rerun.",
            "Generic MCL parser is shared with another audit; this does not independently validate MCL inference.",
            "Sequence differences remain a failed identity gate; transformation descriptions do not identify an actor or causal accuracy effect.",
            "Reference exposure uses full family membership before low-certainty exclusions, not an effect bound.",
            "April 2026 and July 2025 are different retained runs; no metadata or timing transfer between them."])


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo",type=Path,required=True)
    p.add_argument("--base",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    result = audit(a.repo.resolve(),a.base.resolve())
    with a.output.open("x") as stream:
        json.dump(result,stream,indent=2,sort_keys=True)
        stream.write("\n")
    if not result["input_identity_all_runs"]:
        sys.exit(1)
