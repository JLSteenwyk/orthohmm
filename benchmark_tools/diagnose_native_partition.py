"""Read-only coverage diagnosis; never benchmark scoring or terminal admission."""

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path
import sys

from Bio import SeqIO


class Bindings:
    def __init__(self):
        self.files = {}

    def bind(self, path, expected=None):
        path = Path(path)
        if not path.is_absolute() or path.resolve() != path or not path.is_file():
            raise ValueError("Require a direct absolute file")
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        ref = dict(path=str(path), bytes=path.stat().st_size, sha256=digest.hexdigest())
        if expected is not None and ref != expected:
            raise ValueError("Pinned evidence differs: " + str(path))
        if str(path) in self.files and self.files[str(path)] != ref:
            raise ValueError("Evidence changed during diagnosis: " + str(path))
        self.files[str(path)] = ref
        return path

    def read(self, path, expected=None):
        return json.loads(self.bind(path, expected).read_text())

    def finish(self):
        for ref in list(self.files.values()):
            self.bind(ref["path"], ref)
        return list(self.files.values())


def root_groups(path):
    groups, sources = {}, defaultdict(list)
    with Path(path).open() as handle:
        rows = csv.reader(handle, delimiter="\t")
        if next(rows, None) != ["root_hog", "source_family", "genes"]:
            raise ValueError("Invalid root-HOG header")
        for row in rows:
            if len(row) != 3 or not row[0] or not row[1] or row[0] in groups:
                raise ValueError("Invalid or duplicate root-HOG row")
            genes = row[2].split(",")
            if any(not g or g.strip() != g for g in genes):
                raise ValueError("Invalid root-HOG gene")
            groups[row[0]] = genes
            sources[row[1]].extend(genes)
    return groups, dict(sources)


def text_groups(path, named=False):
    groups = {}
    with Path(path).open() as handle:
        for line in handle:
            if not line.strip():
                continue
            if named:
                name, separator, text = line.partition(":")
                if not separator or not name.strip() or name.strip() in groups:
                    raise ValueError("Invalid named group")
                name = name.strip()
            else:
                name, text = str(len(groups)), line
            genes = text.split()
            if not genes:
                raise ValueError("Empty group")
            groups[name] = genes
    return groups


def coverage(groups, owners):
    counts = Counter(g for genes in groups.values() for g in genes)
    missing, extra = sorted(owners.keys() - counts.keys()), sorted(counts.keys() - owners.keys())
    duplicates = {g: n for g, n in sorted(counts.items()) if n != 1}
    taxa = sorted(set(owners.values()))
    per_species = {taxon: dict(input=0, observed_unique=0, memberships=0, missing=0)
                   for taxon in taxa}
    for gene, taxon in owners.items():
        per_species[taxon]["input"] += 1
        per_species[taxon]["observed_unique"] += int(gene in counts)
        per_species[taxon]["memberships"] += counts[gene]
        per_species[taxon]["missing"] += int(gene not in counts)
    return dict(groups=len(groups), unique_genes=len(counts), memberships=sum(counts.values()),
                missing=missing, extra=extra, duplicates=duplicates, per_species=per_species,
                complete_unique_input_partition=not (missing or extra or duplicates))


def canonical_partition_digest(groups):
    # Names are not part of co-membership; preserve duplicates rather than deduplicate.
    digest = hashlib.sha256()
    for genes in sorted(tuple(sorted(genes)) for genes in groups.values()):
        digest.update(json.dumps(genes, separators=(",", ":")).encode() + b"\n")
    return digest.hexdigest()


def source_payload_digest(groups):
    digest = hashlib.sha256()
    for name in sorted(groups):
        digest.update((" ".join(sorted(groups[name])) + "\n").encode())
    return digest.hexdigest()


def diagnose(plan_ref, index, failure_ref):
    from benchmark_tools import score_ygob_groups as frozen_parser
    from benchmark_tools import validate_native_factorial_outputs as frozen_validator

    evidence = Bindings()
    plan = evidence.read(plan_ref["path"], plan_ref)
    failure = evidence.read(failure_ref["path"], failure_ref)
    if type(index) is not int or index < 0 or index >= len(plan["runs"]):
        raise ValueError("Invalid run index")
    run = plan["runs"][index]
    if run["index"] != index or failure["index"] != index or failure["plan"] != plan_ref:
        raise ValueError("Failure does not bind this run")
    if (failure.get("status") != "terminal_factorial_review_failed"
            or failure.get("terminal_reviewed") is not False
            or failure.get("next_identity_authorized") is not False
            or failure.get("accuracy_evaluated") is not False
            or failure.get("automatic_retry") is not False):
        raise ValueError("Require the retained failed review, not a new disposition")
    root = Path(run["output_root"])
    owners, per_file, digest = {}, {}, hashlib.sha256()
    inputs = {Path(ref["path"]).name: ref for ref in run["inputs"]}
    order = run["native_order"]
    if (len(inputs) != len(run["inputs"]) or len(order) != len(set(order))
            or set(order) != set(inputs) or len(inputs) != run["proteomes"]
            or {p.name for p in (root / "input").iterdir()} != set(inputs)):
        raise ValueError("Input inventory differs")
    for name in order:
        evidence.bind(inputs[name]["path"], inputs[name])
        path = evidence.bind(root / "input" / name)
        if any(evidence.files[str(path)][k] != inputs[name][k] for k in ("bytes", "sha256")):
            raise ValueError("Prepared input bytes differ")
        per_file[name] = 0
        for sequence in SeqIO.parse(path, "fasta"):
            if not sequence.id or sequence.id in owners or not sequence.seq:
                raise ValueError("Invalid or duplicate input record")
            owners[sequence.id] = Path(name).stem
            per_file[name] += 1
            digest.update(json.dumps([name, sequence.id], separators=(",", ":")).encode() + b"\n")
        if not per_file[name]:
            raise ValueError("Empty input proteome")
    preparation = evidence.read(root / "preparation.json")
    if (len(owners) != run["genes"] or preparation["genes"] != len(owners)
            or preparation["gene_ownership_sha256"] != digest.hexdigest()
            or preparation["per_species_counts"] != per_file):
        raise ValueError("Prepared gene ownership differs")
    working = root / "native/orthohmm_working_res"
    names_path = evidence.bind(working / "high_sensitivity_checkpoint/gene_names.txt")
    names = names_path.read_text().splitlines()
    checkpoint = evidence.read(working / "high_sensitivity_checkpoint/manifest.json")
    if checkpoint["files"]["gene_names.txt"] != {
            k: evidence.files[str(names_path)][k] for k in ("bytes", "sha256")}:
        raise ValueError("Checkpoint name checksum differs")
    if checkpoint["genes"] != len(names):
        raise ValueError("Checkpoint name count differs")
    root_path = evidence.bind(root / "native/orthohmm_phylogeny/orthohmm_root_hogs.tsv")
    independent, sources = root_groups(root_path)
    clustered = text_groups(evidence.bind(working / "orthohmm_edges_clustered.txt"))
    materialized = text_groups(evidence.bind(root / "native/orthohmm_orthogroups.txt"), named=True)
    frozen = frozen_parser.read_predictions(root_path, "root_hogs")
    stages = dict(checkpoint_names=coverage({"checkpoint": names}, owners),
                  root_hogs=coverage(independent, owners), clustered=coverage(clustered, owners),
                  materialized=coverage(materialized, owners), frozen_parser=coverage(frozen, owners))
    frozen_gate = {}
    try:
        value = frozen_validator.groups(root_path, "root_hogs", owners, frozen_validator.Evidence())
        frozen_gate.update(status="passed", groups=len(value))
    except Exception as error:
        frozen_gate.update(status="failed", error_type=type(error).__name__, error=str(error))
    provenance = evidence.read(root / "native/orthohmm_phylogeny/provenance_manifest.json")
    summary = evidence.read(root / "native/orthohmm_phylogeny/reconciliation_summary.json")
    partitions = {name: canonical_partition_digest(value) for name, value in
                  dict(root_hogs=independent, frozen_parser=frozen, clustered=clustered,
                       materialized=materialized).items()}
    for module in (frozen_parser, frozen_validator):
        evidence.bind(Path(module.__file__).resolve())
    evidence.bind(Path(__file__).resolve())
    return dict(schema="native_partition_diagnosis_v1", status="diagnosis_completed",
                index=index, job_id=failure["job_id"], plan=plan_ref, original_failure=failure_ref,
                input_genes=len(owners), input_proteomes=len(inputs), ownership_sha256=digest.hexdigest(),
                stages=stages, frozen_root_coverage_gate=frozen_gate,
                parsers_identical=frozen == independent, checkpoint_names_lexical=names == sorted(owners),
                partition_digests=partitions, all_output_partitions_identical=len(set(partitions.values())) == 1,
                reconstructed_source_families=len(sources),
                source_family_ids_complete=sorted(sources) == [
                    f"Family{i:07d}" for i in range(summary["candidate_families"])],
                source_payload_sha256=source_payload_digest(sources),
                source_payload_matches_native_input=source_payload_digest(sources) == provenance["input_cluster_sha256"],
                root_count_matches_summary=len(independent) == summary["root_hogs"],
                evidence=evidence.finish(), python=sys.version,
                accuracy_evaluated=False, terminal_reviewed=False, next_identity_authorized=False,
                native_outputs_validated=False, automatic_retry=False, resources_admitted=False,
                historical_failure_cause_established=False,
                limitations=["Current coverage does not establish why the historical review failed.",
                             "No full semantic/resource review, scoring, failure reclassification or retry.",
                             "Source families are reconstructed from native labels, not independent biological truth."])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--plan-sha256", required=True)
    parser.add_argument("--index", type=int, required=True)
    parser.add_argument("--failure", type=Path, required=True)
    parser.add_argument("--failure-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink() or not args.output.parent.is_dir():
        raise ValueError("Require a fresh output in an existing directory")
    refs = Bindings()
    refs.bind(args.plan)
    refs.bind(args.failure)
    plan_ref, failure_ref = refs.files[str(args.plan)], refs.files[str(args.failure)]
    if plan_ref["sha256"] != args.plan_sha256 or failure_ref["sha256"] != args.failure_sha256:
        raise ValueError("Selected plan/failure hash differs")
    result = diagnose(plan_ref, args.index, failure_ref)
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps({k: result[k] for k in ("status", "index", "input_genes", "input_proteomes",
                    "frozen_root_coverage_gate", "all_output_partitions_identical",
                    "source_payload_matches_native_input", "historical_failure_cause_established")}, sort_keys=True))


if __name__ == "__main__":
    main()
