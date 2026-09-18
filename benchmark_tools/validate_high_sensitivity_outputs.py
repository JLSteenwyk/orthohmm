"""Cross-check native high-sensitivity groups, metrics, checkpoint and FASTAs."""

import json
from pathlib import Path

from Bio import SeqIO
from orthohmm.accuracy import load_accuracy_checkpoint
from benchmark_tools.audit_accuracy_checkpoint import audit
from benchmark_tools.audit_qfo_replay_inputs import verify_species_partition
from benchmark_tools.score_ygob_groups import read_predictions, membership
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check


def validate_metrics(metrics, genes, species, groups):
    if metrics.get("status") != "complete":
        raise ValueError("Native metrics are not complete")
    metadata = metrics.get("metadata")
    if not isinstance(metadata, dict):
        raise ValueError("Missing native metadata")
    expected = {"accuracy_profile": "high_sensitivity", "search_mode": "builtin", "clustering": "leiden",
                "cpm_resolution": 0.1, "substitution_matrix": "BLOSUM62", "evalue_threshold": 1e-4,
                "leiden_seed": 4, "cpu_budget": 32}
    for key, value in expected.items():
        if metadata.get(key) != value:
            raise ValueError("Unexpected native high-sensitivity parameter: " + key)
    counts = metrics.get("counts")
    if not isinstance(counts, dict):
        raise ValueError("Missing native counts")
    for key, value in (("genes", genes), ("species", species), ("orthogroups", groups)):
        if type(counts.get(key)) is not int or counts[key] != value:
            raise ValueError("Native output count mismatch: " + key)
    profiles = counts.get("high_sensitivity_profiles")
    if type(profiles) is not int or profiles <= 0:
        raise ValueError("Missing positive HMM profile-build evidence")


def validate_partition(raw_path, groups, names):
    with raw_path.open() as stream:
        raw = {str(i): line.split() for i, line in enumerate(stream) if line.strip()}
    index = membership(raw)
    universe = set(names)
    if not set(index) <= universe:
        raise ValueError("Raw clustering contains foreign genes")
    missing = universe - set(index)
    # Native export preserves clusters and adds each unclustered gene alone.
    expected = {frozenset(genes) for genes in raw.values()}
    expected.update(frozenset([gene]) for gene in missing)
    if {frozenset(genes) for genes in groups.values()} != expected:
        raise ValueError("Native groups differ from raw clusters plus singletons")
    return {"raw_clusters": len(raw), "added_singletons": len(missing)}


def validate(output, metrics_path, inputs, checkpoint_sha):
    output, metrics_path = Path(output), Path(metrics_path)
    inputs = [Path(p) for p in inputs]
    checkpoint = output / "orthohmm_working_res/high_sensitivity_checkpoint"
    group_path = output / "orthohmm_orthogroups.txt"
    raw_path = output / "orthohmm_working_res/orthohmm_edges_clustered.txt"
    checkpoint_paths = sorted(checkpoint.iterdir())
    files = [metrics_path, group_path, raw_path, *inputs, *checkpoint_paths]
    records = [record(p) for p in files]
    numeric = audit(checkpoint, checkpoint_sha)
    names, species, *_ = load_accuracy_checkpoint(checkpoint, verify=False)
    ownership = verify_species_partition(names, species,
        ((str(p), (entry.id for entry in SeqIO.parse(p, "fasta"))) for p in inputs))
    groups = read_predictions(group_path, "named_groups")
    index = membership(groups)
    if set(index) != set(names):
        raise ValueError("Native groups differ from complete checkpoint/FASTA universe")
    partition = validate_partition(raw_path, groups, names)
    metrics = json.loads(metrics_path.read_text())
    validate_metrics(metrics, len(names), len(ownership), len(groups))
    if Path(metrics["metadata"]["output_directory"]).resolve() != output.resolve():
        raise ValueError("Metrics output directory differs")
    if {p.parent.resolve() for p in inputs} != {Path(metrics["metadata"]["fasta_directory"]).resolve()}:
        raise ValueError("Metrics FASTA directory differs")
    for item in records:
        check(item)
    if sorted(checkpoint.iterdir()) != checkpoint_paths:
        raise ValueError("Checkpoint inventory changed during validation")
    return {"status": "high_sensitivity_output_content_verified", "source": record(__file__),
            "checked_records": records, "numeric_checkpoint": numeric, "species_ownership": ownership,
            "partition": partition,
            "genes": len(names), "groups": len(groups), "singleton_groups": sum(len(g) == 1 for g in groups.values()),
            "metrics": record(metrics_path), "native_groups": record(group_path), "accuracy_evaluated": False,
            "limitations": ["Content consistency only; terminal scheduler, source/runtime and command provenance require separate admission.",
                            "Checkpoint hash binds this content; it does not prove search completeness or independent score correctness.",
                            "Profile-build evidence is native instrumentation, not independent rescoring of HMMs."]}
