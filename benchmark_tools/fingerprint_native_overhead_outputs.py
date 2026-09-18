"""Canonical output identities for paired native overhead validation."""

import hashlib
import json
from pathlib import Path

from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.report_ygob_validation import read_checkpoint
from benchmark_tools.score_ygob_groups import membership, read_predictions
from benchmark_tools.simulation_method_outputs import unique_path, orthohmm_pairs, orthofinder_pairs
from benchmark_tools.validate_scaling_outputs import input_universe, relocate_evidence


def digest_rows(rows):
    digest = hashlib.sha256()
    for row in rows:
        digest.update(json.dumps(row, ensure_ascii=True, separators=(",", ":")).encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


def partition_identity(groups, owners):
    if set(membership(groups)) != set(owners):
        raise ValueError("Partition must cover exactly the frozen input universe")
    canonical = sorted(tuple(sorted(genes)) for genes in groups.values())
    return dict(groups=len(canonical), genes=len(owners), sha256=digest_rows(canonical))


def pair_identity(pairs, owners):
    canonical = set()
    for pair in pairs:
        if not isinstance(pair, (list, tuple)) or len(pair) != 2:
            raise ValueError("Require two gene identifiers per native pair")
        a, b = pair
        if (not isinstance(a, str) or not isinstance(b, str) or a not in owners
                or b not in owners or owners[a] == owners[b]):
            raise ValueError("Native pairs must contain known cross-species genes")
        canonical.add(tuple(sorted((a, b))))
    return dict(pairs=len(canonical), sha256=digest_rows(sorted(canonical)))


def fingerprint(run, evidence_roots):
    accessible, _ = relocate_evidence(run, evidence_roots)
    owners, species = input_universe(accessible["dataset"])
    output = Path(accessible["configuration"]["output"])
    method = run["native_method"]
    paths = [Path(row["path"]) for row in accessible["dataset"]["inputs"]]
    identity = dict(input_species_mapping=digest_rows(sorted(owners.items())))
    if method in {"orthohmm_high_sensitivity", "orthohmm_satellite_v2"}:
        groups = output / "orthohmm_orthogroups.txt"
        paths.append(groups)
        roots = output / "orthohmm_phylogeny/orthohmm_root_hogs.tsv"
        pairs = output / "orthohmm_phylogeny/orthohmm_pairwise_orthologs.tsv"
        if method == "orthohmm_satellite_v2":
            paths.extend([roots, pairs])
    elif method == "orthofinder_full":
        mapping = unique_path(output, "**/SequenceIDs.txt")
        clusters = unique_path(output, "**/clusters_OrthoFinder_I*.txt_id_pairs.txt")
        relations = unique_path(output, "**/Orthologues")
        paths.extend([mapping, clusters, *sorted(relations.glob("Orthologues_*/*__v__*.tsv"))])
    else:
        raise ValueError("Unknown native overhead method")
    evidence = [record(path) for path in paths]
    if method.startswith("orthohmm_"):
        identity["orthogroups"] = partition_identity(read_predictions(groups, "named_groups"), owners)
        if method == "orthohmm_satellite_v2":
            identity["root_hogs"] = partition_identity(read_predictions(roots, "root_hogs"), owners)
            identity["native_pairs"] = pair_identity(orthohmm_pairs(pairs, owners), owners)
    else:
        identity["checkpoint_groups"] = partition_identity(read_checkpoint(clusters, mapping, owners), owners)
        identity["native_pairs"] = pair_identity(orthofinder_pairs(relations.parent, owners, species), owners)
    for row in evidence:
        check(row)
    return dict(identity=identity, evidence=evidence, source=record(__file__),
        helpers=[record(Path(__file__).with_name(name)) for name in (
            "score_ygob_groups.py", "report_ygob_validation.py", "simulation_method_outputs.py",
            "validate_scaling_outputs.py", "orthofinder_mcl_to_orthogroups.py",
            "orthofinder_to_pairwise.py", "prepare_ob_candidate_neighborhood.py")],
        scientific_timings_admitted=False,
        limitations=["Canonical output equality, not identical internal computational work or biological correctness.",
            "Group labels, member/row order and native pair orientation do not affect identities.",
            "Pair-set fingerprinting deduplicates; native-format duplicate validation remains a separate gate.",
            "Caller must independently validate command completion, source/runtime, native formats and scheduler evidence."])
