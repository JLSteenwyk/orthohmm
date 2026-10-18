"""Read back original local profiles/trees and build an authenticated cache.

These helpers repeat shared validated algorithms, not independent algorithm
implementations. They never replace a missing/tampered source tree or write
official orthology predictions.
"""

import ast
import json
from pathlib import Path
import shutil

from Bio import SeqIO

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from orthohmm.local_rescue import (
    _constraint_evidence, _context_constraints, assess_joint_outcome, validate_local_profiles,
)
from orthohmm.phylogeny import parse_species_tree
from orthohmm.phylogeny_pipeline import (
    _aligner_command, _family_tree_input_hash, _gene_tokens, _load_reusable_raw_trees,
    _reconcile_family, _tree_command, apply_satellite_membership_constraints,
)


def scientific_compatibility(original, reader):
    """Prove that used original local scientific code has not changed."""
    original, reader = Path(original).resolve(), Path(reader).resolve()
    paths = [p.relative_to(original) for p in sorted((original / "orthohmm").rglob("*.py"))]
    paths += [Path("benchmark_tools") / name for name in ("run_local_rescue_probe.py", "replay_phylogeny.py")]
    checked, differences = [], []
    for relative in paths:
        old, new = original / relative, reader / relative
        before, after = old.read_bytes(), new.read_bytes()
        if relative == Path("orthohmm/search/profile_expansion.py"):
            trees = [ast.parse(payload) for payload in (before, after)]
            for tree in trees:
                tree.body = [node for node in tree.body if not isinstance(node, ast.FunctionDef) or node.name != "expand_profiles"]
            equal = ast.dump(trees[0], include_attributes=False) == ast.dump(trees[1], include_attributes=False)
        else:
            equal = before == after
        if not equal:
            differences.append(str(relative))
        checked.extend((record(old), record(new)))
    extras = {p.relative_to(reader) for p in (reader / "orthohmm").rglob("*.py")} - set(paths)
    if differences or extras - {Path("orthohmm/search/selective_profiles.py")}:
        raise ValueError("Local readback scientific code differs from original probe: " + str(differences or sorted(map(str, extras))))
    return {"status": "original_local_scientific_code_unchanged", "checked": checked,
            "exception": "Only the unused expand_profiles hook and new selective-only module are excluded; original local profile building/scoring and all native reconciliation/writers are unchanged."}


def authenticate_tree(directory, context, sequences, species_tree_pin, config):
    directory = Path(directory)
    family = context["context_id"]
    genes = tuple(sorted(context["genes_a"] + context["genes_b"]))
    checkpoint_path = directory / "checkpoints" / (family + ".json")
    checkpoint = json.loads(checkpoint_path.read_text())
    paths = {key: directory / "gene_trees" / (family + suffix) for key, suffix in (
        ("raw_tree_sha256", ".raw.nwk"), ("rooted_tree_sha256", ".rooted.nwk"),
        ("annotated_tree_sha256", ".reconciled.nwk"))}
    if (checkpoint.get("status") != "complete" or checkpoint.get("schema_version") != 2
            or checkpoint.get("family_id") != family or checkpoint.get("genes") != list(genes)
            or checkpoint.get("tree_input_sha256") != _family_tree_input_hash(family, genes, sequences, config)
            or checkpoint.get("species_tree_sha256") != species_tree_pin["sha256"]
            or checkpoint.get("gene_tree_support_mode") != "legacy"
            or any(record(path)["sha256"] != checkpoint.get(key) for key, path in paths.items())):
        raise ValueError("Local tree checkpoint/input/raw/rooted/reconciled identity differs")
    mapping = _gene_tokens(genes)
    for label in ("candidate_fastas", "alignments"):
        path = directory / label / (family + ".faa")
        rows = list(SeqIO.parse(path, "fasta"))
        if len(rows) != len(mapping) or {r.id for r in rows} != set(mapping):
            raise ValueError("Local tree FASTA/alignment token ownership differs")
        for row in rows:
            expected = sequences[mapping[row.id]]
            normalize = lambda s: str(s).upper().replace("-", "").replace(".", "")
            if (str(row.seq) != expected if label == "candidate_fastas" else normalize(row.seq) != normalize(expected)):
                raise ValueError("Local tree FASTA/alignment changed source residues")
    cache = _load_reusable_raw_trees(directory)
    if set(cache) != {genes} or cache[genes].tree_input_sha256 != checkpoint["tree_input_sha256"]:
        raise ValueError("Original local raw tree is not exactly reusable")
    pins = [record(checkpoint_path), *[record(p) for p in paths.values()],
            *[record(directory / label / (family + ".faa")) for label in ("candidate_fastas", "alignments")]]
    return cache, pins


def read_context(proposal_path, context, gene_names, database, gene_species, sequences,
                 species_tree_pin, proteome_files, config, constraints, baseline_high_pairs,
                 output_directory, *, cpu=16):
    directory = Path(output_directory)
    directory.mkdir(parents=True, exist_ok=False)
    proposal_path = Path(proposal_path)
    proposal_pin = record(proposal_path)
    proposal = json.loads(proposal_path.read_text())
    checked = [proposal_pin, species_tree_pin]
    check(species_tree_pin)
    groups = [context[key] for key in ("genes_a", "genes_b")]
    selected_constraints = _context_constraints(constraints, groups)
    expected_tree = {"path": species_tree_pin["path"], "sha256": species_tree_pin["sha256"]}
    if (proposal["context"] != context or proposal["species_tree"] != expected_tree
            or proposal["constraints"] != _constraint_evidence(selected_constraints)
            or proposal["families_merged"] is not False or proposal["final_native_inference"] is not False
            or "error" in proposal):
        raise ValueError("Original proposal differs from context/species-tree/constraints or is unfinished")
    profiles = validate_local_profiles(context, gene_names, database, gene_species, cpu=cpu)
    if profiles != proposal["profile_validation"]:
        raise ValueError("Original profile construction/scores/spans/decision do not reproduce")
    result = {"context_id": context["context_id"], "proposal": proposal_pin,
              "profile_accepted": profiles["accepted"], "tree_inferred": False, "accepted": False,
              "status": "profile_rejection_read_back", "checked": checked}
    if not profiles["accepted"]:
        if (proposal["status"] != "rejected_by_profiles" or proposal["accepted"] is not False
                or proposal["tree_inferred"] is not False or "joint_admission" in proposal):
            raise ValueError("Profile rejection unexpectedly inferred/admitted a joint tree")
    else:
        if proposal["status"] != "joint_tree_evaluated" or proposal["tree_inferred"] is not True:
            raise ValueError("Require completed source joint-tree evaluation")
        source_directory = proposal_path.parent
        cache, tree_pins = authenticate_tree(source_directory, context, sequences, species_tree_pin, config)
        checked.extend(tree_pins)
        family = context["context_id"]
        commands = {"alignment": _aligner_command(config.aligner, source_directory / "candidate_fastas" / (family + ".faa")),
                    "gene_tree": _tree_command(config.tree_builder, source_directory / "alignments" / (family + ".faa"))}
        if proposal["commands"] != commands:
            raise ValueError("Original local tree commands differ from fixed policies")
        validated = parse_species_tree(species_tree_pin["path"], proteome_files)
        genes = tuple(sorted(context["genes_a"] + context["genes_b"]))
        outcome = _reconcile_family(family, genes, sequences, gene_species, validated.tree,
            species_tree_pin["sha256"], config, directory, cache)
        if not outcome.checkpoint_hit or outcome.remapped_checkpoint_hit:
            raise ValueError("Local readback failed to reuse the exact source raw tree")
        constrained, audit = apply_satellite_membership_constraints([outcome], selected_constraints)
        decision = assess_joint_outcome(constrained[0], *groups, baseline_high_pairs, gene_species)
        if (decision != proposal["joint_admission"] or audit != proposal["constraint_audit"]
                or proposal["accepted"] != decision["accepted"] or proposal["token_to_gene"] != _gene_tokens(genes)):
            raise ValueError("Original joint reconciliation/membership/admission does not reproduce")
        for suffix in ("raw", "rooted", "reconciled"):
            filename = family + "." + suffix + ".nwk"
            if record(directory / "gene_trees" / filename)["sha256"] != record(source_directory / "gene_trees" / filename)["sha256"]:
                raise ValueError("Read-back raw/rooted/reconciled tree differs from original")
        result.update(status="joint_tree_admission_read_back", tree_inferred=True, accepted=decision["accepted"],
            joint_admission=decision, baseline_high_pairs=[list(pair) for pair in baseline_high_pairs],
            original_directory=str(source_directory), genes=list(genes))
    for pin in checked:
        check(pin)
    return result


def build_composite_cache(baseline_phylogeny, contexts, output_directory):
    """Copy baseline cache plus only authenticated accepted raw-tree contexts."""
    baseline, destination = Path(baseline_phylogeny).resolve(), Path(output_directory).resolve()
    if destination.exists() or destination.resolve() == baseline.resolve():
        raise FileExistsError("Require a fresh composite cache namespace")
    accepted = [row for row in contexts if row["accepted"]]
    source_pins = [record(p) for p in sorted(baseline.rglob("*")) if p.is_file()]
    for row in accepted:
        if row["status"] != "joint_tree_admission_read_back" or not row["tree_inferred"] or not row["profile_accepted"]:
            raise ValueError("Composite cache requires independently read-back admissions")
        for pin in row["checked"]:
            check(pin)
    shutil.copytree(baseline, destination, copy_function=shutil.copy2)
    for name in ("gene_trees", "checkpoints"):
        (destination / name).mkdir(exist_ok=True)
    copied = []
    for pin in source_pins:
        source = Path(pin["path"])
        target = destination / source.relative_to(baseline)
        actual = record(target)
        if any(actual[key] != pin[key] for key in ("bytes", "sha256")):
            raise ValueError("Composite baseline copy differs")
        copied.append({"source": pin, "copy": actual})
    for row in accepted:
        family, source = row["context_id"], Path(row["original_directory"])
        allowed = {pin["path"] for pin in row["checked"]}
        for relative in (Path("checkpoints") / (family + ".json"), Path("gene_trees") / (family + ".raw.nwk")):
            pin = record(source / relative)
            if pin["path"] not in allowed or (destination / relative).exists():
                raise ValueError("Unauthenticated local cache file or context-ID collision")
            shutil.copy2(source / relative, destination / relative)
            actual = record(destination / relative)
            if any(actual[key] != pin[key] for key in ("bytes", "sha256")):
                raise ValueError("Composite joint-tree copy differs")
            copied.append({"source": pin, "copy": actual})
    cache = _load_reusable_raw_trees(destination)
    for row in accepted:
        key = tuple(sorted(row["genes"]))
        if key not in cache or cache[key].family_id != row["context_id"]:
            raise ValueError("Accepted joint raw tree is missing/ambiguous in composite cache")
    for row in copied:
        check(row["source"])
        check(row["copy"])
    return {"status": "authenticated_composite_cache_built", "path": str(destination.resolve()),
            "accepted_contexts": [row["context_id"] for row in accepted], "files": copied,
            "scope": "Baseline cache plus accepted raw-tree/checkpoint pairs only; no provisional species-tree cache replacement or official predictions."}
