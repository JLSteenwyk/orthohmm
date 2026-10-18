"""Label-independent, bounded family-HMM and local joint-tree proposals.

Nomination and cross-profile support do not establish orthology. Admitted
joint contexts still require ordinary final native inference and evaluation.
"""

from dataclasses import dataclass
import json
import math
from numbers import Integral
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class _DirectionalSupport:
    score: float
    queries: tuple[int, ...]
    targets: tuple[int, ...]


def _positive_integer(value, name, minimum=1):
    if (isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral)
            or value < minimum):
        raise ValueError(name + " must be a bounded positive integer")
    return int(value)


def _two_witnesses(old, first, second, sentinel):
    return tuple(sorted({*old, *[int(value) for value in (first, second) if value < sentinel]})[:2])


def _partition(clusters, names, species, max_joint_genes):
    count = len(names)
    if count > np.iinfo(np.int32).max or any(
        not isinstance(name, str) or not name for name in names
    ) or len(set(names)) != count:
        raise ValueError("Require unique nonempty global gene names within int32 bounds")
    if species.ndim != 1 or len(species) != count or species.dtype.kind not in "iu" or np.any(species < 0):
        raise ValueError("Require one nonnegative integer species label per gene")
    owner = np.full(count, -1, dtype=np.int32)
    eligible = np.zeros(len(clusters), dtype=bool)
    family_names, excluded = {}, []
    for family, cluster in enumerate(clusters):
        members = np.asarray(cluster)
        if (members.ndim != 1 or not len(members) or members.dtype.kind not in "iu"
                or np.any(members < 0) or np.any(members >= count)
                or len(np.unique(members)) != len(members)):
            raise ValueError("Candidate families must contain unique valid gene indices")
        if np.any(owner[members] >= 0):
            raise ValueError("Candidate families overlap")
        owner[members] = family
        species_count = len(np.unique(species[members]))
        reasons = []
        if len(members) < 3:
            reasons.append("too_few_genes")
        if len(members) > max_joint_genes - 3:
            reasons.append("family_exceeds_joint_budget")
        if species_count < 2:
            reasons.append("too_few_species")
        if reasons:
            excluded.append({"family_id": family, "genes": len(members), "species": species_count, "reasons": reasons})
        else:
            eligible[family] = True
            family_names[family] = tuple(sorted(names[int(index)] for index in members))
    if np.any(owner < 0):
        raise ValueError("Candidate partition does not cover the complete gene table")
    return owner, eligible, family_names, excluded


def _aggregate_chunk(stats, codes, queries, targets, scores, sentinel):
    if not len(codes):
        return
    order = np.argsort(codes, kind="stable")
    keys, q, t, values = codes[order], queries[order], targets[order], scores[order]
    starts = np.r_[0, np.flatnonzero(keys[1:] != keys[:-1]) + 1]
    sizes = np.diff(np.r_[starts, len(keys)])
    maxima = np.maximum.reduceat(values, starts)
    # Two distinct witnesses suffice for the gate; repeated hits add no support.
    first_q, first_t = np.minimum.reduceat(q, starts), np.minimum.reduceat(t, starts)
    second_q = np.minimum.reduceat(np.where(q != np.repeat(first_q, sizes), q, sentinel), starts)
    second_t = np.minimum.reduceat(np.where(t != np.repeat(first_t, sizes), t, sentinel), starts)
    for index, start in enumerate(starts):
        key = int(keys[start])
        previous = stats.get(key, _DirectionalSupport(0., (), ()))
        stats[key] = _DirectionalSupport(
            max(previous.score, float(maxima[index])),
            _two_witnesses(previous.queries, first_q[index], second_q[index], sentinel),
            _two_witnesses(previous.targets, first_t[index], second_t[index], sentinel),
        )


def nominate_local_contexts(
    clusters, gene_names, gene_to_species, hit_queries, hit_targets, hit_scores,
    *, max_contexts=128, max_joint_genes=200, chunk_size=500000,
):
    """Nominate disjoint unique mutual-best family pairs from existing HMM hits.

    Directional witnesses saturate at two distinct genes; they are not exact
    support counts or necessarily reciprocal individual gene-pair matches.
    No benchmark identity or reference annotations are accepted by this API.
    """
    max_contexts = _positive_integer(max_contexts, "Context budget")
    max_joint_genes = _positive_integer(max_joint_genes, "Joint gene budget", minimum=6)
    chunk_size = _positive_integer(chunk_size, "Hit chunk size")
    names = list(gene_names)
    species = np.asarray(gene_to_species)
    if not len(names) and species.shape == (0,):
        species = species.astype(np.int32)
    owner, eligible, family_names, excluded = _partition(clusters, names, species, max_joint_genes)
    queries, targets, scores = (np.asarray(array) for array in (hit_queries, hit_targets, hit_scores))
    if any(array.ndim != 1 for array in (queries, targets, scores)) or not len(queries) == len(targets) == len(scores):
        raise ValueError("HMM hit arrays must be aligned one-dimensional arrays")
    if not len(queries):
        queries, targets = queries.astype(np.int32), targets.astype(np.int32)
    if queries.dtype.kind not in "iu" or targets.dtype.kind not in "iu" or scores.dtype.kind not in "iuf":
        raise ValueError("HMM hit indices must be integers and scores real numbers")
    stats, used_hits = {}, 0
    families = len(clusters)
    for start in range(0, len(scores), chunk_size):
        q, t, s = queries[start:start + chunk_size], targets[start:start + chunk_size], scores[start:start + chunk_size]
        if (np.any(q < 0) or np.any(t < 0) or np.any(q >= len(names)) or np.any(t >= len(names))
                or not np.isfinite(s).all()):
            raise ValueError("Invalid HMM hit index or nonfinite score")
        left, right = owner[q], owner[t]
        keep = (left != right) & (species[q] != species[t]) & eligible[left] & eligible[right] & (s > 0)
        q, t, s, left, right = q[keep], t[keep], s[keep], left[keep], right[keep]
        used_hits += len(s)
        codes = left.astype(np.int64) * families + right
        _aggregate_chunk(stats, codes, q, t, s, len(names))

    unordered = sorted({tuple(sorted(divmod(code, families))) for code in stats})
    supported, pair_audit, neighbors = {}, [], {}
    for left, right in unordered:
        forward, reverse = stats.get(left * families + right), stats.get(right * families + left)
        if (forward is None or reverse is None or any(
            len(witnesses) < 2 for evidence in (forward, reverse) for witnesses in (evidence.queries, evidence.targets)
        )):
            pair_audit.append({"family_a": left, "family_b": right, "status": "insufficient_reciprocal_witnesses"})
            continue
        score = math.sqrt(forward.score) * math.sqrt(reverse.score)
        if not math.isfinite(score):
            raise ValueError("Nonfinite reciprocal family score")
        supported[(left, right)] = score
        neighbors.setdefault(left, []).append((score, right))
        neighbors.setdefault(right, []).append((score, left))
    ranking, best = {}, {}
    for family in sorted(family_names):
        ordered = sorted(neighbors.get(family, []), key=lambda item: (-item[0], item[1]))
        tied = [neighbor for score, neighbor in ordered if score == ordered[0][0]] if ordered else []
        selected = tied[0] if len(tied) == 1 else None
        best[family] = selected
        ranking[family] = {"family_id": family, "best_neighbor": selected,
            "best_score": ordered[0][0] if ordered else None, "best_tied_neighbors": sorted(tied),
            "runner_up_neighbor": ordered[1][1] if len(ordered) > 1 else None,
            "runner_up_score": ordered[1][0] if len(ordered) > 1 else None,
            "score_margin": ordered[0][0] - ordered[1][0] if len(ordered) > 1 else None}
    mutual = []
    for (left, right), score in sorted(supported.items()):
        row = {"family_a": left, "family_b": right, "score": score,
               "joint_genes": len(clusters[left]) + len(clusters[right])}
        if best[left] != right or best[right] != left:
            row["status"] = "not_unique_mutual_best"
        elif row["joint_genes"] > max_joint_genes:
            row["status"] = "joint_gene_budget"
        else:
            row["status"] = "eligible_before_context_budget"
            mutual.append(row)
        pair_audit.append(row)
    mutual.sort(key=lambda row: (-row["score"], *sorted((family_names[row["family_a"]], family_names[row["family_b"]]))))
    contexts = []
    for index, row in enumerate(mutual):
        if index >= max_contexts:
            row["status"] = "context_budget"
            continue
        row["status"] = "selected"
        left, right = row["family_a"], row["family_b"]
        contexts.append({**row, "context_id": f"LocalContext{index + 1:06d}",
            "genes_a": list(family_names[left]), "genes_b": list(family_names[right]),
            "a_to_b_score": stats[left * families + right].score,
            "b_to_a_score": stats[right * families + left].score,
            "family_a_ranking": ranking[left], "family_b_ranking": ranking[right]})
    selected_families = [family for row in contexts for family in (row["family_a"], row["family_b"])]
    if len(selected_families) != len(set(selected_families)):
        raise ValueError("Mutual-best nomination unexpectedly produced overlapping contexts")
    return {
        "schema_version": 1, "status": "nominated_not_profile_or_tree_validated",
        "parameters": {"max_contexts": max_contexts, "max_joint_genes": max_joint_genes,
                       "minimum_family_genes": 3, "minimum_family_species": 2,
                       "minimum_directional_witnesses": 2, "hit_chunk_size": chunk_size},
        "counts": {"genes": len(names), "families": families, "eligible_families": int(eligible.sum()),
                   "input_hits": len(scores), "eligible_positive_cross_hits": used_hits,
                   "directed_family_pairs": len(stats), "reciprocally_supported_pairs": len(supported),
                   "selected_contexts": len(contexts)},
        "excluded_families": excluded, "family_ranking": [ranking[key] for key in sorted(ranking)],
        "pair_audit": sorted(pair_audit, key=lambda row: (row["family_a"], row["family_b"])),
        "directional_support": [{"family_query": code // families, "family_target": code % families,
            "max_score": evidence.score, "query_witnesses": list(evidence.queries),
            "target_witnesses": list(evidence.targets)} for code, evidence in sorted(stats.items())],
        "contexts": contexts, "orthology_inferred": False, "families_merged": False,
    }


def assess_joint_outcome(outcome, genes_a, genes_b, baseline_high_pairs, gene_species):
    """Admit only a constrained reconciled context preserving baseline pairs.

    This gate does not write candidate groups or append native pairs. The
    final ordinary species-tree/reconciliation pass can still change its
    provisional relationships and must be evaluated separately.
    """
    genes_a, genes_b = tuple(genes_a), tuple(genes_b)
    if (not genes_a or not genes_b or any(not isinstance(gene, str) or not gene for gene in (*genes_a, *genes_b))
            or len(set(genes_a)) != len(genes_a) or len(set(genes_b)) != len(genes_b)
            or set(genes_a) & set(genes_b)):
        raise ValueError("Joint admission requires disjoint nonempty unique source families")
    left, right = set(genes_a), set(genes_b)
    universe = left | right
    if (len(outcome.genes) != len(universe) or set(outcome.genes) != universe
            or outcome.reconciliation is None):
        raise ValueError("Joint admission requires the complete reconciled two-family context")
    if any(gene not in gene_species or not isinstance(gene_species[gene], str) or not gene_species[gene]
           for gene in universe):
        raise ValueError("Joint context species ownership is missing or invalid")
    roots = {}
    for index, group in enumerate(outcome.root_groups):
        if not group:
            raise ValueError("Joint root groups must be nonempty")
        for gene in group:
            if gene not in universe or gene in roots:
                raise ValueError("Joint root groups have invalid or duplicate ownership")
            roots[gene] = index
    if set(roots) != universe:
        raise ValueError("Joint root groups do not cover the complete context")

    def pair_key(pair):
        if (len(pair) != 2 or any(not isinstance(gene, str) or gene not in universe for gene in pair) or pair[0] == pair[1]
                or gene_species[pair[0]] == gene_species[pair[1]]):
            raise ValueError("Invalid joint cross-species ortholog pair")
        return tuple(sorted(pair))

    pairs = {pair_key(pair) for pair in outcome.ortholog_pairs}
    if len(pairs) != len(outcome.ortholog_pairs) or any(roots[a] != roots[b] for a, b in pairs):
        raise ValueError("Joint pairs are duplicated or cross root-group boundaries")
    confidence = {}
    for row in outcome.ortholog_pair_confidence:
        if len(row) != 3 or row[2] not in ("high", "medium", "low"):
            raise ValueError("Invalid joint pair confidence")
        pair = pair_key(row[:2])
        if pair in confidence:
            raise ValueError("Duplicate joint pair confidence")
        confidence[pair] = row[2]
    if set(confidence) != pairs:
        raise ValueError("Joint pair/confidence universes differ")
    baseline_rows = tuple(baseline_high_pairs)
    baseline = {pair_key(pair) for pair in baseline_rows}
    if len(baseline) != len(baseline_rows):
        raise ValueError("Baseline high pairs are duplicated")
    if any((a in left) != (b in left) for a, b in baseline):
        raise ValueError("Baseline high pairs cross the original family boundaries")
    lost = baseline - pairs
    downgraded = {pair for pair in baseline & pairs if confidence[pair] != "high"}
    cross = {pair for pair in pairs if confidence[pair] == "high" and ((pair[0] in left) != (pair[1] in left))}
    supported_a = {gene for pair in cross for gene in pair if gene in left}
    supported_b = {gene for pair in cross for gene in pair if gene in right}
    species_a = {gene_species[gene] for gene in supported_a}
    species_b = {gene_species[gene] for gene in supported_b}
    reasons = []
    if lost:
        reasons.append("baseline_high_pair_loss")
    if min(len(supported_a), len(supported_b)) < 2:
        reasons.append("insufficient_novel_high_gene_support")
    if min(len(species_a), len(species_b)) < 2:
        reasons.append("insufficient_novel_high_species_support")
    reconciliation = outcome.reconciliation
    return {
        "accepted": not reasons, "reasons": reasons, "family_id": outcome.family_id,
        "genes": len(universe), "root_groups": len(outcome.root_groups),
        "baseline_high_pairs": len(baseline), "baseline_high_pairs_retained": len(baseline & pairs),
        "baseline_high_pairs_lost": [list(pair) for pair in sorted(lost)],
        "baseline_high_pairs_downgraded": [list(pair) for pair in sorted(downgraded)],
        "novel_cross_high_pairs": [list(pair) for pair in sorted(cross)],
        "supported_genes_a": sorted(supported_a), "supported_genes_b": sorted(supported_b),
        "supported_species_a": sorted(species_a), "supported_species_b": sorted(species_b),
        "duplications": reconciliation.duplication_count, "speciations": reconciliation.speciation_count,
        "uncertain_events": reconciliation.uncertain_count,
        "scope": "Provisional constrained joint-tree admission, not a final pair/group export or accuracy claim",
    }


def _local_families(context, gene_species):
    groups = [tuple(context[key]) for key in ("genes_a", "genes_b")]
    if any(len(group) < 3 or any(not isinstance(gene, str) or not gene for gene in group)
           or len(set(group)) != len(group) for group in groups):
        raise ValueError("Local profiles require complete unique families of at least three genes")
    if set(groups[0]) & set(groups[1]) or sum(map(len, groups)) > 200:
        raise ValueError("Local families must be disjoint and within the joint 200-gene budget")
    for group in groups:
        if any(gene not in gene_species or not isinstance(gene_species[gene], str)
               or not gene_species[gene] for gene in group):
            raise ValueError("Missing local gene species ownership")
        if len({gene_species[gene] for gene in group}) < 2:
            raise ValueError("Local families require at least two species")
    return [tuple(sorted(group)) for group in groups]


def validate_local_profiles(context, gene_names, sequence_database, gene_species, *, cpu=1):
    """Validate a nominated context with reciprocal, unweighted CPU profiles.

    The database must be the complete proteome collection, not a neighborhood
    subset. Its total residue count is used for the existing approximate
    significance calculation. This routine does not infer orthology.
    """
    from .search.evalue import batch_estimate_evalues
    from .search.matrices import get_ka_params
    from .search.profile_expansion import build_cluster_profiles, _pack_profiles
    from .search.profile_spans import profile_alignment_span
    from .search.viterbi import is_c_available, batch_viterbi_multipair_c, batch_viterbi_c

    cpu = _positive_integer(cpu, "CPU count")
    groups = _local_families(context, gene_species)
    names = list(gene_names)
    if (not names or any(not isinstance(gene, str) or not gene for gene in names)
            or len(names) != len(set(names)) or names != list(sequence_database.ids)
            or len(names) > np.iinfo(np.int32).max):
        raise ValueError("Local profiles require the unique complete database-aligned gene table")
    lookup = {gene: index for index, gene in enumerate(names)}
    if any(gene not in lookup for group in groups for gene in group):
        raise ValueError("Local family gene is absent from the complete sequence database")
    lengths, offsets, flat = (np.asarray(array) for array in (
        sequence_database.lengths, sequence_database.offsets, sequence_database.flat_sequences))
    if (lengths.shape != (len(names),) or offsets.shape != lengths.shape or flat.ndim != 1
            or lengths.dtype != np.dtype("int32") or offsets.dtype != np.dtype("int64")
            or flat.dtype != np.dtype("uint8") or np.any(lengths < 0)
            or not np.array_equal(offsets, np.r_[0, np.cumsum(lengths[:-1], dtype=np.int64)])
            or int(lengths.sum(dtype=np.int64)) != len(flat)):
        raise ValueError("Invalid complete sequence database layout for CPU scoring")
    if not is_c_available():
        raise RuntimeError("Local profile protocol requires the original CPU Viterbi library")
    clusters = [[lookup[gene] for gene in group] for group in groups]
    profiles = build_cluster_profiles(clusters, names, sequence_database, "BLOSUM62", cpu)
    result = {
        "schema_version": 1, "genes_a": list(groups[0]), "genes_b": list(groups[1]),
        "accepted": False, "reasons": [], "orthology_inferred": False,
        "parameters": {"matrix": "BLOSUM62", "sequence_weighting": "none", "evalue": 1e-4,
                       "band_width": 64, "coverage_rule": "5*known_matches >= 4*min(model_length,target_length)",
                       "database_genes": len(names), "database_residues": len(flat)},
        "profiles": [], "hits": [], "supported_genes_a": [], "supported_genes_b": [],
        "supported_species_a": [], "supported_species_b": [],
    }
    for index, profile in sorted(profiles.items()):
        if profile.sequence_weighting != "none" or list(profile.member_ids) != list(groups[index]):
            raise ValueError("Local profile sources or weighting differ from the fixed protocol")
        result["profiles"].append({"side": "ab"[index], "members": list(profile.member_ids),
            "length": int(profile.length), "sequence_weighting": profile.sequence_weighting,
            "effective_sequence_count": profile.effective_sequence_count,
            "match_emissions": profile.match_emissions.tolist(), "insert_emissions": profile.insert_emissions.tolist(),
            "transitions": profile.transitions.tolist(), "consensus": profile.consensus.tolist()})
    if set(profiles) != {0, 1}:
        result["reasons"].append("missing_informative_profile")
        return result
    _, packed = _pack_profiles(profiles)
    profile_lengths, profile_offsets, emissions, _, insertions, transitions = packed
    pairs = np.asarray([(index, lookup[gene]) for index in (0, 1) for gene in groups[1 - index]], dtype=np.int32)
    args = (emissions, insertions, transitions, profile_offsets, profile_lengths,
            flat, offsets, lengths, pairs, 64)
    try:
        scores = batch_viterbi_multipair_c(*args, n_threads=cpu)
        result["scorer"] = "original_cpu_multipair"
    except RuntimeError:
        scores = batch_viterbi_c(*args, n_threads=cpu)
        result["scorer"] = "original_cpu_scalar_fallback"
    if (scores.shape != (len(pairs),) or scores.dtype.kind not in "iu" or np.any(scores < 0)):
        raise ValueError("Invalid raw CPU profile scores")
    lam, k = get_ka_params("BLOSUM62")
    evalues = batch_estimate_evalues(scores, profile_lengths[pairs[:, 0]], len(flat), lam, k)
    if not np.isfinite(evalues).all() or np.any(evalues < 0):
        raise ValueError("Invalid approximate profile significance")
    supported = [set(), set()]
    for (profile_id, gene_id), score, evalue in zip(pairs, scores, evalues):
        profile_id, gene_id = int(profile_id), int(gene_id)
        row = {"profile_side": "ab"[profile_id], "target_side": "ab"[1 - profile_id],
               "target_gene": names[gene_id], "raw_score": int(score), "evalue": float(evalue),
               "significant": bool(evalue < 1e-4), "supported": False, "span": None}
        if row["significant"]:
            span = profile_alignment_span(profiles[profile_id], sequence_database.get_sequence(gene_id),
                                          band_width=64, expected_score=int(score))
            row["span"] = span.as_record()
            row["supported"] = bool(span.passes_guard)
            if span.passes_guard:
                supported[1 - profile_id].add(names[gene_id])
        row["reason"] = ("supported" if row["supported"] else
                         "aligned_coverage" if row["significant"] else "profile_significance")
        result["hits"].append(row)
    for index, side in enumerate("ab"):
        species = sorted({gene_species[gene] for gene in supported[index]})
        result["supported_genes_" + side] = sorted(supported[index])
        result["supported_species_" + side] = species
        if len(supported[index]) < 2:
            result["reasons"].append("insufficient_profile_genes_" + side)
        if len(species) < 2:
            result["reasons"].append("insufficient_profile_species_" + side)
    result["accepted"] = not result["reasons"]
    return result


def _context_constraints(constraints, groups):
    sides = [set(group) for group in groups]
    universe, selected = set.union(*sides), []
    for row in constraints:
        source, target = tuple(row["source_genes"]), tuple(row["target_genes"])
        if (not source or not target or any(not isinstance(gene, str) or not gene for gene in (*source, *target))
                or len(set(source)) != len(source) or len(set(target)) != len(target) or set(source) & set(target)):
            raise ValueError("Invalid original family membership constraint")
        genes = set(source) | set(target)
        if not genes & universe:
            continue
        if not any(genes <= side for side in sides):
            raise ValueError("Original membership constraint crosses local source family boundaries")
        selected.append(row)
    return selected


def _constraint_evidence(value):
    """Encode legacy diagnostic nonfinites without modifying inference inputs."""
    if isinstance(value, float) and not math.isfinite(value):
        label = "nan" if math.isnan(value) else "positive_infinity" if value > 0 else "negative_infinity"
        return {"nonfinite_constraint_diagnostic": label}
    if isinstance(value, dict):
        return {key: _constraint_evidence(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_constraint_evidence(item) for item in value]
    return value


def evaluate_local_context(context, gene_names, sequence_database, gene_species, sequences,
                           species_tree_path, proteome_files, config, output_directory,
                           constraints, baseline_high_pairs, *, cpu=1):
    """Validate profiles and infer one fresh, constrained provisional joint tree.

    The caller authenticates that species_tree_path belongs to the unmodified
    baseline. No baseline checkpoint is reused for this fresh joint context.
    Final candidate/native inference and benchmark evaluation remain separate.
    """
    from .phylogeny import parse_species_tree
    from .phylogeny_pipeline import (_atomic_write, _sha256_file, _reconcile_family,
                                     _gene_tokens, _aligner_command, _tree_command,
                                     apply_satellite_membership_constraints)

    groups = _local_families(context, gene_species)
    expected = {"mode": "reconcile", "species_tree_mode": "infer", "species_tree_rooting": "min_variance",
                "root_duplication_rule": "species_overlap", "pair_orthology_rule": "positive_paralogy",
                "satellite_membership_rule": "high_confidence_pair", "gene_tree_support_mode": "legacy"}
    if any(getattr(config, key) != value for key, value in expected.items()):
        raise ValueError("Local rescue requires the fixed baseline phylogeny policies")
    family_id = context["context_id"]
    if (not isinstance(family_id, str) or not family_id.startswith("LocalContext")
            or not family_id[len("LocalContext"):].isdigit()):
        raise ValueError("Invalid local context identifier")
    selected_constraints = _context_constraints(constraints, groups)
    universe = set(groups[0]) | set(groups[1])
    original_pairs = []
    for pair in baseline_high_pairs:
        if len(pair) != 2 or pair[0] == pair[1] or any(gene not in gene_species for gene in pair):
            raise ValueError("Invalid baseline high-confidence pair")
        if set(pair) & universe:
            if not any(set(pair) <= set(group) for group in groups):
                raise ValueError("Baseline pair crosses original candidate family boundaries")
            original_pairs.append(pair)
    from .search.sequences import encode_sequence
    indices = {gene: index for index, gene in enumerate(gene_names)}
    if any(gene not in sequences or not isinstance(sequences[gene], str) or gene not in indices
           or not np.array_equal(encode_sequence(sequences[gene]), sequence_database.get_sequence(indices[gene]))
           for gene in universe):
        raise ValueError("Joint tree and profile sequences differ")
    validated = parse_species_tree(species_tree_path, proteome_files)
    if not {gene_species[gene] for gene in universe} <= set(validated.taxa):
        raise ValueError("Local context species are absent from the baseline species tree")
    directory = Path(output_directory)
    directory.mkdir(parents=True, exist_ok=False)
    result = {"schema_version": 1, "context": context, "profile_validation": None,
              "accepted": False, "status": "profile_validation_running", "tree_inferred": False,
              "species_tree": {"path": str(validated.path.resolve()), "sha256": _sha256_file(validated.path)},
              "constraints": _constraint_evidence(selected_constraints),
              "constraint_evidence_encoding": "Only legacy constraint diagnostic nonfinites are tagged. Inference receives the original unmodified constraints; profile/tree scores remain strict finite JSON.",
              "families_merged": False, "final_native_inference": False}
    receipt = directory / "proposal.json"
    _atomic_write(receipt, json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
    try:
        evidence = validate_local_profiles(context, gene_names, sequence_database, gene_species, cpu=cpu)
    except Exception as error:
        result.update({"status": "profile_validation_failed", "error": {
            "type": type(error).__name__, "message": str(error)}})
        _atomic_write(receipt, json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
        raise
    result.update({"profile_validation": evidence, "status": "rejected_by_profiles"})
    _atomic_write(receipt, json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
    if not evidence["accepted"]:
        return result
    result["status"] = "joint_tree_running"
    result["commands"] = {
        "alignment": _aligner_command(config.aligner, directory / "candidate_fastas" / f"{family_id}.faa"),
        "gene_tree": _tree_command(config.tree_builder, directory / "alignments" / f"{family_id}.faa")}
    _atomic_write(receipt, json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
    genes = tuple(sorted(universe))
    try:
        outcome = _reconcile_family(family_id, genes, sequences, gene_species, validated.tree,
                                    result["species_tree"]["sha256"], config, directory, {})
        result["tree_inferred"] = True
        if outcome.checkpoint_hit or outcome.remapped_checkpoint_hit:
            raise ValueError("Provisional local context unexpectedly reused a gene tree")
        constrained, audit = apply_satellite_membership_constraints([outcome], selected_constraints)
        decision = assess_joint_outcome(constrained[0], *groups, original_pairs, gene_species)
    except Exception as error:
        result.update({"status": "joint_tree_failed", "error": {
            "type": type(error).__name__, "message": str(error)}})
        _atomic_write(receipt, json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
        raise
    result.update({"status": "joint_tree_evaluated", "tree_inferred": True, "accepted": decision["accepted"],
                   "joint_admission": decision, "constraint_audit": audit, "token_to_gene": _gene_tokens(genes)})
    _atomic_write(receipt, json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
    return result


def write_effective_candidates(original_path, gene_names, proposals, output_path, constraints):
    """Join only admitted, disjoint whole-family contexts for final inference.

    This writes candidate groups, never official orthogroups or native pairs.
    With no admitted contexts the original file is preserved byte for byte.
    The caller must authenticate the proposal evidence before using it.
    """
    from .phylogeny_pipeline import _atomic_write

    original_path, output_path = Path(original_path), Path(output_path)
    if output_path.exists() or output_path.is_symlink() or output_path.resolve() == original_path.resolve():
        raise FileExistsError(output_path)
    raw = original_path.read_bytes()
    groups = [tuple(line.split()) for line in raw.decode("utf-8").splitlines() if line.split()]
    names, owner = list(gene_names), {}
    for index, group in enumerate(groups):
        for gene in group:
            if gene in owner:
                raise ValueError("Original candidate partition has duplicate ownership")
            owner[gene] = index
    if (any(not isinstance(gene, str) or not gene for gene in names) or len(set(names)) != len(names)
            or set(owner) != set(names)):
        raise ValueError("Original candidate partition does not match the complete gene table")
    consumed, joins = set(), {}
    accepted = []
    for proposal in proposals:
        if not proposal.get("accepted", False):
            continue
        profile, joint = proposal.get("profile_validation"), proposal.get("joint_admission")
        if (proposal.get("accepted") is not True or proposal.get("status") != "joint_tree_evaluated"
                or proposal.get("tree_inferred") is not True or not isinstance(profile, dict)
                or profile.get("accepted") is not True or not isinstance(joint, dict)
                or joint.get("accepted") is not True):
            raise ValueError("Accepted proposal lacks completed profile and joint-tree admission")
        context = proposal["context"]
        sides = [tuple(context[key]) for key in ("genes_a", "genes_b")]
        indices = []
        for side in sides:
            if not side or len(set(side)) != len(side) or any(gene not in owner for gene in side):
                raise ValueError("Proposal has invalid source family genes")
            index = owner[side[0]]
            if set(side) != set(groups[index]):
                raise ValueError("Proposal does not preserve a complete original source family")
            indices.append(index)
        if indices[0] == indices[1] or any(index in consumed for index in indices):
            raise ValueError("Accepted local contexts overlap")
        consumed.update(indices)
        joins[min(indices)] = tuple(sorted(sides[0] + sides[1]))
        accepted.append({"context_id": context["context_id"], "original_family_indices": indices,
                         "genes": len(joins[min(indices)])})
    effective = [joins[index] if index in joins else group
                 for index, group in enumerate(groups) if index not in consumed or index in joins]
    effective_owner = {gene: index for index, group in enumerate(effective) for gene in group}
    for constraint in constraints:
        source, target = tuple(constraint["source_genes"]), tuple(constraint["target_genes"])
        genes = source + target
        if (not source or not target or len(set(genes)) != len(genes) or any(gene not in owner for gene in genes)
                or len({owner[gene] for gene in genes}) != 1
                or len({effective_owner[gene] for gene in genes}) != 1):
            raise ValueError("Membership constraint is not intact inside one original and effective family")
    if joins:
        _atomic_write(output_path, "".join(" ".join(group) + "\n" for group in effective))
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open("xb") as stream:
            stream.write(raw)
    return {"schema_version": 1, "original_families": len(groups), "effective_families": len(effective),
            "genes": len(owner), "accepted_contexts": accepted, "constraints_preserved": len(constraints),
            "original_bytes_preserved": not joins, "final_native_inference": False,
            "scope": "Candidate partition only; final native inference is required for all official outputs"}
