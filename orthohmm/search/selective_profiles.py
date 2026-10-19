"""Label-free nomination and same-rule decisions for selective profile pilots.

These experimental helpers do not call HMMER, change initial hits, or infer
orthology. Replacement operates on directed profile anchors before ordinary
undirected deduplication, so an unselected reciprocal anchor is not lost.
"""

import hashlib
import json
import math
from numbers import Integral
import os
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from ..accuracy import deduplicate_undirected_edges


MAX_QUERIES = 128
MAX_PROFILES = 128
MAX_EXTERNAL = 4
EVALUE = 0.0001


def _indices(values, size, label):
    array = np.asarray(values)
    if (array.ndim != 1 or (array.size and array.dtype.kind not in "iu")
            or np.any(array < 0) or np.any(array >= size)):
        raise ValueError("Invalid " + label + " indices")
    return array.astype(np.int32, copy=False)


def _numbers(values, length, label, *, nonnegative=False):
    array = np.asarray(values, dtype=np.float64)
    if (array.ndim != 1 or len(array) != length or not np.isfinite(array).all()
            or (nonnegative and np.any(array < 0))):
        raise ValueError("Invalid " + label)
    return array


def _partition(names, clusters, profiles, thresholds):
    if not names or len(set(names)) != len(names) or any(not isinstance(x, str) or not x for x in names):
        raise ValueError("Require unique nonempty gene names")
    owners = np.full(len(names), -1, dtype=np.int32)
    ranks = np.full(len(names), -1, dtype=np.int32)
    for family, group in enumerate(clusters):
        members = _indices(group, len(names), "source member")
        if not len(members) or len(np.unique(members)) != len(members) or np.any(owners[members] != -1):
            raise ValueError("Source families must uniquely partition all genes")
        owners[members], ranks[members] = family, np.arange(len(members))
    if np.any(owners < 0) or set(profiles) != set(thresholds):
        raise ValueError("Incomplete source ownership or self thresholds")
    for family, profile in profiles.items():
        if (type(family) is not int or not 0 <= family < len(clusters)
                or not 3 <= len(clusters[family]) <= 200
                or list(profile.member_ids) != [names[int(g)] for g in clusters[family]]
                or profile.sequence_weighting != "none" or not math.isfinite(float(thresholds[family]))):
            raise ValueError("Require original unweighted source profiles and finite thresholds")
    return owners, ranks


def nominate_queries(gene_names, clusters, profiles, profile_hits, self_thresholds,
                     hit_queries, hit_targets, hit_scores):
    """Select the complete bounded ambiguous pool without reference labels.

    All original significant external competitors are retained. Initial-hit
    arrays are streamed twice; per-member direction records are retained only
    for selected queries, not for the complete all-by-all search.
    """
    names = list(gene_names)
    owners, ranks = _partition(names, clusters, profiles, self_thresholds)
    profile_ids = _indices(profile_hits.profile_cluster_ids, len(clusters), "profile")
    genes = _indices(profile_hits.gene_ids, len(names), "profile query")
    scores = _numbers(profile_hits.scores, len(genes), "profile scores")
    evalues = _numbers(profile_hits.evalues, len(genes), "profile E-values", nonnegative=True)
    if len(profile_ids) != len(genes) or np.any(evalues >= EVALUE) or not set(map(int, profile_ids)) <= set(profiles):
        raise ValueError("Require complete original significant profile-hit rows")
    queries = _indices(hit_queries, len(names), "initial query")
    targets = _indices(hit_targets, len(names), "initial target")
    sequence_scores = _numbers(hit_scores, len(queries), "initial scores")
    if len(targets) != len(queries):
        raise ValueError("Initial hit arrays differ in length")
    external = {}
    for query, family, score in zip(genes, profile_ids, scores):
        query, family = int(query), int(family)
        if owners[query] == family:
            continue
        found = external.setdefault(query, {})
        found[family] = max(found.get(family, -math.inf), float(score))
    family_key = {family: tuple(sorted(profile.member_ids)) for family, profile in profiles.items()}
    audit, potential = {}, {}
    for query, hits in external.items():
        row = {"query_id": query, "query_gene": names[query], "external_profiles": sorted(hits, key=family_key.get)}
        audit[query] = row
        if len(hits) > MAX_EXTERNAL:
            row["status"] = "excluded_complete_competition_over_budget"
            continue
        best = min(hits, key=lambda family: (-hits[family], family_key[family]))
        if len(hits) == 1 and hits[best] >= self_thresholds[best]:
            row["status"] = "not_ambiguous"
            continue
        row["priority"] = abs(hits[best] - self_thresholds[best]) / max(1., abs(self_thresholds[best]))
        row["priority_profile_id"] = best
        potential[query] = row
    potential_mask = np.zeros(len(names), dtype=bool)
    potential_mask[list(potential)] = True

    def directed_hits(mask):
        for start in range(0, len(queries), 500000):
            q, t, s = (array[start:start + 500000] for array in (queries, targets, sequence_scores))
            keep = (q != t) & (s > 0)
            for reverse in (False, True):
                left, right = (t, q) if reverse else (q, t)
                rows = np.flatnonzero(keep & mask[left])
                for index in rows:
                    query, member = int(left[index]), int(right[index])
                    family = int(owners[member])
                    if family in external[query]:
                        yield query, family, member, float(s[index]), reverse

    anchored = {query for query, *_ in directed_hits(potential_mask)}
    for query in potential.keys() - anchored:
        audit[query]["status"] = "excluded_no_positive_sequence_anchor"
    selected, selected_profiles = [], set()
    for query in sorted(anchored, key=lambda q: (audit[q]["priority"], names[q])):
        row = audit[query]
        own = int(owners[query]) if int(owners[query]) in profiles else None
        competitors = set(external[query]) | ({own} if own is not None else set())
        row["own_profile_id"] = own
        if len(selected) >= MAX_QUERIES:
            row["status"] = "excluded_query_budget"
        elif len(selected_profiles | competitors) > MAX_PROFILES:
            row["status"] = "excluded_profile_budget"
        else:
            row["status"] = "selected"
            selected.append(query)
            selected_profiles.update(competitors)
    selected_mask = np.zeros(len(names), dtype=bool)
    selected_mask[selected] = True
    directions = {}
    for query, family, member, score, reverse in directed_hits(selected_mask):
        pair = directions.setdefault((query, family, member), [None, None])
        index = int(reverse)
        pair[index] = score if pair[index] is None else max(pair[index], score)
    anchors = {}
    for (query, family, member), (outgoing, reverse) in directions.items():
        # The historical rule ignores a reverse score if a positive outgoing
        # record for the same pair exists, even when reverse scores higher.
        score = outgoing if outgoing is not None else reverse
        previous = anchors.get((query, family))
        if previous is None or (-score, int(ranks[member])) < (-previous["score"], int(ranks[previous["gene_id"]])):
            anchors[query, family] = {"gene_id": member, "gene": names[member], "score": score}
    result = []
    for query in selected:
        own = audit[query]["own_profile_id"]
        competitors = set(external[query]) | ({own} if own is not None else set())
        result.append({**audit[query], "competitors": [{"profile_id": family,
            "members": list(profiles[family].member_ids), "is_own": family == own,
            "original_hit_score": external[query].get(family),
            "original_self_threshold": float(self_thresholds[family]),
            "anchor": anchors.get((query, family))} for family in sorted(competitors, key=family_key.get)]})
    return {"schema_version": 1, "selected": result, "candidate_digest": query_list_digest(result),
            "audit": [audit[q] for q in sorted(audit, key=lambda q: names[q])],
            "without_external_hits": sorted(set(names) - {names[q] for q in external}),
            "selected_profiles": sorted(selected_profiles, key=family_key.get),
            "orthology_inferred": False}


def query_list_digest(selected):
    """Identify complete nominated comparisons independently of integer IDs."""
    canonical = [{"query": row["query_gene"], "priority": row["priority"],
        "competitors": [{"members": sorted(c["members"]), "original_hit_score": c["original_hit_score"],
        "original_self_threshold": c["original_self_threshold"],
        "is_own": c["is_own"], "anchor": None if c["anchor"] is None else
        {"gene": c["anchor"]["gene"], "score": c["anchor"]["score"]}} for c in row["competitors"]]} for row in selected]
    return hashlib.sha256(json.dumps(canonical, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def decide_query(query, measurements):
    """Choose a unique complete-data winner before significance/anchor gates.

    Measurements use the common full_sequence_score/sequence_evalue schema.
    Missing competitor or calibration data preserve the original directed
    anchor. Malformed scientific values raise instead of becoming rejection.
    """
    competitors = query["competitors"]
    if (not 1 <= len(competitors) <= 5 or len({c["profile_id"] for c in competitors}) != len(competitors)
            or sum(bool(c["is_own"]) for c in competitors) > 1):
        raise ValueError("Invalid complete competitor set")
    values, missing = {}, []
    for competitor in competitors:
        family = competitor["profile_id"]
        measured = measurements.get(family)
        if measured is None:
            missing.append({"profile_id": family, "reason": "profile_missing"})
            continue
        if measured.get("status") not in ("hmmer_profile_scored", "builtin_profile_scored"):
            raise ValueError("Failed or unfinished backend is not a missing measurement")
        if (measured.get("members") != competitor["members"]
                or query["query_gene"] not in measured.get("queries", [])):
            raise ValueError("Measurement source/query ownership differs")
        rows = measured["measurements"]
        required = set(competitor["members"]) | {query["query_gene"]}
        absent = []
        for gene in required:
            row = rows.get(gene)
            if row is None or not row["reported"]:
                absent.append(gene)
                continue
            for key in ("full_sequence_score", "sequence_evalue"):
                value = float(row[key])
                if not math.isfinite(value) or (key == "sequence_evalue" and value < 0):
                    raise ValueError("Nonfinite or invalid backend measurement")
        if absent:
            missing.append({"profile_id": family, "reason": "incomplete_query_or_source", "genes": sorted(absent)})
            continue
        threshold = min(float(rows[gene]["full_sequence_score"]) for gene in competitor["members"])
        if not measured["calibration_complete"] or measured["minimum_member_score"] != threshold:
            raise ValueError("Backend minimum differs from complete source calibration")
        values[family] = {"score": float(rows[query["query_gene"]]["full_sequence_score"]),
                          "evalue": float(rows[query["query_gene"]]["sequence_evalue"]), "threshold": threshold}
    if missing:
        return {"action": "preserve", "reason": "incomplete_measurements", "missing": missing}
    top = max(value["score"] for value in values.values())
    winners = [c for c in competitors if values[c["profile_id"]]["score"] == top]
    if len(winners) != 1:
        return {"action": "reject", "reason": "top_score_tie"}
    winner = winners[0]
    family, anchor = winner["profile_id"], winner["anchor"]
    value = values[family]
    reason = ("own_family_winner" if winner["is_own"] else
              "sequence_evalue_failed" if value["evalue"] >= EVALUE else
              "below_member_threshold" if value["score"] < value["threshold"] else
              "no_positive_sequence_anchor" if anchor is None else None)
    record = {"winner_profile_id": family, **value}
    if reason:
        return {"action": "reject", "reason": reason, **record}
    if (anchor["gene"] not in winner["members"] or anchor["gene"] == query["query_gene"]
            or not math.isfinite(float(anchor["score"])) or anchor["score"] <= 0):
        raise ValueError("Invalid original sequence anchor")
    return {"action": "recruit", "reason": "unique_external_winner_passed", "anchor": dict(anchor), **record}


def replace_profile_anchors(gene_names, legacy, selection, decisions, *, measurement_only=False):
    """Replace selected directed anchors, retaining unselected reciprocal links."""
    names = list(gene_names)
    selected = {row["query_id"]: row for row in selection["selected"]}
    if len(selected) != len(selection["selected"]) or set(decisions) != set(selected):
        raise ValueError("Decisions must cover the complete selected query set")
    for query, row in selected.items():
        if not 0 <= query < len(names) or names[query] != row["query_gene"]:
            raise ValueError("Selected query ownership differs")
    current = {}
    legacy_genes = _indices(legacy.genes, len(names), "legacy query")
    legacy_members = _indices(legacy.members, len(names), "legacy anchor")
    legacy_weights = _numbers(legacy.anchor_scores, len(legacy_genes), "legacy anchor scores")
    if len(legacy_members) != len(legacy_genes):
        raise ValueError("Original directed anchor arrays differ in length")
    for query, member, weight in zip(legacy_genes, legacy_members, legacy_weights):
        query, member, weight = int(query), int(member), float(weight)
        if (query in current or not 0 <= query < len(names) or not 0 <= member < len(names)
                or query == member or not math.isfinite(weight) or weight <= 0):
            raise ValueError("Invalid original directed profile anchors")
        current[query] = (member, weight)
    for query, decision in decisions.items():
        action = decision["action"]
        if action not in ("preserve", "reject", "recruit"):
            raise ValueError("Unknown selective profile decision")
        if action == "recruit":
            anchor = decision["anchor"]
            family = next((c for c in selected[query]["competitors"] if c["profile_id"] == decision["winner_profile_id"]), None)
            if family is None or family["is_own"] or family["anchor"] != anchor:
                raise ValueError("Replacement anchor differs from nominated original support")
            member, weight = anchor["gene_id"], float(anchor["score"])
            if (not 0 <= member < len(names) or names[member] != anchor["gene"] or member == query
                    or not math.isfinite(weight) or weight <= 0):
                raise ValueError("Invalid replacement sequence anchor")
        if measurement_only or action == "preserve":
            continue
        current.pop(query, None)
        if action == "recruit":
            current[query] = (member, weight)
    ordered = sorted(current)
    return deduplicate_undirected_edges(names, ordered, [current[q][0] for q in ordered],
                                        [current[q][1] for q in ordered])


def score_selected_builtin(selection, profiles, gene_names, sequence_database, self_thresholds, cpu):
    """Measure the identical complete target pool with the original CPU scorer."""
    from .evalue import batch_estimate_evalues
    from .matrices import get_ka_params
    from .profile_expansion import _pack_profiles
    from .viterbi import is_c_available, batch_viterbi_multipair_c, batch_viterbi_c

    if isinstance(cpu, bool) or not isinstance(cpu, Integral) or not 1 <= cpu <= 16:
        raise ValueError("Require a bounded positive CPU budget")
    if list(sequence_database.ids) != list(gene_names) or not is_c_available():
        raise ValueError("Require complete original sequence ownership and original CPU scoring library")
    chosen = {p: profiles[p] for p in selection["selected_profiles"]}
    if not chosen:
        return {}
    ids, packed = _pack_profiles(chosen)
    lengths, offsets, emissions, _, inserts, transitions = packed
    gene_ids = {gene: index for index, gene in enumerate(gene_names)}
    query_sets = {int(p): sorted(row["query_gene"] for row in selection["selected"]
        if any(c["profile_id"] == int(p) for c in row["competitors"])) for p in ids}
    target_sets = {int(p): sorted(set(chosen[int(p)].member_ids) | set(query_sets[int(p)])) for p in ids}
    pairs = np.array([(local, gene_ids[g]) for local, p in enumerate(ids) for g in target_sets[int(p)]], dtype=np.int32)
    arguments = (emissions, inserts, transitions, offsets, lengths, sequence_database.flat_sequences,
                 sequence_database.offsets, sequence_database.lengths, pairs, 64)
    try:
        scores = batch_viterbi_multipair_c(*arguments, n_threads=int(cpu))
    except RuntimeError:
        scores = batch_viterbi_c(*arguments, n_threads=int(cpu))
    lam, k = get_ka_params("BLOSUM62")
    evalues = batch_estimate_evalues(scores, lengths[pairs[:, 0]], int(sequence_database.lengths.sum()), lam, k)
    result = {int(p): {"status": "builtin_profile_scored", "members": list(chosen[int(p)].member_ids),
        "queries": query_sets[int(p)], "calibration_complete": True, "measurements": {},
        "score_unit": "raw native integer profile score, not HMMER bits",
        "parameters": {"database_genes": len(gene_names), "database_residues": int(sequence_database.lengths.sum()),
                       "band_width": 64, "matrix": "BLOSUM62", "cpu": int(cpu)}} for p in ids}
    for (local, gene), score, evalue in zip(pairs, scores, evalues):
        result[int(ids[local])]["measurements"][gene_names[int(gene)]] = {
            "reported": True, "full_sequence_score": float(score), "sequence_evalue": float(evalue)}
    for family, measured in result.items():
        if chosen[family].projection_policy != "legacy":
            measured["parameters"]["profile_projection_policy"] = chosen[family].projection_policy
        threshold = min(measured["measurements"][g]["full_sequence_score"] for g in measured["members"])
        if threshold != self_thresholds[family]:
            raise ValueError("Recomputed native source minimum differs from original self threshold")
        measured["minimum_member_score"] = threshold
    return result


def selective_profile_edges(gene_names, clusters, profiles, profile_hits, self_thresholds,
                            hit_queries, hit_targets, hit_scores, sequence_database,
                            *, mode="off", output_directory=None, cpu=1):
    """Run a prespecified experimental arm and return ordinary profile edges.

    Downstream clustering, satellite candidates and final native reconciliation
    remain the caller's unchanged responsibility. An active call requires a
    fresh audit directory; errors retain evidence and are never inactivity.
    """
    from ..phylogeny_pipeline import _atomic_write
    from .profile_expansion import build_strict_profile_edges, _select_strict_profile_anchors, _decode_sequence
    from .selective_hmmer import score_family_hmmer, _record

    if mode not in ("off", "builtin_reselect", "hmmer_score_only", "hmmer_reselect"):
        raise ValueError("Unknown selective profile mode")
    if mode == "off":
        edges = build_strict_profile_edges(gene_names, clusters, profile_hits, self_thresholds,
                                          hit_queries, hit_targets, hit_scores)
        return edges, {"mode": "off", "status": "unchanged", "orthology_inferred": False}
    if (output_directory is None or isinstance(cpu, bool) or not isinstance(cpu, Integral)
            or not 1 <= cpu <= 16 or os.environ.get("CUDA_VISIBLE_DEVICES") != ""):
        raise ValueError("Active selective profile pilot requires fresh evidence, CPU-only visibility and1..16 CPUs")
    directory = Path(output_directory)
    directory.mkdir(parents=True, exist_ok=False)
    report = {"schema_version": 1, "mode": mode, "status": "running", "orthology_inferred": False}

    def save():
        _atomic_write(directory / "selection.json", json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n")

    save()
    try:
        if list(sequence_database.ids) != list(gene_names):
            raise ValueError("Complete sequence database differs from original gene table")
        selection = nominate_queries(gene_names, clusters, profiles, profile_hits, self_thresholds,
                                     hit_queries, hit_targets, hit_scores)
        report["selection"] = selection
        old = _select_strict_profile_anchors(gene_names, clusters, profile_hits, self_thresholds,
                                            hit_queries, hit_targets, hit_scores)
        selected_ids = {row["query_id"] for row in selection["selected"]}
        report["original_selected_anchors"] = [{"query_id": int(q), "query_gene": gene_names[int(q)],
            "profile_id": int(p), "anchor_gene_id": int(member), "anchor_gene": gene_names[int(member)],
            "anchor_score": float(weight), "profile_score": float(score)}
            for q, p, member, weight, score in zip(old.genes, old.profiles, old.members, old.anchor_scores, old.profile_scores)
            if int(q) in selected_ids]
        save()
        if mode == "builtin_reselect":
            measurements = score_selected_builtin(selection, profiles, gene_names, sequence_database, self_thresholds, cpu)
            _atomic_write(directory / "builtin_measurements.json", json.dumps(measurements, indent=2, sort_keys=True, allow_nan=False) + "\n")
            report["measurements"] = {str(p): _record(directory / "builtin_measurements.json") for p in measurements}
        else:
            lookup = {gene: index for index, gene in enumerate(gene_names)}

            def measure(family):
                profile = profiles[family]
                query_names = sorted(row["query_gene"] for row in selection["selected"]
                    if any(c["profile_id"] == family for c in row["competitors"]))
                names = set(profile.member_ids) | set(query_names)
                sequences = {g: _decode_sequence(sequence_database.get_sequence(lookup[g])) for g in names}
                return family, score_family_hmmer(profile, {g: sequences[g] for g in profile.member_ids},
                    {g: sequences[g] for g in query_names}, len(gene_names), directory / f"profile{family:08d}")

            with ThreadPoolExecutor(max_workers=int(cpu)) as executor:
                measurements = dict(executor.map(measure, selection["selected_profiles"]))
            report["measurements"] = {str(p): _record(directory / f"profile{p:08d}" / "measurements.json") for p in measurements}
        decisions = {row["query_id"]: decide_query(row, measurements) for row in selection["selected"]}
        report["decisions"] = [{"query_id": q, "query_gene": gene_names[q], **decision} for q, decision in decisions.items()]
        previous = {int(q): (int(m), float(w)) for q, m, w in zip(old.genes, old.members, old.anchor_scores)}
        changed = []
        for query, decision in decisions.items():
            after = (previous.get(query) if decision["action"] == "preserve" else
                     (decision["anchor"]["gene_id"], decision["anchor"]["score"])
                     if decision["action"] == "recruit" else None)
            if after != previous.get(query):
                changed.append(gene_names[query])
        edges = replace_profile_anchors(gene_names, old, selection, decisions, measurement_only=mode == "hmmer_score_only")
        report.update(status="selective_profile_edges_built", selected_queries=len(selected_ids),
            completely_measured_queries=sum(d["action"] != "preserve" for d in decisions.values()),
            counterfactual_changed_query_anchors=sorted(changed),
            applied_changed_query_anchors=[] if mode == "hmmer_score_only" else sorted(changed),
            final_profile_edges=len(edges), measurement_only=mode == "hmmer_score_only")
    except Exception as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        save()
    return edges, report
