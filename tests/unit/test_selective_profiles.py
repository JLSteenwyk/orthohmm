from copy import deepcopy
from dataclasses import replace
import json
import shutil
from types import SimpleNamespace

import numpy as np
import pytest

from orthohmm.search.profile_expansion import ProfileHits, _select_strict_profile_anchors
from orthohmm.search.selective_profiles import (
    nominate_queries, decide_query, replace_profile_anchors, score_selected_builtin, selective_profile_edges,
)


def fixture(externals=(0, 1), *, own=False):
    names = ["a2", "a1", "a3", "b1", "b2", "b3", "query", "query2", "query3"]
    groups = [[0, 1, 2], [3, 4, 5], [6, 7, 8]] if own else [[0, 1, 2], [3, 4, 5], [6], [7], [8]]
    profiles = {p: SimpleNamespace(member_ids=[names[g] for g in group], sequence_weighting="none")
                for p, group in enumerate(groups) if len(group) >= 3}
    hits = ProfileHits(np.array(externals, dtype=np.int32), np.full(len(externals), 6, dtype=np.int32),
                       np.array([90. - p for p in range(len(externals))]), np.full(len(externals), 1e-10), len(externals))
    return [names, groups, profiles, hits, {p: 100. for p in profiles},
            np.array([6, 6, 6, 0, 1], dtype=np.int32), np.array([0, 1, 3, 6, 6], dtype=np.int32),
            np.array([10., 10., 8., 100., 100.])]


def nominate(args=None):
    return nominate_queries(*(fixture() if args is None else args))


def backend(query, scores=None):
    scores = scores or {c["profile_id"]: 110. - i for i, c in enumerate(query["competitors"])}
    return {c["profile_id"]: {"status": "builtin_profile_scored", "members": list(c["members"]),
        "queries": [query["query_gene"]], "calibration_complete": True, "minimum_member_score": 100.,
        "measurements": {gene: {"reported": True, "full_sequence_score": scores[c["profile_id"]]
            if gene == query["query_gene"] else 100., "sequence_evalue": 1e-10}
            for gene in set(c["members"]) | {query["query_gene"]}}} for c in query["competitors"]}


def legacy(args):
    return _select_strict_profile_anchors(args[0], args[1], args[3], args[4], *args[5:])


def test_complete_competition_priority_anchors_and_readonly_inputs():
    args = fixture()
    before = deepcopy(args)
    result = nominate(args)
    row, = result["selected"]
    assert row["query_gene"] == "query" and row["own_profile_id"] is None
    assert row["priority"] == .1 and row["priority_profile_id"] == 0
    assert [c["profile_id"] for c in row["competitors"]] == [0, 1]
    assert row["competitors"][0]["anchor"] == {"gene_id": 0, "gene": "a2", "score": 10.}
    assert row["competitors"][1]["anchor"]["score"] == 8.
    assert not result["orthology_inferred"] and len(result["candidate_digest"]) == 64
    assert args[:3] == before[:3]
    assert all(np.array_equal(a, b) for a, b in zip(args[5:], before[5:]))
    json.dumps(result, allow_nan=False)


def test_duplicate_profile_rows_use_maximum_and_hits_are_order_invariant():
    args = fixture()
    args[3] = ProfileHits(np.array([1, 0, 0]), np.array([6, 6, 6]), np.array([89., 70., 90.]), np.full(3, 1e-8), 3)
    original = nominate(args)
    args[3] = replace(args[3], profile_cluster_ids=args[3].profile_cluster_ids[::-1],
                      gene_ids=args[3].gene_ids[::-1], scores=args[3].scores[::-1], evalues=args[3].evalues[::-1])
    for index in (5, 6, 7):
        args[index] = args[index][::-1]
    assert nominate(args) == original


def test_global_and_family_order_do_not_change_candidate_content_digest():
    args = fixture(own=True)
    before = nominate(args)
    permutation = [8, 6, 3, 1, 0, 7, 2, 5, 4]
    new_id = {old: new for new, old in enumerate(permutation)}
    family_order = [2, 1, 0]
    family_id = {old: new for new, old in enumerate(family_order)}
    args[0] = [args[0][g] for g in permutation]
    args[1] = [[new_id[g] for g in args[1][old]] for old in family_order]
    args[2] = {family_id[p]: profile for p, profile in args[2].items()}
    args[3] = replace(args[3], gene_ids=np.array([new_id[int(g)] for g in args[3].gene_ids]),
                      profile_cluster_ids=np.array([family_id[int(p)] for p in args[3].profile_cluster_ids]))
    args[4] = {family_id[p]: value for p, value in args[4].items()}
    args[5], args[6] = [np.array([new_id[int(g)] for g in values]) for values in args[5:7]]
    after = nominate(args)
    assert after["candidate_digest"] == before["candidate_digest"]
    assert after["selected"][0]["competitors"][0]["members"] == before["selected"][0]["competitors"][0]["members"]


def test_own_profile_is_complete_competitor_but_small_family_has_none():
    row = nominate(fixture(own=True))["selected"][0]
    assert row["own_profile_id"] == 2
    own, = [c for c in row["competitors"] if c["is_own"]]
    assert own["members"] == ["query", "query2", "query3"] and own["anchor"] is None
    assert own["original_hit_score"] is None
    assert len(nominate()["selected"][0]["competitors"]) == 2


@pytest.mark.parametrize("kind,status", [("one_above", "not_ambiguous"),
    ("one_equal", "not_ambiguous"), ("no_anchor", "excluded_no_positive_sequence_anchor")])
def test_ineligible_queries_preserve_explicit_reasons(kind, status):
    args = fixture(externals=(0,))
    if kind == "one_above":
        args[4][0] = 89.
    elif kind == "one_equal":
        args[4][0] = 90.
    else:
        args[7][:] = 0.
    result = nominate(args)
    assert not result["selected"] and result["audit"][0]["status"] == status


def test_single_below_self_threshold_is_eligible_and_reverse_only_anchor_works():
    args = fixture(externals=(0,))
    args[5], args[6], args[7] = np.array([0]), np.array([6]), np.array([11.])
    row, = nominate(args)["selected"]
    assert row["competitors"][0]["anchor"] == {"gene_id": 0, "gene": "a2", "score": 11.}


def test_unanchored_competitor_is_retained_and_can_block_recruitment():
    args = fixture()
    args[7][2] = 0.
    row, = nominate(args)["selected"]
    assert row["competitors"][1]["anchor"] is None
    decision = decide_query(row, backend(row, {0: 110., 1: 120.}))
    assert decision["action"] == "reject" and decision["reason"] == "no_positive_sequence_anchor"


def large_fixture(*, distinct=False):
    families = 132 if distinct else 2
    names = [f"f{family:03d}_{j}" for family in range(families) for j in range(3)]
    groups = [list(range(3 * p, 3 * p + 3)) for p in range(families)]
    profiles = {p: SimpleNamespace(member_ids=[names[g] for g in group], sequence_weighting="none") for p, group in enumerate(groups)}
    query_count = 66 if distinct else 130
    first = len(names)
    names.extend(f"q{q:03d}" for q in range(query_count))
    groups.extend([[first + q] for q in range(query_count)])
    pids = [2 * q + p if distinct else p for q in range(query_count) for p in (0, 1)]
    qids = [first + q for q in range(query_count) for p in (0, 1)]
    hits = ProfileHits(np.array(pids), np.array(qids), np.tile([90., 80.], query_count), np.full(2 * query_count, 1e-10), len(qids))
    return [names, groups, profiles, hits, {p: 100. for p in profiles},
            np.array(qids), np.array([3 * p for p in pids]), np.full(len(qids), 10.)]


def test_query_budget_keeps_128_complete_queries_and_canonical_name_order():
    result = nominate(large_fixture())
    assert len(result["selected"]) == 128 and len(result["selected_profiles"]) == 2
    assert [r["query_gene"] for r in result["selected"]] == [f"q{q:03d}" for q in range(128)]
    assert [r["status"] for r in result["audit"]][-2:] == ["excluded_query_budget"] * 2


def test_distinct_profile_budget_excludes_whole_query_not_competitors():
    result = nominate(large_fixture(distinct=True))
    assert len(result["selected"]) == 64 and len(result["selected_profiles"]) == 128
    assert all(len(r["competitors"]) == 2 for r in result["selected"])
    assert [r["status"] for r in result["audit"]][-2:] == ["excluded_profile_budget"] * 2


def test_more_than_four_external_profiles_excludes_query_without_truncation():
    args = large_fixture(distinct=True)
    args[3] = ProfileHits(np.arange(5), np.full(5, len(args[2]) * 3), np.full(5, 90.), np.full(5, 1e-9), 5)
    result = nominate(args)
    assert not result["selected"] and len(result["audit"][0]["external_profiles"]) == 5
    assert result["audit"][0]["status"] == "excluded_complete_competition_over_budget"


def test_four_external_profiles_plus_own_retains_all_five_competitors():
    args = large_fixture(distinct=True)
    args[3] = ProfileHits(np.arange(1, 5), np.zeros(4, dtype=int), np.full(4, 90.), np.full(4, 1e-9), 4)
    args[5], args[6], args[7] = np.zeros(4, dtype=int), np.arange(3, 15, 3), np.full(4, 10.)
    row, = nominate(args)["selected"]
    assert row["own_profile_id"] == 0 and len(row["competitors"]) == 5
    assert sum(c["is_own"] for c in row["competitors"]) == 1
    assert [c["profile_id"] for c in row["competitors"]] == list(range(5))


def test_tied_external_priority_uses_canonical_members_not_hit_row_order():
    args = fixture()
    args[3].scores[:] = 90.
    args[4][1] = 200.
    before = nominate(args)
    assert before["selected"][0]["priority"] == .1
    args[3] = replace(args[3], profile_cluster_ids=args[3].profile_cluster_ids[::-1])
    assert nominate(args) == before


@pytest.mark.parametrize("kind", ["duplicates", "overlap", "incomplete", "profile_members", "weighted", "missing_threshold",
    "nan_threshold", "float_query", "bad_target", "nan_score", "nonsignificant", "nonfinite_evalue", "length"])
def test_invalid_nomination_inputs_fail_before_selection(kind):
    args = fixture()
    if kind == "duplicates": args[0][1] = args[0][0]
    elif kind == "overlap": args[1][1][0] = 0
    elif kind == "incomplete": args[1].pop()
    elif kind == "profile_members": args[2][0].member_ids.reverse()
    elif kind == "weighted": args[2][0].sequence_weighting = "henikoff"
    elif kind == "missing_threshold": args[4].pop(0)
    elif kind == "nan_threshold": args[4][0] = float("nan")
    elif kind == "float_query": args[5] = np.array([6.5])
    elif kind == "bad_target": args[6][0] = 100
    elif kind == "nan_score": args[7][0] = float("nan")
    elif kind == "nonsignificant": args[3].evalues[0] = 1e-4
    elif kind == "nonfinite_evalue": args[3].evalues[0] = float("inf")
    else: args[6] = args[6][:-1]
    with pytest.raises(ValueError): nominate(args)


def test_inactive_profile_pool_is_explicit_and_json_safe():
    args = fixture()
    args[2], args[4] = {}, {}
    args[3] = ProfileHits(np.array([], dtype=int), np.array([], dtype=int), np.array([]), np.array([]), 0)
    result = nominate(args)
    assert not result["selected"] and not result["audit"]
    assert result["without_external_hits"] == sorted(args[0])
    json.dumps(result, allow_nan=False)


@pytest.mark.parametrize("kind,reason", [("tie", "top_score_tie"), ("own", "own_family_winner"),
    ("evalue", "sequence_evalue_failed"), ("threshold", "below_member_threshold")])
def test_winner_precedes_gates_and_never_falls_back_to_weaker_family(kind, reason):
    row = nominate(fixture(own=kind == "own"))["selected"][0]
    measurements = backend(row)
    if kind == "tie":
        measurements = backend(row, {0: 110., 1: 110.})
    elif kind == "own":
        measurements = backend(row, {0: 110., 1: 109., 2: 120.})
    elif kind == "evalue":
        measurements[0]["measurements"]["query"]["sequence_evalue"] = 1e-4
    else:
        measurements = backend(row, {0: 99., 1: 98.})
        measurements[1]["minimum_member_score"] = 90.
        for gene in row["competitors"][1]["members"]:
            measurements[1]["measurements"][gene]["full_sequence_score"] = 90.
    decision = decide_query(row, measurements)
    assert decision["action"] == "reject" and decision["reason"] == reason


def test_score_at_member_minimum_passes_but_evalue_must_be_strict():
    row = nominate()["selected"][0]
    decision = decide_query(row, backend(row, {0: 100., 1: 99.}))
    assert decision["action"] == "recruit" and decision["anchor"]["score"] == 10.
    assert decision["score"] == 100.  # Never used as the graph edge weight.


@pytest.mark.parametrize("kind", ["profile", "competitor", "source"])
def test_missing_measurement_preserves_legacy_instead_of_dropping_competitor(kind):
    row = nominate()["selected"][0]
    measurements = backend(row)
    if kind == "profile": measurements.pop(1)
    elif kind == "competitor": measurements[1]["measurements"]["query"] = {"reported": False}
    else: measurements[1]["measurements"]["b2"] = {"reported": False}
    assert decide_query(row, measurements)["action"] == "preserve"


@pytest.mark.parametrize("kind", ["failed", "nan", "negative_evalue", "wrong_members", "wrong_query", "wrong_min", "calibration"])
def test_malformed_backend_is_not_a_negative_hit(kind):
    row = nominate()["selected"][0]
    measurements = backend(row)
    m = measurements[1]
    if kind == "failed": m["status"] = "failed"
    elif kind == "nan": m["measurements"]["query"]["full_sequence_score"] = float("nan")
    elif kind == "negative_evalue": m["measurements"]["query"]["sequence_evalue"] = -1.
    elif kind == "wrong_members": m["members"].reverse()
    elif kind == "wrong_query": m["queries"] = ["other"]
    elif kind == "wrong_min": m["minimum_member_score"] = 50.
    else: m["calibration_complete"] = False
    with pytest.raises(ValueError): decide_query(row, measurements)


def test_replacement_operates_on_directed_anchors_and_preserves_unselected_reciprocal():
    args = fixture()
    selection = nominate(args)
    # Original anchor from query and an independently proposed reverse anchor.
    anchors = SimpleNamespace(genes=np.array([0, 6]), members=np.array([6, 0]), anchor_scores=np.array([5., 10.]))
    decision = {6: {"action": "reject", "reason": "top_score_tie"}}
    edge = replace_profile_anchors(args[0], anchors, selection, decision)
    assert list(zip(edge.sources, edge.targets, edge.weights)) == [(0, 6, 5.)]
    control = replace_profile_anchors(args[0], anchors, selection, decision, measurement_only=True)
    assert list(zip(control.sources, control.targets, control.weights)) == [(0, 6, 10.)]
    decision[6] = {"action": "preserve", "reason": "incomplete_measurements"}
    unchanged = replace_profile_anchors(args[0], anchors, selection, decision)
    assert np.array_equal(unchanged.weights, control.weights)


def test_replacement_uses_original_sequence_weight_and_rejects_changed_anchor():
    args = fixture()
    selection = nominate(args)
    row, = selection["selected"]
    decisions = {6: decide_query(row, backend(row))}
    edge = replace_profile_anchors(args[0], legacy(args), selection, decisions)
    assert list(zip(edge.sources, edge.targets, edge.weights)) == [(0, 6, 10.)]
    decisions[6]["anchor"]["score"] = decisions[6]["score"]
    with pytest.raises(ValueError, match="nominated original support"):
        replace_profile_anchors(args[0], legacy(args), selection, decisions)


def test_legacy_anchor_selected_by_actual_original_helper_is_unchanged_in_measurement_control():
    args = fixture()
    args[4][0] = 90.
    selection = nominate(args)
    old = legacy(args)
    assert list(old.genes) == [6] and list(old.members) == [0] and list(old.anchor_scores) == [10.]
    edge = replace_profile_anchors(args[0], old, selection, {6: {"action": "reject"}}, measurement_only=True)
    assert list(zip(edge.sources, edge.targets, edge.weights)) == [(0, 6, 10.)]


def test_actual_hmmer_backend_through_nomination_decision_and_graph_anchor(tmp_path):
    from orthohmm.search.matrices import get_matrix, get_background_freqs
    from orthohmm.search.msa_profile import build_msa_profile
    from orthohmm.search.selective_hmmer import score_family_hmmer
    if not shutil.which("hmmbuild") or not shutil.which("hmmsearch"):
        pytest.skip("Requires real HMMER3.4 tools")
    args = fixture()
    rng = np.random.default_rng(20261018)
    alphabet = np.array(list("ACDEFGHIKLMNPQRSTVWY"))
    ancestor = rng.choice(alphabet, 180)
    sequences = {}
    for family in (0, 1):
        base = ancestor.copy()
        if family:
            positions = rng.choice(len(base), 60, replace=False)
            base[positions] = rng.choice(alphabet, len(positions))
        for gene in args[2][family].member_ids:
            variant = base.copy()
            positions = rng.choice(len(base), 4, replace=False)
            variant[positions] = rng.choice(alphabet, len(positions))
            sequences[gene] = "".join(variant)
        members = list(args[2][family].member_ids)
        args[2][family] = build_msa_profile([sequences[g] for g in members], members, get_matrix("BLOSUM62"), get_background_freqs("BLOSUM62"))
    selection = nominate(args)
    query, = selection["selected"]
    measurements = {family: score_family_hmmer(profile, {g: sequences[g] for g in profile.member_ids},
        {"query": sequences["a2"]}, 10000, tmp_path / f"family{family}") for family, profile in args[2].items()}
    decision = decide_query(query, measurements)
    assert decision["action"] == "recruit" and decision["winner_profile_id"] == 0
    assert decision["score"] != decision["anchor"]["score"]
    edges = replace_profile_anchors(args[0], legacy(args), selection, {6: decision})
    assert list(zip(edges.sources, edges.targets, edges.weights)) == [(0, 6, 10.)]
    json.dumps({"selection": selection, "decision": decision}, allow_nan=False)


@pytest.fixture
def real_stage(monkeypatch):
    from orthohmm.search.matrices import get_matrix, get_background_freqs
    from orthohmm.search.msa_profile import build_msa_profile
    from orthohmm.search.sequences import encode_sequence, SpeciesSequences
    from orthohmm.search.profile_expansion import compute_profile_self_thresholds
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "")
    monkeypatch.setattr("orthohmm.search.profile_expansion.is_cuda_available", lambda: False)
    args = fixture()
    rng = np.random.default_rng(20261019)
    alphabet = np.array(list("ACDEFGHIKLMNPQRSTVWY"))
    base = rng.choice(alphabet, 160)
    sequences = {}
    for family, profile in args[2].items():
        ancestor = base.copy()
        if family:
            positions = rng.choice(len(base), 55, replace=False)
            ancestor[positions] = rng.choice(alphabet, len(positions))
        for gene in profile.member_ids:
            variant = ancestor.copy()
            positions = rng.choice(len(base), 3, replace=False)
            variant[positions] = rng.choice(alphabet, len(positions))
            sequences[gene] = "".join(variant)
        args[2][family] = build_msa_profile([sequences[g] for g in profile.member_ids], profile.member_ids,
                                           get_matrix("BLOSUM62"), get_background_freqs("BLOSUM62"))
    for gene in args[0][6:]:
        sequences[gene] = sequences["a2"]
    encoded = [encode_sequence(sequences[g]) for g in args[0]]
    lengths = np.array([len(s) for s in encoded], dtype=np.int32)
    offsets = np.concatenate(([0], np.cumsum(lengths[:-1], dtype=np.int64)))
    database = SpeciesSequences("all_genes", args[0], np.concatenate(encoded), offsets, lengths)
    args[4] = compute_profile_self_thresholds(args[2], args[0], database, 1)
    return args, database


def test_real_native_rescoring_reconstructs_source_threshold_and_uses_full_residues(real_stage):
    from orthohmm.search.evalue import estimate_evalue
    from orthohmm.search.matrices import get_ka_params
    args, database = real_stage
    selection = nominate(args)
    scored = score_selected_builtin(selection, args[2], args[0], database, args[4], 1)
    assert len(scored) == 2
    for family, measured in scored.items():
        assert measured["minimum_member_score"] == args[4][family]
        assert set(measured["measurements"]) == set(args[2][family].member_ids) | {"query"}
        value = measured["measurements"]["query"]
        expected = estimate_evalue(int(value["full_sequence_score"]), args[2][family].length,
                                   int(database.lengths.sum()), *get_ka_params("BLOSUM62"))
        assert value["sequence_evalue"] == expected
    assert decide_query(selection["selected"][0], scored)["action"] == "recruit"
    args[4][0] += 1
    with pytest.raises(ValueError, match="original self threshold"):
        score_selected_builtin(selection, args[2], args[0], database, args[4], 1)


def test_off_bypasses_nomination_measurement_and_artifact_creation(tmp_path, monkeypatch):
    args = fixture()
    args[4][0] = 90.
    monkeypatch.setattr("orthohmm.search.selective_profiles.nominate_queries", lambda *a: pytest.fail("Off must bypass selection"))
    directory = tmp_path / "off"
    edges, report = selective_profile_edges(*args, None, output_directory=directory)
    assert report["status"] == "unchanged" and not directory.exists()
    assert list(zip(edges.sources, edges.targets, edges.weights)) == [(0, 6, 10.)]


@pytest.mark.parametrize("mode", ["builtin_reselect", "hmmer_score_only", "hmmer_reselect"])
def test_actual_experimental_stage_records_nonempty_decisions_and_retains_original_scale(real_stage, tmp_path, mode):
    args, database = real_stage
    if mode.startswith("hmmer") and (not shutil.which("hmmbuild") or not shutil.which("hmmsearch")):
        pytest.skip("Requires real HMMER3.4 tools")
    before = deepcopy(args[5:])
    directory = tmp_path / mode
    edges, report = selective_profile_edges(*args, database, mode=mode, output_directory=directory, cpu=2)
    assert report["status"] == "selective_profile_edges_built"
    assert report["selected_queries"] == report["completely_measured_queries"] == 1
    assert not report["orthology_inferred"]
    assert len(report["measurements"]) == 2 and report["decisions"][0]["action"] == "recruit"
    assert all(np.array_equal(a, b) for a, b in zip(args[5:], before))
    assert json.loads((directory / "selection.json").read_text()) == report
    if mode == "hmmer_score_only":
        off, _ = selective_profile_edges(*args, database)
        assert np.array_equal(edges.sources, off.sources) and np.array_equal(edges.targets, off.targets)
        assert edges.weights.tobytes() == off.weights.tobytes()
        assert not report["applied_changed_query_anchors"]
    else:
        assert list(zip(edges.sources, edges.targets, edges.weights)) == [(0, 6, 10.)]
    with pytest.raises(FileExistsError):
        selective_profile_edges(*args, database, mode=mode, output_directory=directory, cpu=2)


def test_tool_failure_is_preserved_not_reinterpreted_as_inactivity(real_stage, tmp_path, monkeypatch):
    args, database = real_stage
    def broken(*a, **kw):
        raise RuntimeError("controlled backend failure")
    monkeypatch.setattr("orthohmm.search.selective_hmmer.score_family_hmmer", broken)
    directory = tmp_path / "failed"
    with pytest.raises(RuntimeError, match="controlled backend failure"):
        selective_profile_edges(*args, database, mode="hmmer_reselect", output_directory=directory)
    report = json.loads((directory / "selection.json").read_text())
    assert report["status"] == "failed" and report["selection"]["selected"]
    assert "decisions" not in report and "completely_measured_queries" not in report


def test_missing_query_measurements_preserve_actual_legacy_stage_edges(real_stage, tmp_path, monkeypatch):
    args, database = real_stage
    # Ensure there is an actual legacy directed edge to preserve.
    args[3] = replace(args[3], scores=np.array([args[4][0] + 5., args[4][1] - 5.]))
    selection = nominate(args)
    measured = score_selected_builtin(selection, args[2], args[0], database, args[4], 1)
    measured[1]["measurements"]["query"] = {"reported": False}
    monkeypatch.setattr("orthohmm.search.selective_profiles.score_selected_builtin", lambda *a: measured)
    edges, report = selective_profile_edges(*args, database, mode="builtin_reselect", output_directory=tmp_path / "missing")
    off, _ = selective_profile_edges(*args, database)
    assert edges.weights.tobytes() == off.weights.tobytes()
    assert edges.sources.tobytes() == off.sources.tobytes() and edges.targets.tobytes() == off.targets.tobytes()
    assert report["decisions"][0]["action"] == "preserve" and report["completely_measured_queries"] == 0
    assert not report["applied_changed_query_anchors"]


@pytest.mark.parametrize("cpu", [0, 17, True, 1.5])
def test_resource_guards_reject_before_creating_artifacts(tmp_path, cpu, monkeypatch):
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "")
    directory = tmp_path / "invalid"
    with pytest.raises(ValueError):
        selective_profile_edges(*fixture(), None, mode="builtin_reselect", output_directory=directory, cpu=cpu)
    assert not directory.exists()


def test_empty_active_stage_is_explicit_not_a_fabricated_accuracy_comparison(real_stage, tmp_path):
    args, database = real_stage
    args[3] = ProfileHits(np.array([], dtype=int), np.array([], dtype=int), np.array([]), np.array([]), 0)
    edges, report = selective_profile_edges(*args, database, mode="builtin_reselect", output_directory=tmp_path / "empty")
    assert len(edges) == report["selected_queries"] == report["completely_measured_queries"] == 0
    assert not report["decisions"] and not report["measurements"]
