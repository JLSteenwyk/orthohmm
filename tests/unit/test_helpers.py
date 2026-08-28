"""Unit tests for orthohmm/helpers.py (pure-function utilities)."""
import multiprocessing
import os

import numpy as np
import pytest

from orthohmm.helpers import (
    IndexedEdgeThresholds,
    IndexedEdges,
    _best_indexed_hits,
    _determine_indexed_network_edges,
    StartStep,
    StopStep,
    SubstitutionMatrix,
    correct_by_phylogenetic_distance,
    determine_edge_thresholds,
    generate_orthogroup_clusters_file,
    generate_orthogroup_files,
    generate_phmmer_cmds,
    get_all_fasta_entries,
    get_best_hits_and_scores,
    get_orthogroup_information,
    get_sequence_lengths,
    get_singletons,
    get_threshold_per_gene,
    merge_with_gene_lengths,
    normalize_by_gene_length,
    process_pair_determine_network_edges,
    process_pair_edge_thresholds,
    read_and_filter_phmmer_output,
    update_progress,
)
from orthohmm.search.engine import IndexedSearchResults, IndexedSpeciesPairResults


# ─── Enums ──────────────────────────────────────────────────────────────


class TestEnums:
    def test_start_step_values(self):
        assert StartStep.search_res.value == "search_res"

    def test_stop_step_values(self):
        assert StopStep.prepare.value == "prepare"
        assert StopStep.infer.value == "infer"
        assert StopStep.write.value == "write"

    def test_substitution_matrix_supports_eleven_matrices(self):
        names = {m.value for m in SubstitutionMatrix}
        for n in ("BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80",
                  "BLOSUM90", "PAM30", "PAM70", "PAM120", "PAM240",
                  "WAG", "LG"):
            assert n in names


# ─── phmmer command generation ──────────────────────────────────────────


class TestGeneratePhmmerCmds:
    def test_pairwise_combinations_count(self):
        # N files -> N² self+cross commands (phmmer is asymmetric per query/target)
        cmds = generate_phmmer_cmds(
            ["a.fa", "b.fa", "c.fa"],
            phmmer="phmmer",
            output_directory="/out",
            fasta_directory="/in",
            cpu=4,
            substitution_matrix=SubstitutionMatrix.blosum62,
        )
        assert len(cmds) == 9  # 3x3 product including self-self

    def test_includes_matrix_and_cpu_flags(self):
        cmds = generate_phmmer_cmds(
            ["a.fa", "b.fa"],
            phmmer="/path/to/phmmer",
            output_directory="/out",
            fasta_directory="/in",
            cpu=8,
            substitution_matrix=SubstitutionMatrix.blosum80,
        )
        first = cmds[0]
        assert "/path/to/phmmer" in first
        assert "--mx BLOSUM80" in first
        assert "--cpu 8" in first
        assert "/out/orthohmm_working_res/a.fa_2_a.fa.phmmerout.txt" in first


# ─── get_sequence_lengths ───────────────────────────────────────────────


class TestGetSequenceLengths:
    def test_two_species_three_genes(self, tmp_path):
        (tmp_path / "h.fa").write_text(">h1\nACDE\n>h2\nACDEFG\n")
        (tmp_path / "m.fa").write_text(">m1\nMM\n")
        gl = get_sequence_lengths(str(tmp_path), ["h.fa", "m.fa"])
        # one row per gene with (spp, name, length)
        assert gl.shape == (3,)
        # The structured array has these fields:
        assert gl.dtype.names == ("spp", "name", "length")
        by_name = {row["name"]: int(row["length"]) for row in gl}
        assert by_name == {"h1": 4, "h2": 6, "m1": 2}

    def test_header_first_whitespace_token_only(self, tmp_path):
        (tmp_path / "h.fa").write_text(">geneA description text\nACDE\n")
        gl = get_sequence_lengths(str(tmp_path), ["h.fa"])
        assert gl[0]["name"] == "geneA"


# ─── get_singletons ─────────────────────────────────────────────────────


class TestGetSingletons:
    def _make_gene_lengths(self, *names):
        return np.array(
            [(n, 100) for n in names],
            dtype=[("name", object), ("length", int)],
        )

    def test_genes_outside_clusters_become_singletons(self):
        gl = self._make_gene_lengths("g1", "g2", "g3", "g4", "g5")
        clusters = [["g1", "g2"]]  # only 2 of 5 covered
        out_clusters, singletons = get_singletons(gl, clusters)
        # Stable ordering keeps orthogroup IDs reproducible across processes.
        assert singletons == [["g3"], ["g4"], ["g5"]]
        # The output clusters list contains the original cluster + 3 singletons
        assert len(out_clusters) == 1 + 3

    def test_all_genes_clustered_yields_no_singletons(self):
        gl = self._make_gene_lengths("g1", "g2", "g3")
        clusters = [["g1", "g2", "g3"]]
        out, singletons = get_singletons(gl, clusters)
        assert singletons == []
        assert out == [["g1", "g2", "g3"]]


# ─── get_all_fasta_entries ──────────────────────────────────────────────


class TestGetAllFastaEntries:
    def test_returns_dict_per_file(self, tmp_path):
        (tmp_path / "h.fa").write_text(">h1\nACDE\n>h2\nFG\nHI\n")
        (tmp_path / "m.fa").write_text(">m1\nKKK\n")
        entries = get_all_fasta_entries(str(tmp_path), ["h.fa", "m.fa"])
        assert set(entries.keys()) == {"h.fa", "m.fa"}
        # multi-line sequence concatenated
        assert entries["h.fa"] == {"h1": "ACDE", "h2": "FGHI"}
        assert entries["m.fa"] == {"m1": "KKK"}

    def test_header_description_dropped(self, tmp_path):
        (tmp_path / "h.fa").write_text(">geneA species:Homo sapiens\nAAA\n")
        entries = get_all_fasta_entries(str(tmp_path), ["h.fa"])
        assert list(entries["h.fa"].keys()) == ["geneA"]


# ─── get_orthogroup_information ─────────────────────────────────────────


class TestGetOrthogroupInformation:
    def _setup_two_species(self):
        gl = np.array(
            [
                ("h.fa", "h1", 100), ("h.fa", "h2", 110),
                ("m.fa", "m1", 105), ("m.fa", "m2", 115),
            ],
            dtype=[("spp", object), ("name", object), ("length", int)],
        )
        entries = {
            "h.fa": {"h1": "A" * 100, "h2": "B" * 110},
            "m.fa": {"m1": "C" * 105, "m2": "D" * 115},
        }
        return gl, entries

    def test_single_copy_detected_when_one_per_species(self):
        gl, entries = self._setup_two_species()
        clusters = [["h1", "m1"], ["h2", "m2"]]
        out_clusters, og_cn, ogs_dat, sc_ogs = get_orthogroup_information(
            files=["h.fa", "m.fa"],
            gene_lengths=gl,
            clustering_res=clusters,
            single_copy_threshold=0.5,
            entries=entries,
        )
        # Both clusters have one gene per species → both are single-copy
        assert "OG0" in sc_ogs
        assert "OG1" in sc_ogs
        # OG IDs prepended to clusters
        assert out_clusters[0][0].startswith("OG0")
        # Copy-number header is first
        assert og_cn["files:"] == ["h.fa", "m.fa"]
        # Each multi-species OG has count "1" per species
        assert og_cn["OG0:"] == ["1", "1"]
        assert og_cn["OG1:"] == ["1", "1"]

    def test_multi_copy_excluded_from_single_copy(self):
        gl, entries = self._setup_two_species()
        # OG0 has two human genes — not single-copy
        clusters = [["h1", "h2", "m1"]]
        _, og_cn, _, sc_ogs = get_orthogroup_information(
            files=["h.fa", "m.fa"],
            gene_lengths=gl,
            clustering_res=clusters,
            single_copy_threshold=0.5,
            entries=entries,
        )
        assert sc_ogs == []
        assert og_cn["OG0:"] == ["2", "1"]

    def test_fasta_data_wraps_at_70_chars(self):
        gl, entries = self._setup_two_species()
        clusters = [["h2"]]  # h2 is length 110 → 2 wrap lines
        _, _, ogs_dat, _ = get_orthogroup_information(
            files=["h.fa", "m.fa"],
            gene_lengths=gl,
            clustering_res=clusters,
            single_copy_threshold=0.5,
            entries=entries,
        )
        og_lines = ogs_dat["OG0"]
        # First line = ">h2", then 2 wrapped sequence lines (70 + 40 chars)
        assert og_lines[0] == ">h2"
        assert len(og_lines[1]) == 70
        assert len(og_lines[2]) == 40


# ─── merge_with_gene_lengths + normalize_by_gene_length ─────────────────


class TestMergeAndNormalize:
    def _setup(self):
        gl = np.array(
            [("h.fa", "h1", 100), ("h.fa", "h2", 200)],
            dtype=[("spp", object), ("name", object), ("length", int)],
        )
        # phmmer-style structured array
        dtype_res = [
            ("target_name", "U50"),
            ("query_name",  "U50"),
            ("evalue",      float),
            ("score",       float),
        ]
        res = np.array(
            [("h1", "h2", 1e-50, 200.0)],
            dtype=dtype_res,
        )
        return res, gl

    def test_merge_brings_lengths_alongside(self):
        res, gl = self._setup()
        merged = merge_with_gene_lengths(res, gl)
        # 6 columns: target, query, evalue, score, t_len, q_len
        assert merged.shape == (1, 6)
        assert merged[0, 0] == "h1"
        assert merged[0, 1] == "h2"
        assert merged[0, 3] == 200.0
        # Lengths come from gl
        assert merged[0, 4] == 100  # target h1 length
        assert merged[0, 5] == 200  # query h2 length

    def test_normalize_divides_score_by_sum_of_lengths(self):
        res, gl = self._setup()
        merged = merge_with_gene_lengths(res, gl)
        normalized = normalize_by_gene_length(merged.copy())
        # 200 / (100 + 200) = 0.6667
        assert float(normalized[0, 3]) == pytest.approx(0.6667, rel=1e-3)


# ─── read_and_filter_phmmer_output ──────────────────────────────────────


def _write_phmmer_file(directory, taxon_a, taxon_b, rows):
    """Helper: write a phmmer-formatted tblout file."""
    wd = os.path.join(directory, "orthohmm_working_res")
    os.makedirs(wd, exist_ok=True)
    path = os.path.join(wd, f"{taxon_a}_2_{taxon_b}.phmmerout.txt")
    with open(path, "w") as f:
        f.write("# phmmer tabular output\n")
        for r in rows:
            # phmmer tblout columns: target_name(0) accession(1) query_name(2)
            # accession(3) evalue(4) score(5) ...
            f.write(
                f"{r['target']}\t-\t{r['query']}\t-\t{r['evalue']}\t{r['score']}\n"
            )
        f.write("# [ok]\n")
    return path


class TestReadAndFilterPhmmerOutput:
    def test_reads_phmmer_file_and_filters_by_evalue(self, tmp_path):
        _write_phmmer_file(str(tmp_path), "a.fa", "b.fa", [
            {"target": "t1", "query": "q1", "evalue": "1e-10", "score": "100"},
            {"target": "t2", "query": "q2", "evalue": "1e-3",  "score": "30"},
            {"target": "t3", "query": "q3", "evalue": "1e-50", "score": "200"},
        ])
        res = read_and_filter_phmmer_output(
            "a.fa", "b.fa", str(tmp_path), evalue_threshold=1e-4,
        )
        # 1e-3 fails (not < 1e-4); other two pass
        assert len(res) == 2
        assert sorted(res["target_name"].tolist()) == ["t1", "t3"]

    def test_search_results_path_returns_pair_results(self, tmp_path):
        from orthohmm.search.engine import SpeciesPairResults
        spr = SpeciesPairResults(
            query_species="b.fa", target_species="a.fa",
            target_names=np.array(["t1", "t2"], dtype="U50"),
            query_names=np.array(["q1", "q2"], dtype="U50"),
            evalues=np.array([1e-10, 1e-3]),
            scores=np.array([100.0, 30.0]),
        )
        search_results = {("a.fa", "b.fa"): spr}
        res = read_and_filter_phmmer_output(
            "a.fa", "b.fa", str(tmp_path), evalue_threshold=1e-4,
            search_results=search_results,
        )
        # results_to_phmmer_format filters by threshold → only t1 passes
        assert len(res) == 1
        assert res["target_name"][0] == "t1"

    def test_search_results_missing_pair_returns_empty(self, tmp_path):
        res = read_and_filter_phmmer_output(
            "a.fa", "b.fa", str(tmp_path), evalue_threshold=1e-4,
            search_results={},
        )
        assert len(res) == 0
        assert res.dtype.names == ("target_name", "query_name", "evalue", "score")


# ─── get_best_hits_and_scores ───────────────────────────────────────────


class TestGetBestHitsAndScores:
    def test_highest_scoring_target_kept_per_query(self):
        # 6-column merged: target, query, evalue, score, t_len, q_len
        merged = np.array([
            ["t1", "q1", 1e-10, 100.0, 100, 100],
            ["t2", "q1", 1e-50, 200.0, 100, 100],  # better hit for q1
            ["t3", "q2", 1e-5,   50.0, 100, 100],
        ], dtype=object)
        best = get_best_hits_and_scores(merged)
        assert best["q1"]["target"] == "t2"
        assert best["q1"]["score"] == 200.0
        assert best["q2"]["target"] == "t3"

    def test_empty_input_yields_empty_mapping(self):
        merged = np.empty((0, 6), dtype=object)
        best = get_best_hits_and_scores(merged)
        assert dict(best) == {}


# ─── get_threshold_per_gene ─────────────────────────────────────────────


class TestGetThresholdPerGene:
    def test_rbh_pair_threshold_is_score_average(self):
        a_to_b = {"a1": {"target": "b1", "score": 0.8}}
        b_to_a = {"b1": {"target": "a1", "score": 0.6}}
        score_a = {"a1": 0.8}
        score_b = {"b1": 0.6}
        thresh = get_threshold_per_gene(a_to_b, b_to_a, score_a, score_b, {})
        # Reciprocal pair: threshold = (0.8 + 0.6) / 2 = 0.7
        assert thresh["a1"] == pytest.approx(0.7)

    def test_non_reciprocal_hit_skipped(self):
        # b1's best is b2, not a1 → not reciprocal
        a_to_b = {"a1": {"target": "b1", "score": 0.8}}
        b_to_a = {"b1": {"target": "b2", "score": 0.6}}
        thresh = get_threshold_per_gene(
            a_to_b, b_to_a, {"a1": 0.8}, {"b1": 0.6}, {},
        )
        assert thresh == {}

    def test_minimum_score_wins_across_pairs(self):
        a_to_b = {"a1": {"target": "b1", "score": 1.0}}
        b_to_a = {"b1": {"target": "a1", "score": 1.0}}
        # Pre-existing threshold higher than current → keep current
        thresh = get_threshold_per_gene(
            a_to_b, b_to_a, {"a1": 0.5}, {"b1": 0.5},
            reciprocal_best_hit_thresholds={"a1": 0.9},
        )
        assert thresh["a1"] == 0.5
        # Pre-existing threshold lower than current → keep existing
        thresh2 = get_threshold_per_gene(
            a_to_b, b_to_a, {"a1": 0.5}, {"b1": 0.5},
            reciprocal_best_hit_thresholds={"a1": 0.1},
        )
        assert thresh2["a1"] == 0.1


# ─── correct_by_phylogenetic_distance ───────────────────────────────────


class TestCorrectByPhylogeneticDistance:
    def test_score_divided_by_mean_rbh_correction(self):
        a_to_b = {
            "a1": {"target": "b1", "score": 1.0},
            "a2": {"target": "b2", "score": 2.0},
        }
        b_to_a = {
            "b1": {"target": "a1", "score": 1.0},
            "b2": {"target": "a2", "score": 2.0},
        }
        # RBH scores = [(1+1)/2, (2+2)/2] = [1, 2] → mean = 1.5
        sc_ab, sc_ba, corr = correct_by_phylogenetic_distance(
            a_to_b, b_to_a, ("a.fa", "b.fa"), {},
        )
        assert corr[frozenset(("a.fa", "b.fa"))] == pytest.approx(1.5)
        assert sc_ab["a1"] == pytest.approx(1.0 / 1.5)
        assert sc_ba["b2"] == pytest.approx(2.0 / 1.5)

    def test_existing_corr_averaged_with_new(self):
        a_to_b = {"a1": {"target": "b1", "score": 2.0}}
        b_to_a = {"b1": {"target": "a1", "score": 2.0}}
        pair = ("a.fa", "b.fa")
        pre = {frozenset(pair): 1.0}  # prior correction = 1.0
        _, _, corr = correct_by_phylogenetic_distance(
            a_to_b, b_to_a, pair, pre,
        )
        # New mean = 2.0; averaged with prior 1.0 → 1.5
        assert corr[frozenset(pair)] == pytest.approx(1.5)


# ─── update_progress ────────────────────────────────────────────────────


class TestUpdateProgress:
    def test_increments_counter_and_writes_percent(self, capsys):
        completed = multiprocessing.Value("i", 0)
        lock = multiprocessing.Lock()
        update_progress(lock, completed, total_tasks=4)
        update_progress(lock, completed, total_tasks=4)
        assert completed.value == 2
        out = capsys.readouterr().out
        # Latest write is 2/4 = 50.00%
        assert "50.00%" in out


# ─── process_pair_edge_thresholds (file-based path) ─────────────────────


class TestProcessPairEdgeThresholds:
    def test_reciprocal_pair_produces_threshold(self, tmp_path):
        # Two species, one RBH pair → threshold should be populated
        _write_phmmer_file(str(tmp_path), "a.fa", "b.fa", [
            {"target": "b1", "query": "a1", "evalue": "1e-50", "score": "100"},
        ])
        _write_phmmer_file(str(tmp_path), "b.fa", "a.fa", [
            {"target": "a1", "query": "b1", "evalue": "1e-50", "score": "100"},
        ])
        gl = np.array(
            [("a.fa", "a1", 100), ("b.fa", "b1", 100)],
            dtype=[("spp", object), ("name", object), ("length", int)],
        )
        thresholds, corr = process_pair_edge_thresholds(
            ("a.fa", "b.fa"), str(tmp_path), gl, evalue_threshold=1e-4,
        )
        # Normalized score = 100/(100+100) = 0.5; correction = 0.5;
        # corrected score = 1.0; threshold (avg of both directions) = 1.0
        assert "a1" in thresholds
        assert corr[frozenset(("a.fa", "b.fa"))] == pytest.approx(0.5)


def test_best_indexed_hits_keeps_first_equal_scoring_target():
    result = IndexedSpeciesPairResults(
        query_species="a.fa",
        target_species="b.fa",
        query_indices=np.array([0, 0, 1], dtype=np.int32),
        target_indices=np.array([1, 0, 0], dtype=np.int32),
        evalues=np.array([1e-20, 1e-30, 1.0]),
        scores=np.array([2.0, 2.0, 9.0]),
    )
    targets, scores = _best_indexed_hits(result, 2, 1e-4)
    assert targets.tolist() == [1, -1]
    assert scores[0] == 2.0
    assert np.isnan(scores[1])


def test_indexed_edge_thresholds_match_legacy_rbh_semantics(tmp_path):
    (tmp_path / "a.fa").write_text(">a0\nAAAA\n>a1\nAAAA\n")
    (tmp_path / "b.fa").write_text(">b0\nAAAA\n>b1\nAAAA\n")
    ab = IndexedSpeciesPairResults(
        query_species="a.fa",
        target_species="b.fa",
        query_indices=np.array([0, 0, 1], dtype=np.int32),
        target_indices=np.array([0, 1, 1], dtype=np.int32),
        evalues=np.array([1e-20, 1e-20, 1e-20]),
        scores=np.array([2.0, 2.0, 4.0]),
    )
    ba = IndexedSpeciesPairResults(
        query_species="b.fa",
        target_species="a.fa",
        query_indices=np.array([0, 1], dtype=np.int32),
        target_indices=np.array([0, 1], dtype=np.int32),
        evalues=np.array([1e-20, 1e-20]),
        scores=np.array([2.0, 4.0]),
    )
    search_results = IndexedSearchResults(
        {("a.fa", "b.fa"): ab, ("b.fa", "a.fa"): ba},
        {"a.fa": ["a0", "a1"], "b.fa": ["b0", "b1"]},
    )

    _, thresholds, corr = determine_edge_thresholds(
        ["a.fa", "b.fa"], str(tmp_path), str(tmp_path), 1, 1e-4,
        search_results=search_results,
    )

    assert isinstance(thresholds, IndexedEdgeThresholds)
    assert thresholds.values == pytest.approx([2 / 3, 4 / 3, 2 / 3, 4 / 3])
    assert corr[frozenset(("a.fa", "b.fa"))] == pytest.approx(3.0)
    assert corr[frozenset(("a.fa",))] == 0.0
    assert corr[frozenset(("b.fa",))] == 0.0


def test_indexed_thresholds_preserve_legacy_float_operation_order(tmp_path):
    for species in ("a", "b"):
        (tmp_path / f"{species}.fa").write_text(
            "".join(f">{species}{idx}\nAAAA\n" for idx in range(3))
        )
    fwd_scores = np.array([
        9.430561055723675,
        5.113275528143616,
        9.762437057077042,
    ])
    rev_scores = np.array([
        0.8083602389560218,
        6.073558319950296,
        3.764865843772726,
    ])

    def pair(query, target, scores):
        indices = np.arange(3, dtype=np.int32)
        return IndexedSpeciesPairResults(
            query_species=query,
            target_species=target,
            query_indices=indices,
            target_indices=indices,
            evalues=np.full(3, 1e-20),
            scores=scores,
        )

    search_results = IndexedSearchResults(
        {
            ("a.fa", "b.fa"): pair("a.fa", "b.fa", fwd_scores),
            ("b.fa", "a.fa"): pair("b.fa", "a.fa", rev_scores),
        },
        {"a.fa": ["a0", "a1", "a2"], "b.fa": ["b0", "b1", "b2"]},
    )
    _, thresholds, corr = determine_edge_thresholds(
        ["a.fa", "b.fa"], str(tmp_path), str(tmp_path), 1, 1e-4,
        search_results=search_results,
    )
    correction = corr[frozenset(("a.fa", "b.fa"))]
    legacy_order = (fwd_scores / correction + rev_scores / correction) / 2.0

    assert np.array_equal(thresholds.values[:3], legacy_order)


# ─── process_pair_determine_network_edges ───────────────────────────────


class TestProcessPairDetermineNetworkEdges:
    def test_hit_above_threshold_yields_edge(self, tmp_path):
        _write_phmmer_file(str(tmp_path), "a.fa", "b.fa", [
            {"target": "b1", "query": "a1", "evalue": "1e-50", "score": "100"},
        ])
        # gene_lengths comes through this function as a {name: length} dict
        gene_lengths = {"a1": 100, "b1": 100}
        pairwise_rbh_corr = {frozenset(("a.fa", "b.fa")): 0.5}
        thresholds = {"a1": 0.5}
        edges = process_pair_determine_network_edges(
            ("a.fa", "b.fa"), str(tmp_path), gene_lengths,
            pairwise_rbh_corr, thresholds, evalue_threshold=1e-4,
        )
        # 100 / (100+100) / 0.5 = 1.0, which is >= threshold 0.5
        assert frozenset(("a1", "b1")) in edges
        assert edges[frozenset(("a1", "b1"))] == pytest.approx(1.0)

    def test_self_hit_excluded(self, tmp_path):
        # query == target → frozenset has length 1, should be skipped
        _write_phmmer_file(str(tmp_path), "a.fa", "a.fa", [
            {"target": "a1", "query": "a1", "evalue": "1e-50", "score": "100"},
        ])
        gene_lengths = {"a1": 100}
        edges = process_pair_determine_network_edges(
            ("a.fa", "a.fa"), str(tmp_path), gene_lengths,
            {frozenset(("a.fa", "a.fa")): 0.5}, {"a1": 0.0},
            evalue_threshold=1e-4,
        )
        assert edges == {}


def test_indexed_network_edges_match_legacy_threshold_semantics(tmp_path):
    (tmp_path / "orthohmm_working_res").mkdir()
    pair = IndexedSpeciesPairResults(
        query_species="a.fa",
        target_species="b.fa",
        query_indices=np.array([0, 1], dtype=np.int32),
        target_indices=np.array([0, 0], dtype=np.int32),
        evalues=np.array([1e-20, 1e-20]),
        scores=np.array([2.0, 0.5]),
        candidate_count=2,
    )
    search_results = IndexedSearchResults(
        {("a.fa", "b.fa"): pair},
        {"a.fa": ["a0", "a1"], "b.fa": ["b0"]},
    )
    gene_lengths = np.array(
        [("a.fa", "a0", 10), ("a.fa", "a1", 10), ("b.fa", "b0", 10)],
        dtype=[("spp", object), ("name", object), ("length", int)],
    )
    edges = _determine_indexed_network_edges(
        search_results,
        gene_lengths,
        {frozenset(("a.fa", "b.fa")): 2.0},
        {"a0": 0.75, "a1": 0.75},
        1e-4,
        str(tmp_path),
    )
    assert isinstance(edges, IndexedEdges)
    assert edges.sources.tolist() == [0]
    assert edges.targets.tolist() == [2]
    assert edges.weights.tolist() == [1.0]
    assert (tmp_path / "orthohmm_working_res" / "orthohmm_edges.txt").read_text() == (
        "a0\tb0\t1.0\n"
    )


# ─── generate_orthogroup_clusters_file ──────────────────────────────────


def _make_minimal_outdir(tmp_path, edges_clustered_lines):
    wd = tmp_path / "orthohmm_working_res"
    wd.mkdir(parents=True, exist_ok=True)
    (wd / "orthohmm_edges_clustered.txt").write_text(
        "\n".join(edges_clustered_lines) + "\n"
    )
    return wd


class TestGenerateOrthogroupClustersFile:
    def test_reads_clusters_and_assigns_og_ids(self, tmp_path):
        _make_minimal_outdir(tmp_path, ["a1 b1", "a2 b2"])
        # Two species with one gene each per cluster
        fasta_dir = tmp_path / "fasta"
        fasta_dir.mkdir()
        (fasta_dir / "a.fa").write_text(">a1\nAAA\n>a2\nAAA\n")
        (fasta_dir / "b.fa").write_text(">b1\nBBB\n>b2\nBBB\n")
        gl = np.array(
            [("a.fa", "a1", 3), ("a.fa", "a2", 3),
             ("b.fa", "b1", 3), ("b.fa", "b2", 3)],
            dtype=[("spp", object), ("name", object), ("length", int)],
        )
        singletons, og_cn, ogs_dat, sc_ogs = generate_orthogroup_clusters_file(
            output_directory=str(tmp_path),
            gene_lengths=gl,
            files=["a.fa", "b.fa"],
            single_copy_threshold=0.5,
            fasta_directory=str(fasta_dir),
        )
        # No singletons — every gene is in a cluster
        assert singletons == []
        # Two single-copy OGs (one gene per species per cluster)
        assert sorted(sc_ogs) == ["OG0", "OG1"]
        # Clusters file written at output_directory root
        assert (tmp_path / "orthohmm_orthogroups.txt").exists()


# ─── generate_orthogroup_files ──────────────────────────────────────────


class TestGenerateOrthogroupFiles:
    def test_writes_expected_files(self, tmp_path):
        # Layout the working subdirectory the writer functions expect
        (tmp_path / "orthohmm_working_res").mkdir(parents=True, exist_ok=True)
        (tmp_path / "orthohmm_orthogroups").mkdir(parents=True, exist_ok=True)
        (tmp_path / "orthohmm_single_copy_orthogroups").mkdir(parents=True, exist_ok=True)

        gl = np.array(
            [("a.fa", "a1", 3), ("b.fa", "b1", 3)],
            dtype=[("spp", object), ("name", object), ("length", int)],
        )
        og_cn = {"files:": ["a.fa", "b.fa"], "OG0:": ["1", "1"]}
        ogs_dat = {"OG0": [">a1", "AAA", ">b1", "BBB"]}
        single_copy_ogs = ["OG0"]
        generate_orthogroup_files(
            output_directory=str(tmp_path),
            gene_lengths=gl,
            og_cn=og_cn,
            ogs_dat=ogs_dat,
            single_copy_ogs=single_copy_ogs,
        )
        # Copy-number and single-copy-name files written
        assert (tmp_path / "orthohmm_gene_count.txt").exists()
        assert (tmp_path / "orthohmm_single_copy_orthogroups.txt").exists()
