"""Unit tests for orthohmm/helpers.py (pure-function utilities)."""
import numpy as np
import pytest

from orthohmm.helpers import (
    StartStep,
    StopStep,
    SubstitutionMatrix,
    generate_phmmer_cmds,
    get_all_fasta_entries,
    get_orthogroup_information,
    get_sequence_lengths,
    get_singletons,
    merge_with_gene_lengths,
    normalize_by_gene_length,
)


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
        # 3 singletons appended (g3, g4, g5), order not guaranteed
        singleton_genes = sorted(s[0] for s in singletons)
        assert singleton_genes == ["g3", "g4", "g5"]
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
