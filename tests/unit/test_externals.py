"""Unit tests for orthohmm/externals.py (clustering + run-completion helpers)."""
import os
import signal
import subprocess

import numpy as np

from orthohmm import externals
from orthohmm.externals import (
    check_if_mcl_command_completed,
    check_if_phmmer_command_completed,
    execute_leiden,
)
from orthohmm.helpers import IndexedEdges


# ─── completion-check helpers ───────────────────────────────────────────


class TestCheckPhmmerCompleted:
    def test_nonexistent_file_returns_false(self, tmp_path):
        assert check_if_phmmer_command_completed(str(tmp_path / "nope")) is False

    def test_incomplete_file_returns_false(self, tmp_path):
        p = tmp_path / "partial.txt"
        p.write_text("# header\nsome data\n")
        assert check_if_phmmer_command_completed(str(p)) is False

    def test_well_formed_phmmer_output_returns_true(self, tmp_path):
        # The completion sentinel phmmer writes is "# [ok]" on the last line.
        p = tmp_path / "done.txt"
        p.write_text("# header\nrow1\nrow2\n# [ok]\n")
        assert check_if_phmmer_command_completed(str(p)) is True


class TestCheckMclCompleted:
    def test_nonexistent_file_returns_false(self, tmp_path):
        assert check_if_mcl_command_completed(str(tmp_path / "nope")) is False

    def test_complete_mcl_output_returns_true(self, tmp_path):
        # MCL signals completion with its citation line on the last line.
        p = tmp_path / "done.txt"
        p.write_text(
            "cluster1 g1 g2\n"
            "cluster2 g3 g4\n"
            "    ( http://link.aip.org/link/?SJMAEL/30/121/1 )\n"
        )
        assert check_if_mcl_command_completed(str(p)) is True


# ─── execute_leiden (in-process Leiden CPM) ─────────────────────────────


def _write_edges_abc(directory, edges):
    """Helper: write an orthohmm_edges.txt under the expected working_res path."""
    wd = os.path.join(directory, "orthohmm_working_res")
    os.makedirs(wd, exist_ok=True)
    with open(os.path.join(wd, "orthohmm_edges.txt"), "w") as f:
        for a, b, w in edges:
            f.write(f"{a}\t{b}\t{w}\n")
    return wd


class TestExecuteLeiden:
    def test_two_disjoint_components_become_two_clusters(self, tmp_path):
        edges = [
            ("a1", "a2", 1.0),
            ("a2", "a3", 1.0),
            # disjoint from above:
            ("b1", "b2", 1.0),
            ("b2", "b3", 1.0),
        ]
        wd = _write_edges_abc(str(tmp_path), edges)
        execute_leiden(0.1, str(tmp_path))
        out = os.path.join(wd, "orthohmm_edges_clustered.txt")
        assert os.path.exists(out)
        clusters = [set(line.split()) for line in open(out) if line.strip()]
        # Should produce exactly two clusters; member sets correspond to the
        # two disjoint components
        assert {frozenset(c) for c in clusters} == {
            frozenset(["a1", "a2", "a3"]),
            frozenset(["b1", "b2", "b3"]),
        }

    def test_auto_resolution_picks_4x_min_weight(self, tmp_path, capsys):
        """Auto formula is γ = 4 × min(positive_weight)."""
        edges = [
            ("a", "b", 0.5),
            ("b", "c", 1.0),
            ("a", "c", 0.25),  # min = 0.25 → γ_auto = 1.0
        ]
        _write_edges_abc(str(tmp_path), edges)
        execute_leiden("auto", str(tmp_path))
        log = capsys.readouterr().out
        # γ printed as "γ=<value>" using 6-significant-digit %g formatting
        assert "γ=1" in log
        out = os.path.join(tmp_path, "orthohmm_working_res",
                           "orthohmm_edges_clustered.txt")
        assert os.path.exists(out)

    def test_output_format_matches_mcl_one_cluster_per_line(self, tmp_path):
        edges = [("g1", "g2", 1.0)]
        _write_edges_abc(str(tmp_path), edges)
        execute_leiden(0.1, str(tmp_path))
        content = open(
            os.path.join(tmp_path, "orthohmm_working_res",
                         "orthohmm_edges_clustered.txt")
        ).read()
        # One cluster on a single whitespace-separated line, no leading tag,
        # alphabetically sorted (matches what helpers.generate_orthogroup_
        # clusters_file consumes via .split())
        assert content.strip() == "g1 g2"

    def test_compact_edges_bypass_abc_reload(self, tmp_path):
        wd = os.path.join(tmp_path, "orthohmm_working_res")
        os.makedirs(wd)
        edges = IndexedEdges(
            gene_names=["a", "b", "c", "d"],
            sources=np.array([0, 2], dtype=np.int32),
            targets=np.array([1, 3], dtype=np.int32),
            weights=np.array([1.0, 1.0]),
        )
        execute_leiden(0.1, str(tmp_path), edges=edges)
        clusters = [
            line.split()
            for line in open(os.path.join(wd, "orthohmm_edges_clustered.txt"))
            if line.strip()
        ]
        assert clusters == [["a", "b"], ["c", "d"]]

    def test_compact_edges_can_include_isolated_genes(self, tmp_path):
        wd = os.path.join(tmp_path, "orthohmm_working_res")
        os.makedirs(wd)
        edges = IndexedEdges(
            gene_names=["a", "b", "isolated"],
            sources=np.array([0], dtype=np.int32),
            targets=np.array([1], dtype=np.int32),
            weights=np.array([1.0]),
        )

        execute_leiden(
            0.1,
            str(tmp_path),
            edges=edges,
            include_isolates=True,
        )

        clusters = [
            line.split()
            for line in open(os.path.join(wd, "orthohmm_edges_clustered.txt"))
            if line.strip()
        ]
        assert clusters == [["a", "b"], ["isolated"]]

    def test_large_compact_edges_use_isolated_worker(self, tmp_path, monkeypatch):
        wd = os.path.join(tmp_path, "orthohmm_working_res")
        os.makedirs(wd)
        edges = IndexedEdges(
            gene_names=["a", "b", "isolated"],
            sources=np.array([0], dtype=np.int32),
            targets=np.array([1], dtype=np.int32),
            weights=np.array([1.0]),
        )
        monkeypatch.setattr(externals, "LEIDEN_ISOLATION_MIN_EDGES", 1)

        execute_leiden(
            0.1,
            str(tmp_path),
            edges=edges,
            include_isolates=True,
            seed=4,
        )

        clusters = [
            line.split()
            for line in open(os.path.join(wd, "orthohmm_edges_clustered.txt"))
            if line.strip()
        ]
        assert clusters == [["a", "b"], ["isolated"]]
        assert not any(
            name.startswith(".leiden-payload-") for name in os.listdir(wd)
        )

    def test_isolated_worker_retries_sigsegv_and_removes_stale_output(
        self,
        tmp_path,
        monkeypatch,
    ):
        wd = os.path.join(tmp_path, "orthohmm_working_res")
        os.makedirs(wd)
        output = os.path.join(wd, "orthohmm_edges_clustered.txt")
        with open(output, "w") as handle:
            handle.write("stale result\n")
        edges = IndexedEdges(
            gene_names=["a", "b"],
            sources=np.array([0], dtype=np.int32),
            targets=np.array([1], dtype=np.int32),
            weights=np.array([1.0]),
        )
        monkeypatch.setattr(externals, "LEIDEN_ISOLATION_MIN_EDGES", 1)
        real_run = subprocess.run
        worker_attempts = 0

        def flaky_run(command, **kwargs):
            nonlocal worker_attempts
            if command[1:3] == ["-m", "orthohmm.leiden_worker"]:
                worker_attempts += 1
                assert not os.path.exists(output)
                if worker_attempts == 1:
                    raise subprocess.CalledProcessError(-signal.SIGSEGV, command)
            return real_run(command, **kwargs)

        monkeypatch.setattr(externals.subprocess, "run", flaky_run)

        execute_leiden(0.1, str(tmp_path), edges=edges, seed=4)

        assert worker_attempts == 2
        assert open(output).read().strip() == "a b"
