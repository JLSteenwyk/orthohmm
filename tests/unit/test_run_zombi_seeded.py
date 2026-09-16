import pytest

from benchmark_tools.run_zombi_seeded import family_lengths, family_seed, length_aware_simulator, seeded_evolver


def test_family_seed_independent_of_output_directory():
    assert family_seed(17, "/first/1_complete.fasta") == family_seed(17, "/second/1_complete.fasta")
    assert family_seed(17, "1_complete.fasta") != family_seed(18, "1_complete.fasta")
    assert family_seed(17, "1_complete.fasta") != family_seed(17, "2_complete.fasta")


@pytest.mark.parametrize("seed,path", [(0, "x"), (-1, "x"), (1, None)])
def test_invalid_seed_inputs(seed, path):
    with pytest.raises(ValueError):
        family_seed(seed, path)


def test_adapter_passes_explicit_seed_without_changing_other_arguments():
    class Original:
        def __call__(self, **kwargs):
            return kwargs
    adapter = seeded_evolver(Original, 12)()
    result = adapter(seqfile="a.fasta", ratefile=None, write_anc=True)
    assert result == dict(seqfile="a.fasta", ratefile=None, write_anc=True, seed=family_seed(12, "a.fasta"))
    assert adapter(seqfile="a.fasta", seed=9)["seed"] == 9


def test_length_rule_is_stable_bounded_and_order_independent():
    import hashlib
    families = [str(i) for i in range(1, 101)]
    result = family_lengths(20261101, families)
    assert result == family_lengths(20261101, reversed(families))
    assert result != family_lengths(20261102, families)
    assert min(result.values()) >= 100 and max(result.values()) <= 500
    assert len(set(result.values())) > 1
    assert result["1"] == 100 + int(hashlib.sha256(b"orthohmm-sim-length-v2:20261101:1").hexdigest()[:16], 16) % 401


@pytest.mark.parametrize("seed,families", [(0, ["1"]), (True, ["1"]), (1, []), (1, ["01"]), (1, ["1", "1"]), (1, ["0"])])
def test_invalid_length_mapping_inputs(seed, families):
    with pytest.raises(ValueError):
        family_lengths(seed, families)


def test_length_adapter_only_changes_size_and_restores_after_error():
    class Original:
        size = 300
        def run(self, tree, output):
            if output == "error":
                raise RuntimeError("native failure")
            return tree, output, self.size
    adapter = length_aware_simulator(Original, {"1": 123})()
    assert adapter.run("/a/1_completetree.nwk", "/out") == ("/a/1_completetree.nwk", "/out", 123)
    assert adapter.size == 300
    with pytest.raises(RuntimeError, match="native failure"):
        adapter.run("1_completetree.nwk", "error")
    assert adapter.size == 300
    with pytest.raises(ValueError, match="no frozen"):
        adapter.run("2_completetree.nwk", "/out")
