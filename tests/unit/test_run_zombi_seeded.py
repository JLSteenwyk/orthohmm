import pytest

from benchmark_tools.run_zombi_seeded import family_seed, seeded_evolver


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
