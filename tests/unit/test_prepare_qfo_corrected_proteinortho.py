import pytest

from benchmark_tools.prepare_qfo_corrected_proteinortho import native_command


def names():
    return [f"species_{i:02d}.fasta" for i in range(78)]


def test_preserves_native_defaults_and_frozen_resources():
    command = native_command("/runtime", "/image.sif", reversed(names()))
    assert command[:8] == ["/runtime", "exec", "--bind", "/mnt", "/image.sif",
                           "proteinortho", "-project=qfo", "-cpus=32"]
    assert command[8:] == names()
    assert not any("step" in arg or "resume" in arg for arg in command)


@pytest.mark.parametrize("replacement", ["../foreign.fasta", "/absolute.fasta", "-flag.fasta", "species.txt"])
def test_rejects_noncanonical_input(replacement):
    files = names()
    files[0] = replacement
    with pytest.raises(ValueError, match="basenames"):
        native_command("/runtime", "/image", files)


@pytest.mark.parametrize("files", [[], names()[:-1], names() + ["extra.fasta"], names()[:-1] + [names()[0]]])
def test_requires_complete_unique_inventory(files):
    with pytest.raises(ValueError, match="78 distinct"):
        native_command("/runtime", "/image", files)
