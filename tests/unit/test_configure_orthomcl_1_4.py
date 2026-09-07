import pytest

from benchmark_tools.configure_orthomcl_1_4 import configure_module


def test_configures_isolated_paths_and_threads(tmp_path):
    module = tmp_path / "orthomcl_module.pm"
    module.write_text(
        'our $BLAST_NOCPU = 8; # worker count\n'
        'our $PATH_TO_ORTHOMCL = "/old/tool/"; # path\n'
        'our $ORTHOMCL_DATA_DIR = $PATH_TO_ORTHOMCL."/sample_data/";\n'
        "our $UNCHANGED = 1;\n"
    )

    configure_module(module, tmp_path / "tool", tmp_path / "tool/data", 32)

    assert module.read_text() == (
        "our $BLAST_NOCPU = 32;\n"
        f'our $PATH_TO_ORTHOMCL = "{tmp_path}/tool/";\n'
        f'our $ORTHOMCL_DATA_DIR = "{tmp_path}/tool/data/";\n'
        "our $UNCHANGED = 1;\n"
    )


def test_rejects_missing_assignment(tmp_path):
    module = tmp_path / "orthomcl_module.pm"
    module.write_text("our $BLAST_NOCPU = 8;\n")

    with pytest.raises(ValueError, match=r"\$PATH_TO_ORTHOMCL assignment"):
        configure_module(module, tmp_path / "tool", tmp_path / "data", 32)
