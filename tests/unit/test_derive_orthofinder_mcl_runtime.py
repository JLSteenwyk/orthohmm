import pytest

from benchmark_tools.derive_orthofinder_mcl_runtime import derive_mcl_runtime


def test_derives_runtime_when_markers_are_in_separate_logs(tmp_path):
    metadata_log = tmp_path / "Log.txt"
    metadata_log.write_text(
        "2026-09-07 03:02:39 : Started OrthoFinder version 3.1.5\n"
        "2026-09-07 04:48:09 : OrthoFinder run completed\n"
    )
    stdout_log = tmp_path / "run.log"
    stdout_log.write_text("2026-09-07 04:24:03 : Ran MCL\n")

    assert derive_mcl_runtime([metadata_log, stdout_log]) == 4884


def test_rejects_missing_mcl_marker(tmp_path):
    log = tmp_path / "Log.txt"
    log.write_text(
        "2026-09-07 03:02:39 : Started OrthoFinder version 3.1.5\n"
    )

    with pytest.raises(ValueError, match="Could not derive"):
        derive_mcl_runtime([log])
