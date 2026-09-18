import json
import os
from pathlib import Path
import shutil
import subprocess

import pytest

from benchmark_tools.probe_orthomcl_bpo_parity import expected_indexes, compare

PERL = shutil.which("perl")
SCRIPT = Path(__file__).resolve().parents[2] / "benchmark_tools/validate_orthomcl_bpo_indexes.pl"


def fixture(tmp_path):
    bpo = tmp_path / "all.bpo"
    bpo.write_text("1;A;10;A;10;0;100;1:1-10:1-10.\n"
                   "2;A;10;B;10;1e-5;80;1:1-10:1-10.\n"
                   "3;B;10;B;10;0;100;1:1-10:1-10.\n")
    return bpo


def test_expected_offsets_and_inclusive_ranges(tmp_path):
    bpo = fixture(tmp_path)
    offsets, ranges = expected_indexes(bpo)
    assert offsets[0] == 0 and offsets[-1] == bpo.stat().st_size
    assert len(offsets) == 4 and ranges == {"A": "1;2", "B": "3;3"}
    assert compare(bpo, bpo, {"offsets": offsets, "query_ranges": ranges})["bpo_byte_identical"] is True


@pytest.mark.parametrize("text", ["", "1;A\n", "2;A;10;A;10;0;100;x\n",
    "1;;10;A;10;0;100;x\n", "1;A;10;A;10;0;100;x\n2;B;10;B;10;0;100;x\n3;A;10;A;10;0;100;x\n"])
def test_invalid_bpo_rejected(tmp_path, text):
    path = tmp_path / "bad.bpo"
    path.write_text(text)
    with pytest.raises(ValueError):
        expected_indexes(path)


def test_byte_difference_and_index_difference_are_distinct(tmp_path):
    bpo = fixture(tmp_path)
    other = tmp_path / "other.bpo"
    other.write_text(bpo.read_text().replace(";80;", ";81;"))
    offsets, ranges = expected_indexes(bpo)
    assert compare(bpo, other, {"offsets": offsets, "query_ranges": ranges})["bpo_byte_identical"] is False
    offsets[-1] += 1
    with pytest.raises(ValueError, match="index differs"):
        compare(bpo, other, {"offsets": offsets, "query_ranges": ranges})


@pytest.mark.skipif(os.environ.get("ORTHOHMM_LEGACY_BLAST_SMOKE") != "1",
                    reason="Opt-in installed OrthoMCL/BioPerl parity probe")
def test_installed_native_bpo_parity(tmp_path):
    from benchmark_tools.probe_orthomcl_bpo_parity import run
    result = run(tmp_path / "probe")
    assert result["status"] == "native_bpo_fixture_parity_verified"
    assert result["content"]["native_pair_records"] == result["content"]["streaming_pair_records"] == 6
    assert result["streaming_index_validation"]["records"] == 6
    assert result["content"]["indexes_verified"] is True
    assert result["accuracy_admitted"] is False


@pytest.mark.skipif(PERL is None, reason="Perl core modules required for native index tests")
@pytest.mark.parametrize("problem", [None, "offset", "negative", "eof", "extra_offset", "short_offsets",
    "start", "end", "missing_query", "extra_query", "index_type", "ranges_type", "record_id", "repeated_query"])
def test_streaming_native_index_validation(tmp_path, problem):
    bpo = fixture(tmp_path)
    offsets, ranges = expected_indexes(bpo)
    if problem == "offset":
        offsets[1] += 1
    elif problem == "negative":
        offsets[0] = -1
    elif problem == "eof":
        offsets[-1] += 1
    elif problem == "extra_offset":
        offsets.append(offsets[-1])
    elif problem == "short_offsets":
        offsets.pop()
    elif problem == "start":
        ranges["A"] = "2;2"
    elif problem == "end":
        ranges["A"] = "1;3"
    elif problem == "missing_query":
        del ranges["B"]
    elif problem == "extra_query":
        ranges["C"] = "4;4"
    elif problem == "index_type":
        offsets = {}
    elif problem == "ranges_type":
        ranges = []
    elif problem == "record_id":
        bpo.write_text(bpo.read_text().replace("3;B;", "4;B;"))
    elif problem == "repeated_query":
        bpo.write_text("1;A;10;A;10;0;100;x\n2;B;10;B;10;0;100;x\n3;A;10;A;10;0;100;x\n")
        lines = bpo.read_bytes().splitlines(keepends=True)
        offsets = [0, len(lines[0]), len(lines[0]) + len(lines[1]), bpo.stat().st_size]
        ranges = {"A": "1;1", "B": "2;2"}
    idx, se = tmp_path / "all.idx", tmp_path / "all.se"
    subprocess.run([PERL, "-MStorable=store", "-MJSON::PP", "-e",
        "local $/; my $v=decode_json(<STDIN>); store($v->[0],$ARGV[0]); store($v->[1],$ARGV[1]);",
        str(idx), str(se)], input=json.dumps([offsets, ranges]), text=True, check=True, capture_output=True)
    done = subprocess.run([PERL, str(SCRIPT), str(bpo), str(idx), str(se)], text=True, capture_output=True)
    if problem:
        assert done.returncode != 0 and done.stderr
        assert not done.stdout
    else:
        assert done.returncode == 0, done.stderr
        report = json.loads(done.stdout)
        assert report["records"] == 3 and report["queries"] == 2
        assert report["offset_entries_including_eof"] == 4
