import pytest

from benchmark_tools import shared_qfo_export as export
from benchmark_tools.admit_qfo_factorial_cell import check_pairs
from benchmark_tools.prepare_ob_candidate_neighborhood import record

DUPLICATE = ("a", "A", "b", "B")
PAIR = b"a\tA\tb\tB\n"
CONFIDENCE = PAIR[:-1] + b"\tmedium\n"
OTHER = b"a\tA\tc\tC\n"
OTHER_CONFIDENCE = OTHER[:-1] + b"\thigh\n"


@pytest.fixture
def files(tmp_path):
    pairs, confidence, output_pairs, output_confidence = [tmp_path / name
        for name in ("pairs", "confidence", "derived_pairs", "derived_confidence")]
    pairs.write_bytes(export.HEADER + PAIR * 2 + OTHER)
    confidence.write_bytes(export.CONFIDENCE_HEADER + CONFIDENCE * 2 + OTHER_CONFIDENCE)
    return pairs, confidence, output_pairs, output_confidence


def normalize(files, **kwargs):
    return export.normalize(*files, kwargs.get("duplicate", DUPLICATE),
                            kwargs.get("rows", 3), kwargs.get("unique", 2))


def test_exact_duplicate_only_derived_files_change(files):
    before = [record(path) for path in files[:2]]
    proof = normalize(files)
    assert [record(path) for path in files[:2]] == before
    assert files[2].read_bytes() == export.HEADER + PAIR + OTHER
    assert files[3].read_bytes() == export.CONFIDENCE_HEADER + CONFIDENCE + OTHER_CONFIDENCE
    assert proof["original_rows"] == 3 and proof["distinct_rows"] == 2
    assert proof["skips"]["pairs"]["line"] == 3
    assert export.verify_projection(proof)["original_files_unchanged"]
    assert check_pairs(files[2], {"a": "A", "b": "B", "c": "C"}, {"a": 0, "b": 0, "c": 0}, 2) == 2


@pytest.mark.parametrize("body,annotated", [
    (PAIR + OTHER, CONFIDENCE + OTHER_CONFIDENCE),
    (PAIR * 3 + OTHER, CONFIDENCE * 3 + OTHER_CONFIDENCE),
    (OTHER + PAIR * 2, OTHER_CONFIDENCE + CONFIDENCE * 2),
    (PAIR * 2, CONFIDENCE + PAIR[:-1] + b"\thigh\n"),
    (PAIR * 2, CONFIDENCE),
    (PAIR, CONFIDENCE * 2),
    (b"b\tB\ta\tA\n" * 2, b"b\tB\ta\tA\thigh\n" * 2),
    (b"a\tA\tb\tA\n" * 2, b"a\tA\tb\tA\thigh\n" * 2),
    (PAIR, b"a\tA\tb\tB\tunknown\n"),
    (PAIR.rstrip(b"\n"), CONFIDENCE.rstrip(b"\n")),
    (b"a\tA\tb\n", b"a\tA\tb\thigh\n"),
    (PAIR, b"a\tA\tc\tB\thigh\n"),
    (b"a \tA\tb\tB\n", b"a \tA\tb\tB\thigh\n"),
])
def test_any_other_export_anomaly_refused(files, body, annotated):
    files[0].write_bytes(export.HEADER + body)
    files[1].write_bytes(export.CONFIDENCE_HEADER + annotated)
    with pytest.raises(ValueError):
        normalize(files)


@pytest.mark.parametrize("change", [{"duplicate": ("a", "A", "c", "C")}, {"rows": 4, "unique": 3},
    {"rows": 3, "unique": 1}, {"rows": True}, {"duplicate": ("a",)}])
def test_unplanned_duplicate_or_counts_refused(files, change):
    with pytest.raises(ValueError):
        normalize(files, **change)


def test_bad_header_refused(files):
    files[0].write_bytes(b"wrong header\n" + PAIR * 2 + OTHER)
    with pytest.raises(ValueError, match="headers"):
        normalize(files)


def test_existing_output_not_overwritten(files):
    files[2].write_bytes(b"existing\n")
    with pytest.raises(ValueError, match="untouched"):
        normalize(files)
    assert files[2].read_bytes() == b"existing\n"


def test_original_cannot_be_used_as_output(files):
    original = files[0].read_bytes()
    with pytest.raises(ValueError, match="separate"):
        normalize((files[0], files[1], files[0], files[3]))
    assert files[0].read_bytes() == original


def test_native_ownership_and_candidate_gates_still_apply(files):
    normalize(files)
    with pytest.raises(ValueError):
        check_pairs(files[2], {"a": "A", "b": "B"}, {"a": 0, "b": 0}, 2)
    with pytest.raises(ValueError):
        check_pairs(files[2], {"a": "A", "b": "B", "c": "C"}, {"a": 0, "b": 1, "c": 0}, 2)


@pytest.mark.parametrize("field,value", [("status", "unverified"), ("removed_rows", 2),
    ("removed_rows", True), ("original_rows", 4), ("native_inference_rerun", True)])
def test_invalid_projection_proof_refused(files, field, value):
    proof = normalize(files)
    proof[field] = value
    with pytest.raises(ValueError):
        export.verify_projection(proof)


def test_rehashed_but_changed_derived_content_fails_byte_comparison(files):
    proof = normalize(files)
    files[2].write_bytes(files[2].read_bytes().replace(b"\tc\t", b"\td\t"))
    proof["derived"]["pairs"] = record(files[2])
    with pytest.raises(ValueError, match="differs beyond"):
        export.verify_projection(proof)


def test_wrong_skip_offset_cannot_remove_a_distinct_prediction(files):
    proof = normalize(files)
    proof["skips"]["pairs"]["offset"] += len(PAIR)
    with pytest.raises(ValueError, match="identical adjacent"):
        export.verify_projection(proof)


def test_omitting_first_twin_is_rejected_without_a_matching_predecessor(files):
    proof = normalize(files)
    proof["skips"]["pairs"]["offset"] -= len(PAIR)
    with pytest.raises(ValueError, match="identical adjacent"):
        export.verify_projection(proof)


def test_original_change_after_normalization_invalidates_proof(files):
    proof = normalize(files)
    files[0].write_bytes(files[0].read_bytes() + OTHER)
    with pytest.raises(ValueError):
        export.verify_projection(proof)
