import pytest

from benchmark_tools import prepare_ob_sequence_features as features


def test_balanced_and_concentrated_sequences():
    balanced = features.sequence_features(features.CANONICAL * 5)
    assert balanced["normalized_entropy20"] == pytest.approx(1.)
    assert balanced["largest_canonical_fraction"] == .05
    assert balanced["short_lt100"] is False
    assert balanced["composition_concentrated"] is False
    concentrated = features.sequence_features("a" * 99)
    assert concentrated["normalized_entropy20"] == 0
    assert concentrated["largest_canonical_fraction"] == 1
    assert concentrated["short_lt100"] is True
    assert concentrated["composition_concentrated"] is True


@pytest.mark.parametrize("sequence", ["", "---...**", "X" * 100])
def test_no_canonical_sequence_has_missing_entropy(sequence):
    row = features.sequence_features(sequence)
    assert row["normalized_entropy20"] is None
    assert row["composition_concentrated"] is None
    assert row["largest_canonical_fraction"] is None


def test_ambiguous_letters_and_symbols_remain_explicit():
    row = features.sequence_features("aCXXUO*-?.")
    assert row["raw_length"] == 10
    assert row["residue_length"] == 6
    assert row["canonical_length"] == 2
    assert row["noncanonical_letter_count"] == 4
    assert row["canonical_fraction"] == pytest.approx(1/3)
    assert row["gap_symbols"] == 2
    assert row["stop_symbols"] == 1
    assert row["other_symbols"] == 1
    assert row["composition_concentrated"] is None


@pytest.mark.parametrize("sequence", ["A" * 19, "A" * 20 + "X" * 3])
def test_insufficient_composition_evidence_is_not_negative(sequence):
    assert features.sequence_features(sequence)["composition_concentrated"] is None


def test_non_ascii_rejected():
    with pytest.raises(ValueError, match="Non-ASCII"):
        features.sequence_features("A\u00e9")


def test_no_overwrite_before_manifest_access(tmp_path, monkeypatch):
    monkeypatch.setattr(features, "read_frozen", lambda *a: pytest.fail("Read manifest before overwrite guard"))
    with pytest.raises(FileExistsError):
        features.prepare(tmp_path, tmp_path)


@pytest.mark.parametrize("problem", ["duplicate_gene", "changed_input", "incomplete"])
def test_bad_input_cannot_produce_completed_manifest(tmp_path, monkeypatch, problem):
    inputs = []
    for i in range(12):
        path = tmp_path / f"species{i:02d}.fa"
        gene = "g0" if problem == "duplicate_gene" and i == 1 else f"g{i}"
        path.write_text(f">{gene}\nACDE\n")
        inputs.append(features.file_provenance(path))
    monkeypatch.setattr(features, "read_frozen", lambda *a: {"fasta_inputs": inputs})
    if problem == "changed_input":
        (tmp_path / "species00.fa").write_text(">g0\nAAAA\n")
    output = tmp_path / "out"
    with pytest.raises(ValueError):
        features.prepare(tmp_path, output)
    assert not (output / "manifest.json").exists()
