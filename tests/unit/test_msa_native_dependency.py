import pytest

from orthohmm.search import msa_center_star


def test_missing_native_alignment_has_actionable_error(monkeypatch):
    monkeypatch.setattr(msa_center_star, "_pa_lib", None)
    error = OSError("missing library or shared dependency")

    def fail_load(_path):
        raise error

    monkeypatch.setattr(msa_center_star.ctypes.cdll, "LoadLibrary", fail_load)
    with pytest.raises(RuntimeError, match="MSA-profile expansion requires") as caught:
        msa_center_star._load_pair_align()
    assert "OpenMP" in str(caught.value)
    assert "Numba search fallback does not" in str(caught.value)
    assert caught.value.__cause__ is error
    assert msa_center_star._pa_lib is None


def test_cached_alignment_library_is_unchanged(monkeypatch):
    library = object()
    monkeypatch.setattr(msa_center_star, "_pa_lib", library)
    assert msa_center_star._load_pair_align() is library
