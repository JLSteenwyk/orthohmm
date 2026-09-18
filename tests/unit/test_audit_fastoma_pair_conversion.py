import gzip

import pytest

from benchmark_tools.audit_fastoma_pair_conversion import compare


@pytest.mark.parametrize("retained", ["a\tb\n", "", "a\tb\na\tb\n", "b\ta\n", "a\tc\n"])
def test_complete_exact_comparison(tmp_path, retained):
    native, old = tmp_path / "native.gz", tmp_path / "retained.tsv"
    native.write_bytes(gzip.compress(b"b\ta\n"))
    old.write_text(retained)
    if retained == "a\tb\n":
        assert compare(native, old, {"a": "s1", "b": "s2"}) == 1
    else:
        with pytest.raises(ValueError):
            compare(native, old, {"a": "s1", "b": "s2"})
