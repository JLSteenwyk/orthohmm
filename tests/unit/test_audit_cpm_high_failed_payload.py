import hashlib

import numpy as np
import pytest

from benchmark_tools.audit_cpm_high_failed_payload import arrays


@pytest.mark.parametrize("problem", [None, "bounds", "nan", "dtype", "shape", "names"])
def test_payload_array_audit(tmp_path, problem):
    (tmp_path / "gene_names.txt").write_text("a\na\n" if problem == "names" else "a\nb\n")
    src = np.array([0, 1, 1], dtype=np.int32)
    dst = np.array([1, 0, 1], dtype=np.int32)
    weight = np.array([0, -.1, 1.5], dtype=np.float64)
    if problem == "bounds":
        dst[0] = 2
    elif problem == "nan":
        weight[0] = np.nan
    elif problem == "dtype":
        src = src.astype(np.int64)
    elif problem == "shape":
        weight = weight[:2]
    for name, value in zip(("sources", "targets", "weights"), (src, dst, weight)):
        np.save(tmp_path / (name + ".npy"), value)
    if problem:
        with pytest.raises(ValueError):
            arrays(tmp_path, 2)
    else:
        result = arrays(tmp_path, 2)
        assert result["constructor_int32_bytes_sha256"] == hashlib.sha256(np.column_stack((src, dst)).tobytes()).hexdigest()
        assert result["negative_weights"] == result["zero_weights"] == result["self_edges"] == 1
        assert result["vertices"] == 2 and result["edges"] == 3
