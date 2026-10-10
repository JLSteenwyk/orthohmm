from copy import deepcopy
import importlib.util
from pathlib import Path
import sys

import pytest

from benchmark_tools import probe_native_fas_context as base


@pytest.fixture
def corrected(monkeypatch):
    monkeypatch.setitem(sys.modules, "probe_native_fas_context", base)
    path = Path(base.__file__).with_name("probe_native_fas_context_lowercase.py")
    spec = importlib.util.spec_from_file_location("context_lowercase_test", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_only_annotation_tool_keys_are_normalized(corrected):
    original, owners = base.fixtures()
    expected = deepcopy(original)
    for doc in expected.values():
        for gene, features in doc["feature"].items():
            doc["feature"][gene] = {key.lower(): value for key, value in features.items()}
    actual, actual_owners = corrected.fixtures()
    assert actual == expected and actual_owners == owners
    assert original == base.fixtures()[0]
    assert actual["T1"]["feature"]["C"]["pfam"] == original["T1"]["feature"]["C"]["Pfam"]
    assert all(key == key.lower() for doc in actual.values() for gene in doc["feature"].values() for key in gene)


def test_original_kernels_are_bound_and_reused(corrected):
    assert base.fingerprint(base.__file__)["sha256"] == corrected.BASE_SHA
    assert corrected.base is base
    assert corrected.base.contexts is base.contexts
    assert corrected.base.execute_context is base.execute_context
    assert corrected.base.compare is base.compare


def test_main_reuses_limits_and_restores_original_probe(corrected, monkeypatch):
    original = base.probe
    def main():
        assert base.probe is corrected.probe
        raise FileExistsError("occupied output")
    monkeypatch.setattr(base, "main", main)
    with pytest.raises(FileExistsError):
        corrected.main()
    assert base.probe is original
