from pathlib import Path
from types import SimpleNamespace

import pytest

from tests.conftest import _swiss_source_options, pytest_configure


def config(options):
    markers = []
    return SimpleNamespace(getoption=options.get,
                           addinivalue_line=lambda key, value: markers.append((key, value)),
                           markers=markers)


@pytest.mark.parametrize("kind", ["duplication", "fragment", "descriptive", "identity"])
def test_default_retains_original_exporter_arguments(kind):
    assert _swiss_source_options(config({}), kind) == {}


@pytest.mark.parametrize("kind", ["duplication", "fragment"])
def test_explicit_pair_is_forwarded_without_repinning_or_reading(kind):
    path, digest = Path("not-yet-transferred/bindings.json"), "1" * 64
    options = {f"--swiss-{kind}-bindings": path, f"--swiss-{kind}-bindings-sha256": digest}
    assert _swiss_source_options(config(options), kind) == dict(source_bindings=path, source_bindings_sha=digest)
    other = "fragment" if kind == "duplication" else "duplication"
    assert _swiss_source_options(config(options), other) == {}
    assert _swiss_source_options(config(options), "descriptive") == {}
    assert _swiss_source_options(config(options), "identity") == {}


@pytest.mark.parametrize("kind", ["duplication", "fragment"])
@pytest.mark.parametrize("partial", ["path", "digest"])
def test_partial_pair_is_usage_error_before_collection(kind, partial):
    name = f"--swiss-{kind}-bindings" + ("-sha256" if partial == "digest" else "")
    options = {name: "value"}
    with pytest.raises(pytest.UsageError, match="required together"):
        _swiss_source_options(config(options), kind)
    with pytest.raises(pytest.UsageError, match="required together"):
        pytest_configure(config(options))


def test_both_panel_pairs_are_independent_and_markers_preserved():
    options = {"--swiss-duplication-bindings": Path("dup.json"),
               "--swiss-duplication-bindings-sha256": "2" * 64,
               "--swiss-fragment-bindings": Path("frag.json"),
               "--swiss-fragment-bindings-sha256": "3" * 64}
    selected = config(options)
    pytest_configure(selected)
    assert _swiss_source_options(selected, "duplication")["source_bindings_sha"] == "2" * 64
    assert _swiss_source_options(selected, "fragment")["source_bindings_sha"] == "3" * 64
    assert selected.markers == [("markers", "integration: mark as integration test"),
                                ("markers", "slow: mark as slow test")]


def test_wrong_binding_digest_is_not_replaced_by_a_computed_one():
    selected = config({"--swiss-fragment-bindings": Path("missing.json"),
                       "--swiss-fragment-bindings-sha256": "0" * 64})
    assert _swiss_source_options(selected, "fragment")["source_bindings_sha"] == "0" * 64


def test_unknown_exporter_kind_rejected():
    with pytest.raises(ValueError, match="Unknown"):
        _swiss_source_options(config({}), "typo")
