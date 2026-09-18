import gzip
import json

import pytest

from benchmark_tools.prepare_qfo_corrected_fastoma_pairs import validate_admission, convert


def admission(root):
    pairs = {"path": str(root / "benchmarks/results/qfo_corrected_fastoma_v1/output/orthologs.tsv.gz"),
             "bytes": 10, "sha256": "native"}
    inputs = [{"path": str(root / f"benchmarks/work/qfo_corrected_fastoma_inputs_20260918/proteome/sp{i}.fa"),
               "bytes": 1, "sha256": "input"} for i in range(78)]
    return {"status": "corrected_fastoma_native_evidence_admitted", "accuracy_evaluated": False,
            "publication_ready": False, "content": {"input_proteins": 984137, "species": 78, "native_pair_rows": 3},
            "scheduler": {"State": "COMPLETED", "ExitCode": "0:0", "AllocCPUS": "180", "ReqMem": "720G", "NodeList": "bizon"},
            "native_pairs": pairs, "input_fastas": inputs, "checked_records": [pairs, *inputs]}


def test_valid_admission(tmp_path):
    assert validate_admission(admission(tmp_path), tmp_path) == 3


@pytest.mark.parametrize("change", ["status", "accuracy", "publication", "genes", "species", "count",
                                   "bool_count", "cpu", "memory", "failed", "pair_path", "input_path",
                                   "missing_input", "duplicate_input", "unbound_input", "unbound_pairs"])
def test_invalid_admission_rejected(tmp_path, change):
    data = admission(tmp_path)
    if change == "status":
        data["status"] = "pending"
    elif change == "accuracy":
        data["accuracy_evaluated"] = True
    elif change == "publication":
        data["publication_ready"] = True
    elif change in ("genes", "species"):
        data["content"]["input_proteins" if change == "genes" else "species"] -= 1
    elif change in ("count", "bool_count"):
        data["content"]["native_pair_rows"] = 0 if change == "count" else True
    elif change in ("cpu", "memory", "failed"):
        key, value = {"cpu": ("AllocCPUS", "32"), "memory": ("ReqMem", "700G"), "failed": ("State", "FAILED")}[change]
        data["scheduler"][key] = value
    elif change == "pair_path":
        data["native_pairs"]["path"] = "/old/pairs.gz"
    elif change == "input_path":
        data["input_fastas"][0]["path"] = "/old/sp.fa"
    elif change == "missing_input":
        data["input_fastas"].pop()
    elif change == "duplicate_input":
        data["input_fastas"][-1] = data["input_fastas"][0]
    elif change == "unbound_input":
        data["checked_records"].pop()
    else:
        data["checked_records"].pop(0)
    with pytest.raises(ValueError):
        validate_admission(data, tmp_path)


def files(tmp_path, rows="b\ta\na\tb\nc\tb\n", mapping=("a", "b", "c"), count=3):
    native = tmp_path / "native.gz"
    native.write_bytes(gzip.compress(rows.encode()))
    reference = tmp_path / "mapping.json.gz"
    reference.write_bytes(gzip.compress(json.dumps({"mapping": dict.fromkeys(mapping, 1)}).encode()))
    data = {"native_pairs": {"path": str(native)}, "content": {"native_pair_rows": count}}
    return data, {"path": str(reference)}


def test_distinct_conversion_retains_native_semantics(tmp_path):
    data, mapping = files(tmp_path)
    counts, raw, filtered = convert(data, {"a": "s1", "b": "s2", "c": "s1"}, tmp_path, mapping)
    assert counts == {"native_rows": 3, "distinct_pairs": 2, "duplicate_relations": 1}
    assert raw.read_text() == "a\tb\nb\tc\n" == filtered.read_text()
    assert not (tmp_path / "pairs.tsv").exists()


@pytest.mark.parametrize("failure", ["mapping_loss", "count", "unknown", "malformed"])
def test_invalid_conversion_keeps_only_partial_files(tmp_path, failure):
    args = {"mapping_loss": {"mapping": ("a", "b")}, "count": {"count": 4},
            "unknown": {"rows": "a\tforeign\n"}, "malformed": {"rows": "a\tb\textra\n"}}[failure]
    data, mapping = files(tmp_path, **args)
    with pytest.raises(ValueError):
        convert(data, {"a": "s1", "b": "s2", "c": "s1"}, tmp_path, mapping)
    assert not (tmp_path / "pairs.tsv").exists()


def test_existing_partial_never_overwritten(tmp_path):
    data, mapping = files(tmp_path)
    partial = tmp_path / "pairs.partial.tsv"
    partial.write_text("preserve")
    with pytest.raises(FileExistsError):
        convert(data, {"a": "s1", "b": "s2", "c": "s1"}, tmp_path, mapping)
    assert partial.read_text() == "preserve"
