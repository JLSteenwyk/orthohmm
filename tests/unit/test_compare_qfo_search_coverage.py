import copy

import numpy as np
import pytest

from benchmark_tools import compare_qfo_search_coverage as module
from orthohmm.accuracy import write_accuracy_checkpoint


def inputs(tmp_path, extra_top=False, wrong_names=False, wrong_species=False):
    metadata = {"a": {"species": "s1"}, "b": {"species": "s1"}, "c": {"species": "s2"}}
    checkpoints = {}
    for label, codes in (("hmm", [0, 1, 3, 8]), ("all_hits", [0, 1, 3, 5, 8]),
                         ("top100", [0, 1, 7] if extra_top else [0, 1, 8])):
        codes = np.asarray(codes)
        names = ["a", "b", "z"] if wrong_names and label == "top100" else sorted(metadata)
        species = [0, 1, 1] if wrong_species and label == "top100" else [3, 3, 0]
        path = write_accuracy_checkpoint(str(tmp_path / label), names, species,
            codes // 3, codes % 3, np.ones(len(codes)))
        checkpoints[label] = module.record(path / "manifest.json")
    return checkpoints, metadata


def test_real_three_checkpoint_comparison(tmp_path):
    checkpoints, metadata = inputs(tmp_path)
    result = module.compare(checkpoints, metadata)
    assert result["overlaps"]["hmm_vs_all_hits"]["all"]["intersection"] == 4
    assert result["overlaps"]["hmm_vs_all_hits"]["all"]["second_only"] == 1
    assert result["overlaps"]["all_hits_vs_top100"]["all"]["second_only"] == 0
    assert result["searches"]["all_hits"]["cross_species_hits"] == 1
    assert result["species_labels"] == ["s1", "s2"]


@pytest.mark.parametrize("problem", ["extra_top", "wrong_names", "wrong_species"])
def test_incompatible_checkpoints_fail(tmp_path, problem):
    checkpoints, metadata = inputs(tmp_path, **{problem: True})
    with pytest.raises(ValueError):
        module.compare(checkpoints, metadata)


def test_pending_admissions_do_not_read_files(tmp_path, monkeypatch):
    def pending(*args):
        raise ValueError("Pending")
    monkeypatch.setattr(module, "completed", pending)
    monkeypatch.setattr(module, "record", lambda *a: pytest.fail("Read unadmitted input"))
    with pytest.raises(ValueError, match="Pending"):
        module.run(tmp_path, "1", "2", tmp_path / "output")
    assert not (tmp_path / "output").exists()


@pytest.mark.parametrize("problem", [None, "source", "status", "conversion", "hits", "genes", "cap", "variant_hits", "accuracy"])
def test_numeric_admission_identity(problem):
    source, conversion_record = {"path": "/source"}, {"path": "/conversion"}
    variants = {label: {"status": "checkpoint_matches_reconstructed_source_hits", "cap": cap,
        "genes": 984137, "hits": count, "accuracy_evaluated": False}
        for label, cap, count in (("all_hits", None, 200), ("top100", 100, 100))}
    conversion = {"hits": 200, "variants": {label: {"audit": {"summary": {"hits": row["hits"]}}}
                                            for label, row in variants.items()}}
    admission = {"status": "corrected_qfo_numeric_source_equivalence_admitted", "numeric_equivalence": True,
        "accuracy_evaluated": False, "publication_ready": False, "source": copy.deepcopy(source),
        "conversion": copy.deepcopy(conversion_record), "genes": 984137, "proteomes": 78,
        "hits": 200, "variants": variants}
    if problem in ("source", "conversion"):
        admission[problem] = {}
    elif problem == "status":
        admission["status"] = "failed"
    elif problem in ("hits", "genes"):
        admission[problem] = 1
    elif problem == "cap":
        admission["variants"]["top100"]["cap"] = 100.0
    elif problem == "variant_hits":
        admission["variants"]["top100"]["hits"] = 99
    elif problem == "accuracy":
        admission["accuracy_evaluated"] = True
    if problem:
        with pytest.raises(ValueError):
            module.validate_numeric(admission, conversion, conversion_record, source)
    else:
        module.validate_numeric(admission, conversion, conversion_record, source)
