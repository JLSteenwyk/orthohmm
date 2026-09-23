import hashlib

import pytest

from benchmark_tools.audit_unisave_fragment_source import choose_version, inspect_entry


def row(version, sequence=1, first="12-Aug-2020", last="12-Aug-2020"):
    return dict(accession="P12345", entryVersion=version, sequenceVersion=sequence,
                firstReleaseDate=first, lastReleaseDate=last)


def test_baseline_annotation_not_current():
    rows = [row(3, first="02-Dec-2020", last="02-Dec-2020"), row(2),
            row(1, first="17-Jun-2020", last="17-Jun-2020")]
    selected, kind = choose_version(dict(results=rows), "P12345", 1)
    assert selected["entryVersion"] == 2
    assert kind == "baseline_release"


def test_corrected_sequence_earliest_later_version():
    rows = [row(4, 2, "10-Feb-2021", "10-Feb-2021"), row(3, 2, "02-Dec-2020", "02-Dec-2020"), row(2)]
    selected, kind = choose_version(dict(results=rows), "P12345", 2)
    assert selected["entryVersion"] == 3
    assert kind == "later_sequence_version"


@pytest.mark.parametrize("rows,version", [([], 1), ([row(1), row(1)], 1),
    ([row(1)], 2), ([row(1, first="17-Jun-2020", last="17-Jun-2020")], 1),
    ([dict(row(1), accession="OTHER")], 1)])
def test_ambiguous_or_missing_history(rows, version):
    with pytest.raises(ValueError):
        choose_version(dict(results=rows), "P12345", version)


def entry(flag):
    return ("ID   TEST_HUMAN              Reviewed;         3 AA.\n"
        "AC   P12345;\nDT   01-JAN-2000, integrated into UniProtKB/Swiss-Prot.\n"
        "DT   01-JAN-2000, sequence version 1.\nDT   12-AUG-2020, entry version 2.\n"
        "DE   RecName: Full=Test;\n" + ("DE   Flags: Fragment;\n" if flag else "") +
        "OS   Homo sapiens (Human).\nOC   Eukaryota.\nOX   NCBI_TaxID=9606;\n"
        "SQ   SEQUENCE   3 AA;  100 MW;  0000000000000000 CRC64;\n     AAA\n//\n")


@pytest.mark.parametrize("flag", [False, True])
def test_exact_identity_and_fragment_flag(flag):
    args = ["P12345", 2, 1, hashlib.sha256(b"AAA").hexdigest(), "9606"]
    result = inspect_entry(entry(flag), *args)
    assert result["fragment_flag"] is flag
    assert "does not prove completeness" in result["limitation"]
    for i, value in enumerate(["OTHER", 3, 2, hashlib.sha256(b"BBB").hexdigest(), "8364"]):
        wrong = args.copy()
        wrong[i] = value
        with pytest.raises(ValueError, match="exact input identity"):
            inspect_entry(entry(flag), *wrong)


@pytest.mark.parametrize("feature", ["NON_TER", "NON_CONS"])
def test_incomplete_sequence_feature_retained(feature):
    text = entry(False).replace("SQ   SEQUENCE", f"FT   {feature:<15} 1\nSQ   SEQUENCE")
    result = inspect_entry(text, "P12345", 2, 1, hashlib.sha256(b"AAA").hexdigest(), "9606")
    assert result["incomplete_sequence_features"][0]["type"] == feature
    assert result["fragment_flag"] is False
