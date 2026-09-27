import pytest

from benchmark_tools.audit_ob_orthomcl_provenance import duration, parse_timestamp, mapping, read_index, flag


def test_distinct_years_and_intervals():
    april = duration("Mon Apr 20 04:03:24 PM EDT 2026","Thu Apr 23 07:12:39 PM EDT 2026")
    july = duration("Fri Jul 25 01:27:53 PM EDT 2025","Tue Jul 29 08:42:55 PM EDT 2025")
    assert april["seconds"] == 270555
    assert july["seconds"] == 371702
    assert july["start"].startswith("2025-")


@pytest.mark.parametrize("date", ["Mon Jul 25 01:27:53 PM EDT 2025", "Fri Jul 25 01:27:53 PM EST 2025"])
def test_reject_inconsistent_or_implicit_timezone(date):
    with pytest.raises(ValueError):
        parse_timestamp(date)


def test_species_alias_matched_by_complete_gene_set(tmp_path):
    p = tmp_path / "all.gg"
    p.write_text("Canis: a b\nother: c\n")
    result = mapping(p,{"dog.fa":{"a","b"},"other.fa":{"c"}})
    assert result["Canis"]["input_species"] == "dog.fa"


@pytest.mark.parametrize("text", ["s: a a\n", "s: a\ns: b\n", "s: a\n", "s: z\n"])
def test_bad_species_mapping(tmp_path,text):
    p = tmp_path / "all.gg"
    p.write_text(text)
    with pytest.raises(ValueError):
        mapping(p,{"one":{"a"},"two":{"b"}})


@pytest.mark.parametrize("text", ["0 a\n0 b\n", "0 a\n1 a\n", "1 a\n", "0 z\n"])
def test_bad_index(tmp_path,text):
    p = tmp_path / "idx"
    p.write_text(text)
    with pytest.raises(ValueError):
        read_index(p,{"a","b"})


def test_duplicate_command_flag_rejected():
    with pytest.raises(ValueError):
        flag(["blastall","-a","8","-a","32"],"-a")
