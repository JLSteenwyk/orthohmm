import pytest
import xml.etree.ElementTree as ET

from benchmark_tools.summarize_cpm_memcheck import digest, summarize


XML = """<valgrindoutput><protocoltool>memcheck</protocoltool>
<status><state>RUNNING</state></status>
<error><unique>0x0</unique><kind>InvalidRead</kind><what>read</what>
<stack><frame><fn>startup</fn><obj>python</obj></frame></stack>
<auxwhat>allocation</auxwhat><stack><frame><fn>malloc</fn></frame></stack></error>
<error><unique>0x1</unique><kind>Leak_DefinitelyLost</kind>
<stack><frame><fn>malloc</fn><obj>libc</obj></frame></stack></error>
<status><state>FINISHED</state></status>
<errorcounts><pair><unique>0x0</unique><count>3</count></pair></errorcounts>
</valgrindoutput>"""


def test_retains_stacks_and_separates_records_from_occurrences(tmp_path):
    path = tmp_path / "input.xml"
    path.write_text(XML)
    result = summarize(path, digest(path))
    assert result["error_records_by_kind"] == {"InvalidRead": 1, "Leak_DefinitelyLost": 1}
    assert result["reported_occurrences"] == {"0x0": 3}
    assert len(result["nonleak_errors"][0]["stacks"]) == 2
    assert result["nonleak_errors"][0]["auxiliary"] == ["allocation"]
    assert result["leak_top_frames"][0]["records"] == 1
    assert not result["accuracy_admitted"] and not result["publication_ready"]


@pytest.mark.parametrize("old,new", [
    ("memcheck", "other"), ("FINISHED", "RUNNING"),
    ("<unique>0x1", "<unique>0x0"), ("<count>3", "<count>0"),
    ("<pair><unique>0x0", "<pair><unique>0x9"),
])
def test_rejects_invalid_reports(tmp_path, old, new):
    path = tmp_path / "input.xml"
    path.write_text(XML.replace(old, new))
    with pytest.raises(ValueError):
        summarize(path, digest(path))


def test_rejects_wrong_digest(tmp_path):
    path = tmp_path / "input.xml"
    path.write_text(XML)
    with pytest.raises(ValueError, match="Changed"):
        summarize(path, "wrong")


def test_rejects_truncated_xml(tmp_path):
    path = tmp_path / "input.xml"
    path.write_text(XML[:-20])
    with pytest.raises(ET.ParseError):
        summarize(path, digest(path))


def test_rechecks_identity_after_extraction(tmp_path, monkeypatch):
    from benchmark_tools import summarize_cpm_memcheck as module
    path = tmp_path / "input.xml"
    path.write_text(XML)
    checks = iter(["original", "changed"])
    monkeypatch.setattr(module, "digest", lambda p: next(checks))
    with pytest.raises(ValueError, match="during extraction"):
        module.summarize(path, "original")


def test_fatal_signal_is_retained_not_admitted(tmp_path):
    path = tmp_path / "input.xml"
    path.write_text(XML.replace("</valgrindoutput>", "<fatal_signal/></valgrindoutput>"))
    result = summarize(path, digest(path))
    assert result["fatal_signal_records"] == 1
    assert not result["accuracy_admitted"]
