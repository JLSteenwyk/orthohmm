from pathlib import Path
import re
import shutil
import subprocess
import json
import zipfile

import pytest


ROOT = Path(__file__).resolve().parents[2] / "manuscript/hmm_frontier_20261020"
FIRST_PERSON = re.compile(r"\b(?:I|me|my|mine|myself|we|us|our|ours|ourselves)\b", re.I)


def manuscript():
    return (ROOT / "manuscript.md").read_text()


def test_manuscript_and_captions_have_no_first_person():
    assert FIRST_PERSON.findall(manuscript()) == []


def test_no_parenthetical_em_or_en_dashes():
    assert not any(character in manuscript() for character in ("\u2013", "\u2014"))


def test_running_prose_has_no_colons():
    body = manuscript().split("---", 2)[2].split("## References")[0]
    lines = [line for line in body.splitlines() if not line.startswith("#")]
    prose = re.sub(r"https?://[^>\s]+", "", "\n".join(lines))
    assert ":" not in prose


def test_byline_has_only_confirmed_author_and_no_inferred_affiliation():
    metadata = manuscript().split("---", 2)[1]
    assert 'author: "Jacob L. Steenwyk"' in metadata
    assert len(re.findall(r"^author:", metadata, re.M)) == 1
    assert not re.search(r"\b(?:affiliation|email|institute|corresponding)\b", metadata, re.I)
    assert "Jacob L. Steenwyk is the sole author." in manuscript()
    assert "Author order" not in manuscript()


def test_every_numerical_token_and_figure_link_matches_original_draft():
    with zipfile.ZipFile(ROOT / "OrthoHMM_manuscript_20261020.zip") as archive:
        original = archive.read("OrthoHMM_manuscript_20261020/manuscript.md").decode()
    assert re.findall(r"[-+]?\d+(?:\.\d+)?", manuscript()) == re.findall(r"[-+]?\d+(?:\.\d+)?", original)
    assert re.findall(r"!\[[^\]]*\]\(([^)]+)\)", manuscript()) == re.findall(r"!\[[^\]]*\]\(([^)]+)\)", original)
    assert set(re.findall(r"@([a-zA-Z0-9_-]+)", manuscript())) == set(re.findall(r"@([a-zA-Z0-9_-]+)", original))


@pytest.mark.skipif(shutil.which("pandoc") is None, reason="pandoc required")
def test_pandoc_interprets_exactly_one_author():
    document = json.loads(subprocess.check_output(["pandoc", "--from=markdown", "--to=json", str(ROOT / "manuscript.md")], text=True))
    author = document["meta"]["author"]
    assert author["t"] == "MetaInlines"
    assert " ".join(n["c"] for n in author["c"] if n["t"] == "Str") == "Jacob L. Steenwyk"
