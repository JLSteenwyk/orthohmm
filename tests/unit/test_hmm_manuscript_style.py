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


def without_generation_timing(text):
    for heading in ("Small-dataset generation timing", "Generation timing protocol"):
        text = re.sub(r"### " + heading + r"\n.*?(?=\n### |\n## )", "", text, flags=re.S)
    return re.sub(r"::: \{\.figure-sheet\}\n## Figure 6\..*?\n:::\n", "", text, flags=re.S)


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
    def scientific_numbers(text):
        text = re.sub(r"@([a-zA-Z0-9_-]+)", "", text)
        text = re.sub(r"assets_v\d+/", "assets/", text)
        text = text.replace("OrthoHMM2", "OrthoHMM")
        return re.findall(r"[-+]?\d+(?:\.\d+)?", text)

    retained = without_generation_timing(manuscript())
    assert scientific_numbers(retained) == scientific_numbers(original)
    links = lambda text: [Path(p).name for p in re.findall(r"!\[[^\]]*\]\(([^)]+)\)", text)]
    assert links(retained) == links(original)
    original_citations = set(re.findall(r"@([a-zA-Z0-9_-]+)", original))
    assert set(re.findall(r"@([a-zA-Z0-9_-]+)", manuscript())) == original_citations | {
        "qfo2024", "schlicker2006", "Dosch2023FAS"}


def test_original_orthohmm_preprint_is_credited_without_renaming_its_authors_or_title():
    references = json.loads((ROOT / "references.csl.json").read_text())
    preprint = next(r for r in references if r["id"] == "orthohmm2024preprint")
    with zipfile.ZipFile(ROOT / "OrthoHMM_manuscript_20261020.zip") as archive:
        original = json.loads(archive.read("OrthoHMM_manuscript_20261020/references.csl.json"))
    assert preprint == next(r for r in original if r["id"] == preprint["id"])
    introduction = manuscript().split("## Introduction")[1].split("## Results")[0]
    assert "original OrthoHMM" in introduction
    assert "[@orthohmm2024preprint]" in introduction
    assert "OrthoHMM2 extends that framework" in introduction
    assert 'title: "OrthoHMM2: HMM-centered orthology inference"' in manuscript()


def test_all_citations_resolve_to_unique_bibliography_entries():
    references = json.loads((ROOT / "references.csl.json").read_text())
    ids = [r["id"] for r in references]
    dois = [r["DOI"].lower() for r in references]
    assert len(ids) == len(set(ids))
    assert len(dois) == len(set(dois))
    assert set(re.findall(r"@([a-zA-Z0-9_-]+)", manuscript())) == set(ids)


def test_citation_counts_include_repeated_mentions_and_grouped_citations():
    mentions = re.findall(r"@([a-zA-Z0-9_-]+)", manuscript())
    groups = re.findall(r"\[[^\]]*?@[^\]]*\]", manuscript())
    assert len(set(mentions)) == 17
    assert len(mentions) == 20
    assert len(groups) == 14


def test_journal_references_include_renderable_pages_or_article_locators():
    references = json.loads((ROOT / "references.csl.json").read_text())
    for entry in references:
        if entry["type"] == "article-journal":
            assert entry.get("volume") and entry.get("page"), entry["id"]
    by_id = {r["id"]: r for r in references}
    for key, locator in {"orthofinder2019": "238", "orthofinder2026correction": "2110",
                         "proteinortho2023": "1322477", "sonicparanoid2024": "195",
                         "Dosch2023FAS": "btad226"}.items():
        assert by_id[key]["page"] == locator


def test_qfo_consortium_order_matches_publisher_not_crossref_order():
    references = json.loads((ROOT / "references.csl.json").read_text())
    entry = next(r for r in references if r["id"] == "qfo2016")
    assert entry["author"][0]["family"] == "Altenhoff"
    assert entry["author"][21]["family"] == "Muffato"
    assert entry["author"][22] == {"literal": "Quest for Orthologs consortium"}
    assert entry["author"][23]["family"] == "Gabald\u00f3n"


def test_functional_similarity_attribution_identifies_qfo_adaptation():
    assert "QfO adaptation of Schlicker similarity" in manuscript()


def test_citation_audit_covers_every_bibliography_entry():
    report = json.loads((ROOT / "validation/CITATION_AUDIT.json").read_text())
    references = json.loads((ROOT / "references.csl.json").read_text())
    assert {r["id"] for r in report["references"]} == {r["id"] for r in references}
    assert all(r["doi_and_title_verified"] and r["claim_fit_reviewed"] for r in report["references"])
    assert all(r["registry_record_url"] and r["claim_scope"] for r in report["references"])


def test_hmm2_figure_labels_preserve_all_scientific_snapshot_values():
    old = json.loads((ROOT / "assets_v3/analysis.json").read_text())
    new = json.loads((ROOT / "assets_v4/analysis.json").read_text())
    for method in old["methods"]:
        method["label"] = method["label"].replace("OrthoHMM ", "OrthoHMM2 ")
    assert old == new
    for name in ("scores_all_29.tsv", "plot_coordinates.tsv", "pareto_membership.tsv", "timing_attempts.tsv"):
        assert (ROOT / "assets_v3" / name).read_bytes() == (ROOT / "assets_v4" / name).read_bytes()


def test_generation_timing_claims_match_all_nine_fresh_observations():
    from benchmark_tools.plot_orthohmm_generations import summarize
    data = json.loads((ROOT / "assets_speed/timing.json").read_text())
    assert data["summary"] == summarize(data)
    assert len(data["runs"]) == 9
    assert data["cpu_budget"] == 16
    assert sum(r["proteins"] for r in data["input"]) == 6432
    section = manuscript().split("### Small-dataset generation timing\n")[1].split("\n### ")[0]
    for row in data["summary"]:
        assert f'{row["elapsed_s"]["median"]:.1f}' in section
    for row in data["summary"][1:]:
        assert f'{row["original_to_modern_elapsed_ratio"]:.2f}' in section
    assert "assets_speed/fig6_generations.png" in manuscript()


def test_generation_timing_outputs_cover_the_same_unchanged_proteins():
    import hashlib
    from benchmark_tools.compare_orthohmm_generations import partition, sequences
    assets = ROOT / "assets_speed"
    data = json.loads((assets / "timing.json").read_text())
    expected = [name for _, ids in sequences(assets / "input") for name in ids]
    for ref in data["input"]:
        path = assets / "input" / Path(ref["path"]).name
        assert hashlib.sha256(path.read_bytes()).hexdigest() == ref["sha256"]
    for row in data["runs"]:
        path = assets / "records" / f'{row["index"]:02d}_{row["method"]}_groups.txt'
        assert partition(path, expected) == row["partition"]
    assert all(r["all_sequences_retained"] for r in data["input_source_integrity"])


def test_generation_phylogeny_really_infers_and_reconciles_without_checkpoint_reuse():
    assets = ROOT / "assets_speed"
    data = json.loads((assets / "timing.json").read_text())
    for row in data["runs"]:
        if row["method"] != "phylogeny":
            continue
        counts = json.loads((assets / "records" / f'{row["index"]:02d}_phylogeny_metrics.json').read_text())["counts"]
        assert counts["phylogeny_reconciled_families"] > 0
        assert counts["phylogeny_species_tree_families"] > 0
        assert counts["phylogeny_checkpoint_hits"] == 0
        assert counts["phylogeny_species_tree_checkpoint_hit"] is False


@pytest.mark.skipif(shutil.which("pandoc") is None, reason="pandoc required")
def test_pandoc_interprets_exactly_one_author():
    document = json.loads(subprocess.check_output(["pandoc", "--from=markdown", "--to=json", str(ROOT / "manuscript.md")], text=True))
    author = document["meta"]["author"]
    assert author["t"] == "MetaInlines"
    assert " ".join(n["c"] for n in author["c"] if n["t"] == "Str") == "Jacob L. Steenwyk"
