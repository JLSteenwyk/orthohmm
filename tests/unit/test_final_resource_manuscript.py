import json
from pathlib import Path
import re


RESULTS = Path(__file__).resolve().parents[2] / 'benchmark_tools/results'
OLD = RESULTS / 'PUBLICATION_MAIN_TEXT_20261003.md'
MAIN = RESULTS / 'PUBLICATION_MAIN_TEXT_20261004.md'


def between(text, start, end):
    return text.split(start, 1)[1].split(end, 1)[0]


def test_final_resource_section_is_exact_generated_text():
    text = MAIN.read_text()
    section = (RESULTS / 'threadripper_shared_resource_section_20261004_v27/resource_section.md').read_text()
    report = json.loads((RESULTS / 'threadripper_shared_panel_snapshot_20261004_v27/panel.json').read_text())
    assert section.strip() in text
    assert report['all_planned_attempts_reviewed'] is True
    assert 'all 27 attempts, with 24 eligible measurements' in text
    assert 'this is not 27 eligible measurements' in text
    assert 'Both retain null native endpoints' in text
    assert 'None of these failures was retried' in text
    assert 'unknown and potentially method dependent' in text
    assert '### Shared-Host Resource Panel (Interim)' not in text
    assert 'panel is in progress' not in text
    assert 'No submission-ready release or archival DOI is claimed' in text


def test_resource_update_does_not_change_scientific_sections_or_citations():
    old, new = OLD.read_text(), MAIN.read_text()
    for start, end in (
        ('## Methods\n', '### Shared-Host Resource Measurement\n'),
        ('### Shared-Host Resource Measurement\n', '## Results\n'),
        ('## Results\n', '### Shared-Host Resource Panel'),
        ('## Discussion And Limitations\n', 'Recovered OrthoMCL results retain'),
    ):
        assert between(old, start, end) == between(new, start, end)
    identifiers = lambda text: set(re.findall(r'@([A-Za-z0-9_]+)', text))
    assert identifiers(new) == identifiers(old)
    assert len(identifiers(new)) == 18
    assert '7f3a9e4' in new
    assert 'VGNC dependence' in new and 'Valid paired uncertainty for GO/EC, FAS' in new
    assert 'rather than family-disjoint validation' in new


def test_final_resource_links_and_scope_are_present():
    text = MAIN.read_text()
    for name in (
        'threadripper_shared_panel_snapshot_20261004_v27/panel.json',
        'threadripper_shared_panel_snapshot_20261004_v27/attempts.tsv',
        'threadripper_shared_panel_snapshot_20261004_v27/cells.tsv',
        'threadripper_shared_resource_figure_20261004_v26/shared_threadripper_resources.pdf',
        'threadripper_shared_resource_section_20261004_v27/manifest.json',
        'FINAL_RESOURCE_REPORTING_COMPONENT_20261004.md',
        'PUBLICATION_EXPORT_INTEGRATION_20261004.md',
    ):
        assert name in text
        assert (RESULTS / name).is_file()
    assert 'No background overhead is subtracted' in text
    assert 'Preparation, conversion and' in text
    assert 'not pure algorithm RSS' in text
    assert 'taxon-count scaling from proteome composition' in text
    assert 'Whole-study executable/versioned' in text
    assert 'archival packaging remains open' in text
