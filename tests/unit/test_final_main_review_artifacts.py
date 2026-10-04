import hashlib
import json
from pathlib import Path

import fitz


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / 'benchmark_tools/results'


def test_final_main_review_is_closed_and_has_bounded_scope():
    receipt = json.loads((RESULTS / 'publication_main_final_visual_review_20261004.json').read_text())
    assert receipt['page_count'] == 12
    assert receipt['all_twelve_pages_inspected'] is True
    assert receipt['bounds_violations'] == []
    assert receipt['observed_clipping_or_overlap'] is False
    assert receipt['resource_tables_readable'] is True
    assert len(receipt['citation_ids']) == 18
    assert receipt['untracked_targets'] == []
    assert receipt['figures_linked_not_embedded'] is True
    assert receipt['other_linked_figures_newly_revalidated'] is False
    for key in ('native_inference_or_scoring_rerun', 'benchmark_scores_or_defaults_changed',
                'new_study_archive_built', 'public_release_or_deposition_executed', 'publication_ready'):
        assert receipt[key] is False
    assert receipt['direct_input_hashes_and_git_blobs_checked'] is True
    assert all(row['revision'] == receipt['render_time_source_commit']
               for row in receipt['checked_committed_direct_inputs'])
    artifacts = receipt['closed_review_artifacts']
    assert len(artifacts) == len({row['path'] for row in artifacts}) == 17
    for row in artifacts:
        path = ROOT / row['path']
        assert path.resolve().is_relative_to(ROOT)
        raw = path.read_bytes()
        assert len(raw) == row['bytes']
        assert hashlib.sha256(raw).hexdigest() == row['sha256']


def test_final_pdf_keeps_all_displayed_resource_values_and_caveats():
    with fitz.open(RESULTS / 'publication_main_print_20261004/document.pdf') as document:
        assert len(document) == 12
        text = ' '.join(' '.join(page.get_text().split()) for page in document)
    section = (RESULTS / 'threadripper_shared_resource_section_20261004_v27/resource_section.md').read_text()
    rows = [line.split('|')[1:-1] for line in section.splitlines()
            if line.startswith('| Ortho')]
    assert len(rows) == 27
    values = [row[-1].strip() for row in rows]
    assert text.count('Unavailable') == values.count('Unavailable') == 9
    for value in values:
        if value != 'Unavailable':
            assert value in text
    for phrase in (
        '27 of 27 reviewed attempts', '25 with measured native resources',
        '24 eligible observations', 'unknown and potentially method dependent',
        'not isolated comparative timing', 'No background overhead is subtracted',
        'this is not 27 eligible measurements', 'Unreviewed indices: none',
        'Original TreeFam-A family mappings', 'not submission-ready',
    ):
        assert phrase.lower() in text.lower()
