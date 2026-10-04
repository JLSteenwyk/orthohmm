import json
from pathlib import Path
from statistics import median

import matplotlib.pyplot as plt
import pytest

from benchmark_tools.results import plot_shared_prepared_resources_20261004 as plotter

RESULTS = Path(__file__).resolve().parents[2] / 'benchmark_tools/results'


@pytest.fixture
def checkpoint():
    candidates = []
    for path in RESULTS.glob('threadripper_shared_resource_figure_20261004_v*/manifest.json'):
        manifest = json.loads(path.read_text())
        table = Path(manifest['source_results']['path'])
        if table.parent.name.startswith('threadripper_shared_panel_snapshot_20261004_'):
            candidates.append((manifest['reviewed_attempts'], path, table, manifest))
    assert candidates
    _, path, table, manifest = max(candidates, key=lambda item: item[0])
    return json.loads(table.read_text()), path.parent, table, manifest


def test_latest_arithmetic_preserves_failure_endpoints_and_complete_cells(checkpoint):
    report, _, _, _ = checkpoint
    measured = plotter.validate(report)
    reviewed = [r for r in report['runs'] if r['status'] != 'not_yet_reviewed']
    eligible = [r for r in reviewed if r['comparative_timing_eligible'] is True]
    assert report['reviewed_attempts'] == len(reviewed)
    assert [r['index'] for r in reviewed] == list(range(len(reviewed)))
    assert report['resource_reviewed_attempts'] == len(measured)
    assert report['eligible_attempts'] == len(eligible)
    assert report['excluded_attempts'] == [r['index'] for r in reviewed if r['comparative_timing_eligible'] is False]
    assert {0, 17, 20}.issubset(report['excluded_attempts'])
    assert report['pre_native_aborted_indices'] == [17, 20]
    assert report['runs'][0]['resources'] is not None
    assert report['runs'][17]['resources'] is None
    assert report['runs'][20]['resources'] is None
    assert report['all_planned_attempts_reviewed'] is (len(reviewed) == 27)
    assert report['uncontended_timing'] is False
    assert report['publication_ready'] is False
    for cell in report['cells']:
        rows = [r for r in eligible if (r['method'], r['proteomes']) == (cell['method'], cell['proteomes'])]
        assert cell['eligible_repeats'] == len(rows)
        for metric, summary in cell['resources'].items():
            if len(rows) == 3:
                values = [r['resources'][metric] for r in rows]
                assert summary == dict(median=median(values), minimum=min(values), maximum=max(values))
            else:
                assert summary is None


def test_latest_terminal_row_and_handoff_remain_bound(checkpoint):
    report, _, _, _ = checkpoint
    row = report['runs'][report['reviewed_attempts'] - 1]
    summary_path = RESULTS / f"threadripper_shared_attempt_{row['job_id']}.json"
    summary = json.loads(summary_path.read_text())
    assert summary['index'] == row['index']
    assert summary['resources'] == row['resources']
    assert summary['native_outputs']['input_genes'] == {4: 73266, 8: 165168, 12: 251378}[row['proteomes']]
    assert summary['native_outputs']['accuracy_evaluated'] is False
    eligible = row['comparative_timing_eligible']
    if eligible:
        assert summary['original_environment_failures'] == dict(processes={}, pressure={})
    assert len(summary['reviews']) == 4
    for category, pin in summary['reviews'].items():
        review = plotter.executor.read(pin)
        assert review['category'] == category
        assert review['decision'] == ('failed' if category == 'environment' and not eligible else 'passed')
        assert (review['job_id'], review['index']) == (row['job_id'], row['index'])
        plotter.executor.check(review['source'])
    handoff = plotter.executor.read(summary['prepared_handoff'])
    assert handoff['worker_terminal_reviewed'] is True
    assert handoff['prepared_before_release_request'] is True
    assert handoff['native_gate_go'] is True
    assert handoff['request_to_bound_review_seconds'] < 20
    for pin in handoff['evidence']:
        plotter.executor.check(pin)


def test_latest_figure_provenance_bounds_pixels_and_coverage(checkpoint):
    import fitz
    import numpy as np

    report, directory, table, manifest = checkpoint
    assert manifest['source_results'] == plotter.executor.record(table)
    assert manifest['plotter'] == plotter.executor.record(plotter.__file__)
    assert manifest['eligible_attempts'] == report['eligible_attempts']
    for pin in manifest['outputs']:
        plotter.executor.check(pin)
    with fitz.open(directory / 'shared_threadripper_resources.pdf') as document:
        assert len(document) == 1
        page = document[0]
        text = page.get_text()
        assert f"{report['reviewed_attempts']}/27 attempts reviewed" in text
        assert 'unknown and potentially method dependent' in text
        for block in page.get_text('dict')['blocks']:
            for line in block.get('lines', []):
                for span in line['spans']:
                    assert page.rect.contains(fitz.Rect(span['bbox']))
        pixmap = page.get_pixmap(alpha=False)
        pixels = np.frombuffer(pixmap.samples, dtype=np.uint8).reshape(pixmap.height, pixmap.width, 3) / 255.
    figure = plotter.plot(report)
    complete = sum(c['eligible_repeats'] == 3 for c in report['cells'])
    try:
        for axis in figure.axes:
            assert sum(len(c.get_offsets()) for c in axis.collections) == report['resource_reviewed_attempts']
            assert len(axis.lines) == 2 * complete
            left, bottom, width, height = axis.get_position().bounds
            crop = pixels[int((1-bottom-height)*pixels.shape[0]):int((1-bottom)*pixels.shape[0]),
                          int(left*pixels.shape[1]):int((left+width)*pixels.shape[1])]
            for color in ('#007d83', '#a66b0b', '#755297'):
                rgb = np.array([int(color[i:i+2], 16) for i in (1, 3, 5)]) / 255.
                assert np.sum(np.max(np.abs(crop-rgb), axis=2) < .04) > 5
    finally:
        plt.close(figure)
