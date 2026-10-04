import json
from pathlib import Path
from statistics import median

import matplotlib.pyplot as plt

from benchmark_tools.results import plot_shared_prepared_resources_20261004 as plotter

RESULTS = Path(__file__).resolve().parents[2] / 'benchmark_tools/results'
TABLE = RESULTS / 'threadripper_shared_panel_snapshot_20261004_v22/panel.json'
FIGURE = RESULTS / 'threadripper_shared_resource_figure_20261004_v21'


def test_actual_repeat21_preserves_prefix_and_completes_only_high8_cell():
    value = json.loads(TABLE.read_text())
    prior = json.loads((RESULTS / 'threadripper_shared_panel_snapshot_20261004_v21/panel.json').read_text())
    summary = json.loads((RESULTS / 'threadripper_shared_attempt_22417.json').read_text())
    assert value['runs'][:21] == prior['runs'][:21]
    assert (value['reviewed_attempts'], value['resource_reviewed_attempts'], value['eligible_attempts']) == (22, 20, 19)
    assert value['excluded_attempts'] == [0, 17, 20]
    assert value['pre_native_aborted_indices'] == [17, 20]
    assert len(plotter.validate(value)) == 20
    assert [row['index'] for row in value['runs'] if row['status'] == 'not_yet_reviewed'] == list(range(22, 27))
    row = value['runs'][21]
    assert row['job_id'] == summary['job_id'] == 22417
    assert row['resources'] == summary['resources']
    assert row['comparative_timing_eligible'] is True
    assert summary['native_outputs']['input_genes'] == 165168
    assert summary['native_outputs']['accuracy_evaluated'] is False
    complete = [cell for cell in value['cells'] if cell['eligible_repeats'] == 3]
    assert {(cell['method'], cell['proteomes']) for cell in complete} == {
        ('orthofinder_3_1_5_full', 4), ('orthohmm_high_sensitivity', 8)}
    high = next(cell for cell in complete if cell['proteomes'] == 8)
    assert high['reviewed_repeats'] == 3 and high['pending_indices'] == []
    for metric in plotter.tables.METRICS:
        observations = [value['runs'][index]['resources'][metric] for index in (5, 13, 21)]
        assert high['resources'][metric] == dict(median=median(observations), minimum=min(observations), maximum=max(observations))
    for cell in value['cells']:
        if cell['method'] != high['method'] or cell['proteomes'] != high['proteomes']:
            assert cell == next(old for old in prior['cells'] if (old['method'], old['proteomes']) == (cell['method'], cell['proteomes']))


def test_actual_terminal_handoff_and_four_review_pins_remain_bound():
    summary = json.loads((RESULTS / 'threadripper_shared_attempt_22417.json').read_text())
    executor = plotter.executor
    assert summary['original_environment_failures'] == dict(processes={}, pressure={})
    assert summary['shared_host_resources_reviewed'] is True
    for category, pin in summary['reviews'].items():
        review = executor.read(pin)
        assert review['decision'] == 'passed' and review['category'] == category
        assert review['index'] == 21 and review['job_id'] == 22417
        executor.check(review['source'])
    handoff = executor.read(summary['prepared_handoff'])
    assert handoff['worker_terminal_reviewed'] is True
    assert handoff['prepared_before_release_request'] is True
    assert handoff['request_to_bound_review_seconds'] < 20
    assert handoff['native_gate_go'] is True
    for pin in handoff['evidence']:
        executor.check(pin)


def test_actual_repeat21_figure_provenance_bounds_pixels_and_complete_cells():
    import fitz
    import numpy as np

    value = json.loads(TABLE.read_text())
    manifest = json.loads((FIGURE / 'manifest.json').read_text())
    assert manifest['source_results'] == plotter.executor.record(TABLE)
    assert manifest['plotter'] == plotter.executor.record(plotter.__file__)
    assert manifest['eligible_attempts'] == 19 and manifest['excluded_indices'] == [0, 17, 20]
    assert manifest['pre_native_aborted_indices'] == [17, 20]
    for pin in manifest['outputs']:
        plotter.executor.check(pin)
    with fitz.open(FIGURE / 'shared_threadripper_resources.pdf') as document:
        assert len(document) == 1
        page = document[0]
        text = page.get_text()
        assert '22/27 attempts reviewed' in text and '20 with measured resources' in text
        assert 'unknown and potentially method dependent' in text
        for block in page.get_text('dict')['blocks']:
            for line in block.get('lines', []):
                for span in line['spans']:
                    assert page.rect.contains(fitz.Rect(span['bbox']))
        pixmap = page.get_pixmap(alpha=False)
        pixels = np.frombuffer(pixmap.samples, dtype=np.uint8).reshape(pixmap.height, pixmap.width, 3) / 255.
    figure = plotter.plot(value)
    for axis in figure.axes:
        assert sum(len(collection.get_offsets()) for collection in axis.collections) == 20
        assert len(axis.lines) == 4
        left, bottom, width, height = axis.get_position().bounds
        crop = pixels[int((1-bottom-height)*pixels.shape[0]):int((1-bottom)*pixels.shape[0]),
            int(left*pixels.shape[1]):int((left+width)*pixels.shape[1])]
        for color in ('#007d83', '#a66b0b', '#755297'):
            rgb = np.array([int(color[index:index+2], 16) for index in (1, 3, 5)]) / 255.
            assert np.sum(np.max(np.abs(crop-rgb), axis=2) < .04) > 5
    plt.close(figure)
