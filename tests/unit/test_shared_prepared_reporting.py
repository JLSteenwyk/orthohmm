import copy
import json
from pathlib import Path

import matplotlib.pyplot as plt
import pytest

from benchmark_tools.results import plot_shared_prepared_resources_20261004 as plotter

RESULTS = Path(__file__).resolve().parents[2] / 'benchmark_tools/results'


def snapshot():
    return json.loads((RESULTS / 'threadripper_shared_panel_snapshot_20261004_v21/panel.json').read_text())


def test_both_actual_aborts_remain_missing_with_unchanged_measured_prefix():
    value = snapshot()
    prior = json.loads((RESULTS / 'threadripper_shared_panel_snapshot_20261004_v20/panel.json').read_text())
    assert value['runs'][:20] == prior['runs'][:20]
    assert (value['reviewed_attempts'], value['resource_reviewed_attempts'], value['eligible_attempts']) == (21, 19, 18)
    assert value['excluded_attempts'] == [0, 17, 20]
    assert value['pre_native_aborted_indices'] == [17, 20]
    assert len(plotter.validate(value)) == 19
    failure = json.loads((RESULTS / 'threadripper_shared_prenative_failure_22416.json').read_text())
    assert value['runs'][20]['job_id'] == failure['job_id'] == 22416
    assert value['runs'][20]['resources'] is None and failure['resources'] is None
    assert value['runs'][20]['status'] == 'reviewed_pre_native_abort'
    assert value['runs'][20]['whole_run_maximum_foreign_average_cores'] is None
    assert [row['index'] for row in value['runs'] if row['status'] == 'not_yet_reviewed'] == list(range(21, 27))
    for method in ('orthohmm_high_sensitivity', 'orthohmm_satellite_v2'):
        cell = next(cell for cell in value['cells'] if cell['method'] == method and cell['proteomes'] == 4)
        assert cell['reviewed_repeats'] == 3 and cell['eligible_repeats'] == 2
        assert cell['pending_indices'] == []
        assert all(summary is None for summary in cell['resources'].values())
    of4 = next(cell for cell in value['cells'] if cell['method'] == 'orthofinder_3_1_5_full' and cell['proteomes'] == 4)
    assert of4 == next(cell for cell in prior['cells'] if cell['method'] == of4['method'] and cell['proteomes'] == 4)


@pytest.mark.parametrize('change', ['eligible_abort', 'zero_endpoints', 'unattempted_abort',
    'resource_count', 'review_count', 'excluded_count', 'fabricated_summary'])
def test_reporting_rejects_rehabilitated_or_imputed_failure(change):
    value = copy.deepcopy(snapshot())
    if change == 'eligible_abort': value['runs'][20]['comparative_timing_eligible'] = True
    elif change == 'zero_endpoints':
        value['runs'][20]['resources'] = dict(wall_seconds=0, cpu_seconds=0, peak_memory_bytes=0)
    elif change == 'unattempted_abort': value['runs'][20]['status'] = 'not_yet_reviewed'
    elif change == 'resource_count': value['resource_reviewed_attempts'] = 20
    elif change == 'review_count': value['reviewed_attempts'] = 20
    elif change == 'excluded_count': value['excluded_attempts'] = [0, 17]
    else:
        cell = next(cell for cell in value['cells'] if cell['method'] == 'orthohmm_satellite_v2' and cell['proteomes'] == 4)
        cell['resources']['wall_seconds'] = dict(median=0, minimum=0, maximum=0)
    with pytest.raises(ValueError):
        plotter.validate(value)


def test_actual_two_abort_figure_bounds_pixels_and_source_identity():
    import fitz
    import numpy as np

    value = snapshot()
    directory = RESULTS / 'threadripper_shared_resource_figure_20261004_v20'
    manifest = json.loads((directory / 'manifest.json').read_text())
    assert manifest['source_results'] == plotter.executor.record(RESULTS / 'threadripper_shared_panel_snapshot_20261004_v21/panel.json')
    assert manifest['plotter'] == plotter.executor.record(plotter.__file__)
    assert manifest['pre_native_aborted_indices'] == [17, 20]
    assert manifest['excluded_indices'] == [0, 17, 20]
    for pin in manifest['outputs']:
        plotter.executor.check(pin)
    with fitz.open(directory / 'shared_threadripper_resources.pdf') as document:
        assert len(document) == 1
        page = document[0]
        text = page.get_text()
        assert '21/27 attempts reviewed' in text and '19 with measured resources' in text
        assert '2 pre-native aborts' in text and 'unknown and potentially method dependent' in text
        for block in page.get_text('dict')['blocks']:
            for line in block.get('lines', []):
                for span in line['spans']:
                    assert page.rect.contains(fitz.Rect(span['bbox']))
        pixmap = page.get_pixmap(alpha=False)
        pixels = np.frombuffer(pixmap.samples, dtype=np.uint8).reshape(pixmap.height, pixmap.width, 3) / 255.
    figure = plotter.plot(value)
    for axis in figure.axes:
        assert sum(len(collection.get_offsets()) for collection in axis.collections) == 19
        assert len(axis.lines) == 2
        left, bottom, width, height = axis.get_position().bounds
        crop = pixels[int((1-bottom-height)*pixels.shape[0]):int((1-bottom)*pixels.shape[0]),
            int(left*pixels.shape[1]):int((left+width)*pixels.shape[1])]
        for color in ('#007d83', '#a66b0b', '#755297'):
            rgb = np.array([int(color[index:index+2], 16) for index in (1, 3, 5)]) / 255.
            assert np.sum(np.max(np.abs(crop-rgb), axis=2) < .04) > 5
    plt.close(figure)
