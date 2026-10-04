from copy import deepcopy
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import pytest

from benchmark_tools.prepare_scaling_inputs import planned_runs
from benchmark_tools.results import report_shared_prenative_panel_20261004 as tables
from benchmark_tools.results import plot_shared_prenative_resources_20261004 as plotter


def reviewed(index, eligible=True):
    return dict(planned_runs()[index], job_id=100 + index,
        execution_scope=tables.executor.SHARED_SCOPE, uncontended_timing=False,
        scientific_timings_admitted=False, resource_scopes=tables.SCOPES,
        primary_resources_replayed=True, shared_host_resources_reviewed=eligible,
        original_environment_protocol_passed=eligible,
        status='shared_attempt_independently_reviewed' if eligible else 'shared_attempt_reviewed_with_environment_failure',
        scheduler_state='COMPLETED' if eligible else 'FAILED', scheduler_exit_code='0:0' if eligible else '1:0',
        resources=dict(wall_seconds=10. + index, cpu_seconds=100. + index, peak_memory_bytes=1024 + index),
        preflight_foreign_average_cores=50., whole_run_maximum_foreign_average_cores=60.)


def aborted(index=17):
    return dict(planned_runs()[index], job_id=100 + index,
        execution_scope=tables.executor.SHARED_SCOPE, uncontended_timing=False,
        scientific_timings_admitted=False, resources=None, comparative_timing_eligible=False,
        status='pre_native_infrastructure_failure_reviewed', native_outcome='not_started',
        automatic_retry=False, next_submission_authorized=False,
        scheduler_state='FAILED', scheduler_exit_code='1:0', preflight_foreign_average_cores=42.)


def report(count=18):
    return tables.summarize(planned_runs(), [reviewed(index, index != 0)
        for index in range(count) if index != 17], [aborted()] if count >= 18 else [])


def test_abort_retained_without_endpoint_or_false_complete_cell():
    value = report()
    assert value['reviewed_attempts'] == 18
    assert value['resource_reviewed_attempts'] == 17
    assert value['eligible_attempts'] == 16
    assert value['excluded_attempts'] == [0, 17]
    assert value['pre_native_aborted_indices'] == [17]
    assert value['runs'][17]['status'] == 'reviewed_pre_native_abort'
    assert value['runs'][17]['resources'] is None
    assert value['runs'][17]['whole_run_maximum_foreign_average_cores'] is None
    assert value['runs'][18]['status'] == 'not_yet_reviewed'
    assert all(metric is None for cell in value['cells'] for metric in cell['resources'].values())
    measured = plotter.validate(value)
    assert len(measured) == 17 and all(row['index'] != 17 for row in measured)


def test_complete_attempt_inventory_does_not_invent_missing_repeat():
    value = report(27)
    assert value['all_planned_attempts_reviewed']
    assert not value['all_cells_have_three_eligible_repeats']
    assert value['eligible_attempts'] == 25
    assert sum(cell['summary_status'] == 'three_eligible_repeats' for cell in value['cells']) == 7
    failed_cell = next(cell for cell in value['cells'] if cell['pre_native_aborted_indices'])
    assert failed_cell['eligible_repeats'] == 2 and failed_cell['resources']['wall_seconds'] is None
    assert len(plotter.validate(value)) == 26


@pytest.mark.parametrize('field,value', [
    ('native_outcome', 'exited_zero'), ('resources', dict(wall_seconds=0, cpu_seconds=0, peak_memory_bytes=0)),
    ('comparative_timing_eligible', True), ('automatic_retry', True), ('next_submission_authorized', True),
    ('scheduler_state', 'COMPLETED'), ('scheduler_exit_code', '0:0'), ('scientific_timings_admitted', True),
    ('uncontended_timing', True), ('index', True), ('job_id', 100), ('repeat', False),
    ('preflight_foreign_average_cores', float('nan')),
    ('whole_run_maximum_foreign_average_cores', float('inf'))])
def test_abort_cannot_borrow_success_or_impute_endpoints(field, value):
    abort = aborted()
    abort[field] = value
    with pytest.raises(ValueError):
        tables.summarize(planned_runs(), [reviewed(index, index != 0) for index in range(17)], [abort])


def test_prefix_duplicate_and_wrong_input_order_fail():
    with pytest.raises(ValueError):
        tables.summarize(planned_runs(), [reviewed(index) for index in range(16)], [aborted()])
    with pytest.raises(ValueError):
        tables.summarize(planned_runs(), [reviewed(index) for index in range(18)], [aborted()])
    changed = deepcopy(planned_runs())
    changed.reverse()
    with pytest.raises(ValueError): tables.summarize(changed, [], [])


@pytest.mark.parametrize('field,value', [('preflight_foreign_average_cores', True),
    ('whole_run_maximum_foreign_average_cores', None), ('preflight_foreign_average_cores', -1.)])
def test_measured_contention_requires_finite_original_annotations(field, value):
    row = reviewed(0)
    row[field] = value
    with pytest.raises(ValueError): tables.summarize(planned_runs(), [row], [])


@pytest.mark.parametrize('target,key,value', [
    ('abort', 'status', 'not_yet_reviewed'), ('abort', 'comparative_timing_eligible', True),
    ('abort', 'resources', dict(wall_seconds=0, cpu_seconds=0, peak_memory_bytes=0)),
    ('top', 'resource_reviewed_attempts', 18), ('top', 'reviewed_attempts', 17),
    ('top', 'excluded_attempts', [0]), ('top', 'pre_native_aborted_indices', []),
    ('top', 'publication_ready', True), ('top', 'uncontended_timing', True),
    ('pending', 'resources', dict(wall_seconds=0))])
def test_plot_recomputes_counts_rows_and_nulls(target, key, value):
    changed = report()
    selected = {'top': changed, 'abort': changed['runs'][17], 'pending': changed['runs'][18]}[target]
    selected[key] = value
    with pytest.raises(ValueError): plotter.validate(changed)


def test_export_distinguishes_abort_from_unreviewed_without_zero(tmp_path):
    output = tmp_path / 'table'
    tables.previous.export(output, report())
    with (output / 'attempts.tsv').open() as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    assert rows[17]['job_id'] == '117' and rows[17]['status'] == 'reviewed_pre_native_abort'
    assert rows[17]['comparative_timing_eligible'] == 'False' and rows[17]['wall_seconds'] == ''
    assert rows[18]['job_id'] == '' and rows[18]['status'] == 'not_yet_reviewed'
    with pytest.raises(FileExistsError): tables.previous.export(output, report())


@pytest.mark.parametrize('count', [0, 18, 27])
def test_no_abort_symbol_or_layout_clipping(count):
    figure = plotter.plot(report(count))
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    for text in figure.findobj(plt.Text):
        if text.get_visible() and text.get_text():
            bounds = text.get_window_extent(renderer)
            assert bounds.x0 >= 0 and bounds.y0 >= 0
            assert bounds.x1 <= figure.bbox.x1 and bounds.y1 <= figure.bbox.y1
    measured = max(count - (1 if count >= 18 else 0), 0)
    for axis in figure.axes:
        assert sum(len(collection.get_offsets()) for collection in axis.collections) == measured
    plt.close(figure)


@pytest.mark.parametrize('version,reviewed_count,resource_count,eligible_count', [
    ('v18', 18, 17, 16), ('v19', 19, 18, 17)])
def test_retained_actual_snapshot_and_pdf_pixels(tmp_path, version, reviewed_count, resource_count, eligible_count):
    import fitz
    import numpy as np

    results = Path(__file__).resolve().parents[2] / 'benchmark_tools/results'
    path = results / f'threadripper_shared_panel_snapshot_20261004_{version}/panel.json'
    value = json.loads(path.read_text())
    original = json.loads((results / 'threadripper_shared_panel_snapshot_20261003_v17/panel.json').read_text())
    assert value['runs'][:17] == original['runs'][:17]
    assert value['runs'][17]['job_id'] == 22413
    assert value['runs'][17]['resources'] is None
    assert len(plotter.validate(value)) == resource_count
    output = tmp_path / 'figure'
    manifest = plotter.export(path, plotter.executor.record(path)['sha256'], output)
    assert manifest['reviewed_attempts'] == reviewed_count and manifest['resource_reviewed_attempts'] == resource_count
    assert manifest['eligible_attempts'] == eligible_count and manifest['excluded_indices'] == [0, 17]
    assert manifest['pre_native_aborted_indices'] == [17]
    for pin in manifest['outputs']: plotter.executor.check(pin)
    with fitz.open(output / 'shared_threadripper_resources.pdf') as document:
        assert len(document) == 1
        page = document[0]
        text = page.get_text()
        assert f'{reviewed_count}/27 attempts reviewed' in text and f'{resource_count} with measured resources' in text
        assert '1 pre-native aborts' in text and 'not zero or unattempted runs' in text
        assert 'unknown and potentially method dependent' in text
        for block in page.get_text('dict')['blocks']:
            for line in block.get('lines', []):
                for span in line['spans']: assert page.rect.contains(fitz.Rect(span['bbox']))
        pixmap = page.get_pixmap(alpha=False)
        pixels = np.frombuffer(pixmap.samples, dtype=np.uint8).reshape(pixmap.height, pixmap.width, 3) / 255.
    figure = plotter.plot(value)
    for axis in figure.axes:
        left, bottom, width, height = axis.get_position().bounds
        crop = pixels[int((1-bottom-height)*pixels.shape[0]):int((1-bottom)*pixels.shape[0]),
            int(left*pixels.shape[1]):int((left+width)*pixels.shape[1])]
        for color in ('#007d83', '#a66b0b', '#755297'):
            rgb = np.array([int(color[index:index+2], 16) for index in (1, 3, 5)]) / 255.
            assert np.sum(np.max(np.abs(crop-rgb), axis=2) < .04) > 5
    plt.close(figure)


def test_first_complete_actual_cell_uses_all_three_reviewed_repeats():
    results = Path(__file__).resolve().parents[2] / 'benchmark_tools/results'
    value = json.loads((results / 'threadripper_shared_panel_snapshot_20261004_v19/panel.json').read_text())
    prior = json.loads((results / 'threadripper_shared_panel_snapshot_20261004_v18/panel.json').read_text())
    outcome = json.loads((results / 'threadripper_shared_attempt_22414.json').read_text())
    assert value['runs'][:18] == prior['runs'][:18]
    assert value['runs'][18]['job_id'] == outcome['job_id'] == 22414
    assert value['runs'][18]['resources'] == outcome['resources']
    assert outcome['native_outputs']['input_genes'] == 73266
    assert outcome['native_outputs']['checkpoint_groups'] == 24052
    assert outcome['native_outputs']['native_pair_rows'] == 88890
    assert outcome['native_outputs']['accuracy_evaluated'] is False
    assert (value['reviewed_attempts'], value['resource_reviewed_attempts'], value['eligible_attempts']) == (19, 18, 17)
    assert value['excluded_attempts'] == [0, 17] and value['pre_native_aborted_indices'] == [17]
    assert value['runs'][17]['resources'] is None and value['runs'][19]['resources'] is None
    cell = next(cell for cell in value['cells'] if cell['method'] == 'orthofinder_3_1_5_full' and cell['proteomes'] == 4)
    assert cell['eligible_repeats'] == 3 and cell['summary_status'] == 'three_eligible_repeats'
    assert cell['excluded_indices'] == [] and cell['pending_indices'] == []
    earliest = json.loads((results / 'threadripper_shared_panel_snapshot_20261003_v3/panel.json').read_text())
    summaries = [earliest['runs'][2],
        json.loads((results / 'threadripper_shared_attempt_22406.json').read_text()), outcome]
    assert [summary['job_id'] for summary in summaries] == [22398, 22406, 22414]
    for metric in tables.METRICS:
        ordered = sorted(summary['resources'][metric] for summary in summaries)
        assert cell['resources'][metric] == dict(median=ordered[1], minimum=ordered[0], maximum=ordered[2])
    assert sum(cell['summary_status'] == 'three_eligible_repeats' for cell in value['cells']) == 1
    figure = plotter.plot(value)
    for axis in figure.axes:
        assert sum(len(collection.get_offsets()) for collection in axis.collections) == 18
        assert len(axis.lines) == 2
    plt.close(figure)
