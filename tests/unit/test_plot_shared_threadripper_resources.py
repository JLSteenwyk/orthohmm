from copy import deepcopy
import json
from pathlib import Path

import matplotlib.pyplot as plt
import pytest

from benchmark_tools.results import plot_shared_threadripper_resources_20261003 as module
from benchmark_tools.results.report_shared_threadripper_panel_20261003 import summarize
from benchmark_tools.prepare_scaling_inputs import planned_runs


def fixture(index, eligible=True):
    return dict(planned_runs()[index], job_id=100 + index, execution_scope='shared_host_matched_resources',
        uncontended_timing=False, scientific_timings_admitted=False, resource_scopes=module.tables.SCOPES,
        primary_resources_replayed=True, shared_host_resources_reviewed=eligible,
        original_environment_protocol_passed=eligible,
        status='shared_attempt_independently_reviewed' if eligible else 'shared_attempt_reviewed_with_environment_failure',
        scheduler_state='COMPLETED' if eligible else 'FAILED', scheduler_exit_code='0:0' if eligible else '1:0',
        resources=dict(wall_seconds=60. + index, cpu_seconds=3600. + index, peak_memory_bytes=1024**3 + index),
        preflight_foreign_average_cores=50., whole_run_maximum_foreign_average_cores=60.)


def test_empty_figure_never_creates_zero_observations_or_summaries():
    report = summarize(planned_runs(), [])
    figure = module.plot(report)
    assert len(figure.axes) == 3
    assert all(not axis.collections and not axis.lines for axis in figure.axes)
    assert '(partial)' in ' '.join(text.get_text() for text in figure.texts)
    plt.close(figure)


def test_partial_figure_keeps_exclusion_and_has_no_incomplete_medians():
    report = summarize(planned_runs(), [fixture(0, False), fixture(1)])
    figure = module.plot(report)
    for axis in figure.axes:
        assert sum(len(points.get_offsets()) for points in axis.collections) == 2
        assert len(axis.lines) == 0
        assert axis.get_ylim()[0] == 0
        assert axis.get_ylim()[1] > max(points.get_offsets()[0, 1] for points in axis.collections)
    assert list(figure.axes[0].collections[0].get_offsets()[0]) == pytest.approx([-.26, 1.])
    assert list(figure.axes[1].collections[0].get_offsets()[0]) == pytest.approx([-.26, 1.])
    assert list(figure.axes[2].collections[0].get_offsets()[0]) == pytest.approx([-.26, 1.])
    text = ' '.join(item.get_text() for item in figure.texts)
    assert '1 eligible | 1 excluded' in text
    assert 'not confidence intervals' in text and 'potentially method dependent' in text
    assert '4 proteomes: 0/3 eligible, 1 excluded' in text
    plt.close(figure)


def test_complete_panel_retains_every_attempt_and_only_eight_complete_cells():
    report = summarize(planned_runs(), [fixture(index, index != 0) for index in range(27)])
    figure = module.plot(report)
    for axis in figure.axes:
        assert sum(len(points.get_offsets()) for points in axis.collections) == 27
        assert len(axis.lines) == 16
    assert '(partial)' not in ' '.join(text.get_text() for text in figure.texts)
    assert not report['all_cells_have_three_eligible_repeats']
    plt.close(figure)


@pytest.mark.parametrize('problem', ['identity', 'missing', 'impute', 'excluded', 'count', 'scope', 'isolation',
    'admission', 'median', 'contended', 'nan', 'ready'])
def test_modified_table_or_unsupported_claim_is_refused(problem):
    report = deepcopy(summarize(planned_runs(), [fixture(0, False), fixture(1)]))
    if problem == 'identity': report['runs'][1]['proteomes'] = 8
    elif problem == 'missing': report['runs'].pop()
    elif problem == 'impute': report['runs'][2]['resources'] = dict(wall_seconds=0.)
    elif problem == 'excluded': report['runs'][0]['status'] = 'reviewed_shared_observation'
    elif problem == 'count': report['eligible_attempts'] = 2
    elif problem == 'scope': report['primary_scopes']['cpu_seconds'] = 'process_rss'
    elif problem == 'isolation': report['uncontended_timing'] = True
    elif problem == 'admission': report['scientific_timings_admitted'] = True
    elif problem == 'median': report['cells'][0]['resources']['wall_seconds'] = dict(median=1.)
    elif problem == 'contended': report['runs'][0]['preflight_foreign_average_cores'] = -1.
    elif problem == 'nan': report['runs'][0]['whole_run_maximum_foreign_average_cores'] = float('nan')
    else: report['publication_ready'] = True
    with pytest.raises(ValueError): module.validate(report)


def test_export_refuses_missing_evidence_and_existing_output(tmp_path):
    report = deepcopy(summarize(planned_runs(), []))
    report['source'] = dict(path=str(tmp_path / 'absent'), bytes=1, sha256='absent')
    report['plan'] = report['source']
    report['evidence'] = []
    path = tmp_path / 'report.json'
    path.write_text(json.dumps(report))
    with pytest.raises(FileNotFoundError):
        module.export(path, module.executor.record(path)['sha256'], tmp_path / 'figure')
    assert not (tmp_path / 'figure').exists()
    with pytest.raises(FileExistsError): module.export(path, 'unused', tmp_path)


def test_retained_actual_table_is_consistent_with_plot():
    path = Path(__file__).resolve().parents[2] / 'benchmark_tools/results/threadripper_shared_panel_snapshot_20261003_v2/panel.json'
    report = json.loads(path.read_text())
    rows = module.validate(report)
    assert [row['job_id'] for row in rows] == [22396, 22397]
    assert [row['comparative_timing_eligible'] for row in rows] == [False, True]


@pytest.mark.parametrize('count', [2, 27])
def test_figure_text_fits_for_partial_and_full_panels(count):
    report = summarize(planned_runs(), [fixture(index, index != 0) for index in range(count)])
    figure = module.plot(report)
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    for text in figure.findobj(plt.Text):
        if not text.get_visible() or not text.get_text():
            continue
        bounds = text.get_window_extent(renderer)
        assert bounds.x0 >= 0 and bounds.y0 >= 0
        assert bounds.x1 <= figure.bbox.x1 and bounds.y1 <= figure.bbox.y1
    plt.close(figure)


@pytest.mark.parametrize('version,reviewed,eligible', [('v3', 3, 2), ('v4', 4, 3), ('v5', 5, 4), ('v6', 6, 5), ('v7', 7, 6), ('v8', 8, 7)])
def test_actual_export_pdf_text_pixels_and_provenance(tmp_path, version, reviewed, eligible):
    import fitz
    import numpy as np

    path = Path(__file__).resolve().parents[2] / f'benchmark_tools/results/threadripper_shared_panel_snapshot_20261003_{version}/panel.json'
    report = json.loads(path.read_text())
    output = tmp_path / 'figure'
    manifest = module.export(path, module.executor.record(path)['sha256'], output)
    assert manifest['reviewed_attempts'] == reviewed and manifest['eligible_attempts'] == eligible
    assert manifest['excluded_indices'] == [0]
    assert not manifest['uncontended_timing'] and not manifest['publication_ready']
    for pin in manifest['outputs']:
        module.executor.check(pin)
    assert json.loads((output / 'manifest.json').read_text()) == manifest
    with fitz.open(output / 'shared_threadripper_resources.pdf') as document:
        assert len(document) == 1
        page = document[0]
        text = page.get_text()
        assert '(partial)' in text and f'{reviewed}/27 attempts reviewed' in text
        assert 'unknown and potentially method dependent' in text
        assert 'excluded raw values' in text
        for block in page.get_text('dict')['blocks']:
            for line in block.get('lines', []):
                for span in line['spans']:
                    assert page.rect.contains(fitz.Rect(span['bbox']))
        pixmap = page.get_pixmap(alpha=False)
        pixels = np.frombuffer(pixmap.samples, dtype=np.uint8).reshape(pixmap.height, pixmap.width, 3) / 255.
    figure = module.plot(report)
    for axis in figure.axes:
        left, bottom, width, height = axis.get_position().bounds
        crop = pixels[int((1 - bottom - height) * pixels.shape[0]):int((1 - bottom) * pixels.shape[0]),
            int(left * pixels.shape[1]):int((left + width) * pixels.shape[1])]
        for color in ('#a66b0b', '#755297'):
            rgb = np.array([int(color[index:index + 2], 16) for index in (1, 3, 5)]) / 255.
            assert np.sum(np.max(np.abs(crop - rgb), axis=2) < .04) > 5
    plt.close(figure)
