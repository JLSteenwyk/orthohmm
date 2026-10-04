from copy import deepcopy
import json
from pathlib import Path

import pytest

from benchmark_tools.results import render_shared_resource_section_20261004 as section

RESULTS = Path(__file__).resolve().parents[2] / 'benchmark_tools/results'


@pytest.fixture
def report():
    return json.loads((RESULTS / 'threadripper_shared_panel_snapshot_20261004_v24/panel.json').read_text())


def test_actual_partial_section_keeps_missing_cells_and_units(report):
    text = section.render(report)
    assert '(Interim)' in text and '24 of 27 reviewed attempts' in text
    assert '22 with measured native resources and 21 eligible observations' in text
    assert '4 method/size cells' in text
    assert text.count('Unavailable') == 15
    assert '1238.9352 [1221.4444, 1303.6330]' in text
    assert '34395.7481 [33878.9946, 36610.5696]' in text
    value = report['cells'][1]['resources']['peak_memory_bytes']
    assert f"{value['median']/1024**3:.3f}" in text
    assert 'Excluded attempt indices: 0, 17, 20' in text
    assert 'Unreviewed indices: 24, 25, 26' in text
    assert 'not zero measurements' in text and 'not isolated' in text
    assert 'unknown and potentially method dependent' in text
    assert 'not taxon-invariant scaling' in text


@pytest.mark.parametrize('key', ['reviewed_attempts', 'eligible_attempts', 'resource_reviewed_attempts'])
def test_corrupt_counts_are_not_reported(report, key):
    report[key] += 1
    with pytest.raises(ValueError): section.render(report)


def test_corrupt_summary_is_not_reported(report):
    report['cells'][1]['resources']['wall_seconds']['median'] += .1
    with pytest.raises(ValueError): section.render(report)


def test_completed_synthetic_panel_does_not_impute_failed_cells(report):
    # Fill future rows only as a test fixture; no synthetic row is persisted as evidence.
    changed = deepcopy(report)
    for row in changed['runs'][24:]:
        old = next(r for r in changed['runs'][:24]
            if (r['method'], r['proteomes']) == (row['method'], row['proteomes']))
        row.update({k: deepcopy(old[k]) for k in ('status', 'resources',
            'comparative_timing_eligible', 'preflight_foreign_average_cores',
            'whole_run_maximum_foreign_average_cores')}, job_id=100000+row['index'])
    reviews, aborts = [], []
    for row in changed['runs']:
        identity = {k: row[k] for k in section.plotter.tables.IDENTITY}
        common = dict(identity, job_id=row['job_id'], execution_scope=section.executor.SHARED_SCOPE,
            uncontended_timing=False, scientific_timings_admitted=False, resources=row['resources'],
            preflight_foreign_average_cores=row['preflight_foreign_average_cores'],
            whole_run_maximum_foreign_average_cores=row['whole_run_maximum_foreign_average_cores'])
        if row['status'] == 'reviewed_pre_native_abort':
            aborts.append(dict(common, status='pre_native_infrastructure_failure_reviewed',
                native_outcome='not_started', comparative_timing_eligible=False, automatic_retry=False,
                next_submission_authorized=False, scheduler_state='FAILED', scheduler_exit_code='1:0'))
        else:
            eligible = row['comparative_timing_eligible']
            reviews.append(dict(common, resource_scopes=section.plotter.tables.SCOPES,
                primary_resources_replayed=True, shared_host_resources_reviewed=eligible,
                original_environment_protocol_passed=eligible,
                status='shared_attempt_independently_reviewed' if eligible else 'shared_attempt_reviewed_with_environment_failure',
                scheduler_state='COMPLETED' if eligible else 'FAILED', scheduler_exit_code='0:0' if eligible else '1:0'))
    complete = section.plotter.tables.summarize(section.plotter.planned_runs(), reviews, aborts)
    text = section.render(complete)
    assert '(Interim)' not in text and '27 of 27 reviewed attempts' in text
    assert '6 method/size cells' in text and text.count('Unavailable') == 9
    assert 'Unreviewed indices: none' in text
    assert 'Excluded attempt indices: 0, 17, 20' in text


def test_no_overwrite_or_bad_pin_creates_output(tmp_path, report):
    existing = tmp_path / 'existing'
    existing.mkdir()
    with pytest.raises(FileExistsError): section.build(tmp_path / 'missing', '0'*64, existing)
    table = tmp_path / 'table.json'
    table.write_text(json.dumps(report))
    fresh = tmp_path / 'fresh'
    with pytest.raises(ValueError): section.build(table, '0'*64, fresh)
    assert not fresh.exists()


def test_real_build_binds_section_and_table(tmp_path, report):
    table = RESULTS / 'threadripper_shared_panel_snapshot_20261004_v24/panel.json'
    pin = section.executor.record(table)
    receipt = section.build(table, pin['sha256'], tmp_path / 'section')
    assert receipt['table'] == pin
    assert receipt['native_inference_repeated'] is False and receipt['raw_audit_repeated'] is False
    section.executor.check(receipt['section'])
    assert Path(receipt['section']['path']).read_text() == section.render(report)
