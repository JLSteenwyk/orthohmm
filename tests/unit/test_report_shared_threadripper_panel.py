from copy import deepcopy
import csv
import json
from pathlib import Path

import pytest

from benchmark_tools.results.report_shared_threadripper_panel_20261003 import SCOPES, export, summarize
from benchmark_tools.results import report_shared_threadripper_panel_20261003 as reporter
from benchmark_tools.prepare_scaling_inputs import planned_runs


def reviewed(index, eligible=True):
    return dict(planned_runs()[index], job_id=100 + index, execution_scope='shared_host_matched_resources',
        uncontended_timing=False, scientific_timings_admitted=False, resource_scopes=SCOPES,
        primary_resources_replayed=True, shared_host_resources_reviewed=eligible,
        original_environment_protocol_passed=eligible,
        status='shared_attempt_independently_reviewed' if eligible else 'shared_attempt_reviewed_with_environment_failure',
        scheduler_state='COMPLETED' if eligible else 'FAILED', scheduler_exit_code='0:0' if eligible else '1:0',
        resources=dict(wall_seconds=10. + index, cpu_seconds=100. + index, peak_memory_bytes=1024 + index),
        preflight_foreign_average_cores=50., whole_run_maximum_foreign_average_cores=60.)


def test_empty_panel_has_nulls_not_estimates_or_completed_claims():
    result = summarize(planned_runs(), [])
    assert len(result['runs']) == 27 and len(result['cells']) == 9
    assert result['reviewed_attempts'] == 0
    assert not result['all_planned_attempts_reviewed'] and not result['all_cells_have_three_eligible_repeats']
    assert all(row['resources'] is None for row in result['runs'])
    assert all(value is None for cell in result['cells'] for value in cell['resources'].values())


def test_excluded_attempt_keeps_raw_values_but_never_supplies_complete_cell():
    result = summarize(planned_runs(), [reviewed(0, False)])
    assert result['runs'][0]['resources']['wall_seconds'] == 10.
    assert result['runs'][0]['comparative_timing_eligible'] is False
    assert result['excluded_attempts'] == [0] and result['eligible_attempts'] == 0
    assert result['cells'][0]['excluded_indices'] == [0]
    assert result['cells'][0]['resources']['wall_seconds'] is None


def test_complete_panel_with_one_exclusion_has_eight_complete_cells_not_nine():
    rows = [reviewed(index, index != 0) for index in range(27)]
    report = summarize(planned_runs(), rows)
    assert report['all_planned_attempts_reviewed'] and not report['all_cells_have_three_eligible_repeats']
    assert report['eligible_attempts'] == 26
    assert report['cells'][0]['eligible_repeats'] == 2
    assert report['cells'][0]['resources']['wall_seconds'] is None
    assert sum(cell['summary_status'] == 'three_eligible_repeats' for cell in report['cells']) == 8


def test_complete_eligible_summaries_use_all_repeats_and_preserve_maximum():
    report = summarize(planned_runs(), [reviewed(index) for index in range(27)])
    cell = report['cells'][0]
    assert cell['resources']['wall_seconds'] == dict(median=21., minimum=10., maximum=29.)
    assert report['all_cells_have_three_eligible_repeats']
    assert not report['scientific_timings_admitted'] and not report['publication_ready']


@pytest.mark.parametrize('field,value', [('index', True), ('method', 'historical_dgx'), ('job_id', False),
    ('execution_scope', 'isolated_controlled'), ('uncontended_timing', True),
    ('scientific_timings_admitted', True), ('primary_resources_replayed', False),
    ('original_environment_protocol_passed', False), ('scheduler_state', 'RUNNING'),
    ('scheduler_exit_code', '1:0'), ('shared_host_resources_reviewed', 1)])
def test_identity_isolation_live_or_admission_drift_rejected(field, value):
    row = reviewed(0)
    row[field] = value
    with pytest.raises(ValueError):
        summarize(planned_runs(), [row])


@pytest.mark.parametrize('field,value', [('wall_seconds', float('nan')), ('cpu_seconds', float('inf')),
    ('cpu_seconds', True), ('wall_seconds', 0), ('peak_memory_bytes', '1024'),
    ('peak_memory_bytes', 1024.), ('peak_memory_bytes', 129 * 1024**3), ('wall_seconds', 85831)])
def test_invalid_values_never_become_figures_or_cells(field, value):
    row = reviewed(0)
    row['resources'][field] = value
    with pytest.raises(ValueError):
        summarize(planned_runs(), [row])


def test_scope_and_missing_prefix_are_rejected():
    row = reviewed(0)
    row['resource_scopes'] = dict(SCOPES, peak_memory_bytes='process_rss')
    with pytest.raises(ValueError): summarize(planned_runs(), [row])
    with pytest.raises(ValueError): summarize(planned_runs(), [reviewed(1)])
    with pytest.raises(ValueError): summarize(planned_runs(), [reviewed(0), reviewed(0)])
    changed = deepcopy(planned_runs())
    changed[0]['repeat'] = False
    with pytest.raises(ValueError): summarize(changed, [])


def test_export_preserves_exclusion_nulls_and_refuses_overwrite(tmp_path):
    output = tmp_path / 'panel'
    report = summarize(planned_runs(), [reviewed(0, False)])
    artifacts = export(output, report)
    assert len(artifacts) == 3
    assert json.loads((output / 'panel.json').read_text()) == report
    with (output / 'attempts.tsv').open() as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    assert rows[0]['comparative_timing_eligible'] == 'False' and rows[0]['wall_seconds'] == '10.0'
    assert rows[1]['wall_seconds'] == '' and rows[1]['comparative_timing_eligible'] == ''
    with (output / 'cells.tsv').open() as stream:
        cells = list(csv.DictReader(stream, delimiter='\t'))
    assert cells[0]['wall_seconds_median'] == '' and cells[0]['eligible_repeats'] == '0'
    with pytest.raises(FileExistsError): export(output, report)


def evidence_fixture(eligible=True):
    row = reviewed(0, eligible)
    row['reviews'] = dict.fromkeys(reporter.REVIEW_CATEGORIES, 'fixture')
    plan = dict(sha256='fixture_plan')
    session = dict(schema='threadripper_panel_session_v1', index=0, job_id=row['job_id'],
        plan_sha256=plan['sha256'], phase='terminal', native_outcome='exited_zero', reviews=row['reviews'])
    raw = (f"JobId={row['job_id']} JobState={row['scheduler_state']} Partition=gpu NodeList=bizon NumNodes=1 "
        "NumCPUs=64 NumTasks=1 CPUs/Task=64 OverSubscribe=OK MinMemoryNode=128G Requeue=0 Restarts=0 "
        f"Command={reporter.ROOT}/benchmark_tools/run_threadripper_shared_scaling.sh WorkDir={reporter.ROOT} "
        f"TimeLimit={reporter.executor.TIME_LIMIT} ExitCode={row['scheduler_exit_code']}")
    controller = dict(command=['scontrol', 'show', 'job', str(row['job_id']), '--oneliner'], returncode=0, stdout=raw)
    preflight = dict(index=0, job_id=row['job_id'], execution_scope=reporter.executor.SHARED_SCOPE,
        decision='passed', uncontended_timing=False, observed_foreign_average_cores=50.)
    replay = dict(execution_scope=reporter.executor.SHARED_SCOPE, uncontended_timing=False, scientific_timings_admitted=False,
        processes=dict(sampled_process_policy_satisfied=eligible, maximum_observed_foreign_average_cores=60.),
        pressure=dict(sampled_pressure_evidence_satisfied=True))
    return row, session, controller, plan, replay, preflight


@pytest.mark.parametrize('eligible', [False, True])
def test_terminal_and_environment_metadata_follow_original_evidence(eligible):
    row, session, controller, plan, replay, preflight = evidence_fixture(eligible)
    reporter.terminal_matches(row, session, controller, plan)
    reporter.environment_matches(row, replay, preflight)


@pytest.mark.parametrize('target,key,value', [
    ('session', 'phase', 'running'), ('session', 'job_id', 999),
    ('controller', 'returncode', 1), ('controller', 'command', ['sacct']),
    ('row', 'scheduler_state', 'FAILED'), ('preflight', 'observed_foreign_average_cores', 0.),
    ('row', 'whole_run_maximum_foreign_average_cores', 0.), ('row', 'shared_host_resources_reviewed', False),
    ('replay', 'uncontended_timing', True)])
def test_summary_cannot_relabel_retained_terminal_or_contention(target, key, value):
    row, session, controller, plan, replay, preflight = evidence_fixture()
    dict(row=row, session=session, controller=controller, replay=replay, preflight=preflight)[target][key] = value
    with pytest.raises(ValueError):
        reporter.terminal_matches(row, session, controller, plan)
        reporter.environment_matches(row, replay, preflight)


@pytest.mark.parametrize('version,count', [('v4', 4), ('v5', 5), ('v6', 6), ('v7', 7)])
def test_retained_partial_snapshots_replay_tables_without_original_evidence(tmp_path, version, count):
    results = Path(__file__).resolve().parents[2] / 'benchmark_tools/results'
    snapshot = results / f'threadripper_shared_panel_snapshot_20261003_{version}'
    report = json.loads((snapshot / 'panel.json').read_text())
    outcome = json.loads((results / 'threadripper_shared_attempt_22399.json').read_text())
    row = report['runs'][3]
    assert {key: row[key] for key in reporter.IDENTITY} == planned_runs()[3]
    assert row['job_id'] == outcome['job_id'] == 22399
    assert row['resources'] == outcome['resources'] == dict(
        wall_seconds=1999.550829241, cpu_seconds=56008.34709, peak_memory_bytes=8243466240)
    assert row['whole_run_maximum_foreign_average_cores'] == outcome['whole_run_maximum_foreign_average_cores']
    assert report['reviewed_attempts'] == count and report['eligible_attempts'] == count - 1
    assert report['excluded_attempts'] == [0]
    assert all(row['resources'] is None for row in report['runs'][count:])
    if count >= 5:
        orthofinder = json.loads((results / 'threadripper_shared_attempt_22400.json').read_text())
        assert report['runs'][4]['job_id'] == orthofinder['job_id'] == 22400
        assert report['runs'][4]['resources'] == orthofinder['resources'] == dict(
            wall_seconds=1139.939834238, cpu_seconds=17734.408768, peak_memory_bytes=10987712512)
        assert orthofinder['native_outputs']['input_genes'] == 165168
        assert orthofinder['native_outputs']['checkpoint_groups'] == 33013
        assert orthofinder['native_outputs']['native_pair_rows'] == 597451
    if count >= 6:
        high = json.loads((results / 'threadripper_shared_attempt_22401.json').read_text())
        assert report['runs'][5]['job_id'] == high['job_id'] == 22401
        assert report['runs'][5]['resources'] == high['resources'] == dict(
            wall_seconds=1238.935206005, cpu_seconds=34395.748088, peak_memory_bytes=8293474304)
        assert high['native_outputs']['input_genes'] == 165168
        assert high['native_outputs']['orthogroups'] == 58278
        assert all(cell['eligible_repeats'] == 1 for cell in report['cells'] if cell['proteomes'] == 8)
    if count == 7:
        full = json.loads((results / 'threadripper_shared_attempt_22402.json').read_text())
        assert report['runs'][6]['job_id'] == full['job_id'] == 22402
        assert report['runs'][6]['resources'] == full['resources'] == dict(
            wall_seconds=1896.732060134, cpu_seconds=39246.594856, peak_memory_bytes=11897724928)
        assert full['native_outputs']['input_genes'] == 251378
        assert full['native_outputs']['checkpoint_groups'] == 34230
        assert full['native_outputs']['native_pair_rows'] == 1487084
        assert sum(cell['eligible_repeats'] for cell in report['cells'] if cell['proteomes'] == 12) == 1
    assert all(value is None for cell in report['cells'] for value in cell['resources'].values())
    assert not report['all_planned_attempts_reviewed'] and not report['publication_ready']
    output = tmp_path / 'partial_attempt_replay'
    export(output, report)
    assert all((output / name).read_bytes() == (snapshot / name).read_bytes()
               for name in ('panel.json', 'attempts.tsv', 'cells.tsv'))
