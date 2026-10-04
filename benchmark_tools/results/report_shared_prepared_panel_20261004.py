"""Report the frozen panel including pre-native aborts with unavailable endpoints."""

import argparse
import json
import math
from pathlib import Path
from statistics import median
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from benchmark_tools.results import report_shared_threadripper_panel_20261003 as previous
from benchmark_tools import run_threadripper_scaling as executor
from benchmark_tools.prepare_scaling_inputs import METHODS, SIZES, planned_runs

METRICS, IDENTITY, SCOPES = previous.METRICS, previous.IDENTITY, previous.SCOPES


def summarize(runs, reviews, aborts):
    expected = planned_runs()
    identities = [{key:run[key] for key in IDENTITY} for run in runs]
    if identities != expected or any(type(row[key]) is not type(want[key])
            for row,want in zip(identities, expected) for key in IDENTITY):
        raise ValueError('Changed frozen identities or order')
    by_index, jobs = {}, set()
    for review, is_abort in [(row, False) for row in reviews] + [(row, True) for row in aborts]:
        index, job = review['index'], review['job_id']
        if type(index) is not int or not 0 <= index < 27 or index in by_index or type(job) is not int or job <= 0 or job in jobs:
            raise ValueError('Duplicate or invalid reviewed identity/job')
        executor.expect(review, {key:runs[index][key] for key in IDENTITY})
        executor.expect(review, dict(execution_scope=executor.SHARED_SCOPE,
            uncontended_timing=False, scientific_timings_admitted=False))
        if is_abort:
            executor.expect(review, dict(status='pre_native_infrastructure_failure_reviewed',
                native_outcome='not_started', resources=None, comparative_timing_eligible=False,
                automatic_retry=False, next_submission_authorized=False,
                scheduler_state='FAILED', scheduler_exit_code='1:0'))
            status, eligible = 'reviewed_pre_native_abort', False
        else:
            eligible = review['shared_host_resources_reviewed']
            if type(eligible) is not bool or review['original_environment_protocol_passed'] is not eligible:
                raise ValueError('Inconsistent resource eligibility')
            executor.expect(review, dict(resource_scopes=SCOPES, primary_resources_replayed=True,
                status='shared_attempt_independently_reviewed' if eligible else 'shared_attempt_reviewed_with_environment_failure',
                scheduler_state='COMPLETED' if eligible else 'FAILED', scheduler_exit_code='0:0' if eligible else '1:0'))
            resources = review['resources']
            if set(resources) != set(METRICS) or any(type(v) not in (int,float) or not math.isfinite(v) or v <= 0 for v in resources.values()):
                raise ValueError('Invalid measured endpoints')
            if type(resources['peak_memory_bytes']) is not int or resources['peak_memory_bytes'] > 128 * 1024**3 or resources['wall_seconds'] > 85830:
                raise ValueError('Resource endpoint exceeds frozen envelope or has wrong units')
            status = 'reviewed_shared_observation' if eligible else 'reviewed_excluded_attempt'
        for key in ('preflight_foreign_average_cores', 'whole_run_maximum_foreign_average_cores'):
            value = review.get(key)
            if is_abort and key == 'whole_run_maximum_foreign_average_cores' and value is None:
                continue
            if type(value) not in (int, float) or not math.isfinite(value) or value < 0:
                raise ValueError('Invalid contention annotation')
        by_index[index] = dict(status=status, job_id=job, comparative_timing_eligible=eligible,
            resources=review['resources'], preflight_foreign_average_cores=review['preflight_foreign_average_cores'],
            whole_run_maximum_foreign_average_cores=review.get('whole_run_maximum_foreign_average_cores'))
        jobs.add(job)
    if set(by_index) != set(range(len(by_index))):
        raise ValueError('Reviewed prefix has a gap')
    rows = [dict(identity, **by_index.get(identity['index'], dict(status='not_yet_reviewed', job_id=None,
        comparative_timing_eligible=None, resources=None, preflight_foreign_average_cores=None,
        whole_run_maximum_foreign_average_cores=None))) for identity in identities]
    cells = []
    for method in METHODS:
        for size in SIZES:
            selected = [row for row in rows if (row['method'],row['proteomes']) == (method,size)]
            eligible = [row for row in selected if row['comparative_timing_eligible'] is True]
            excluded = [row['index'] for row in selected if row['comparative_timing_eligible'] is False]
            aborted = [row['index'] for row in selected if row['status'] == 'reviewed_pre_native_abort']
            complete = len(eligible) == 3
            cells.append(dict(method=method, proteomes=size, planned_repeats=3,
                reviewed_repeats=len(eligible)+len(excluded), eligible_repeats=len(eligible),
                excluded_indices=excluded, pre_native_aborted_indices=aborted,
                pending_indices=[row['index'] for row in selected if row['status'] == 'not_yet_reviewed'],
                summary_status='three_eligible_repeats' if complete else 'incomplete_eligible_repeats',
                resources={key:dict(median=median([row['resources'][key] for row in eligible]),
                    minimum=min(row['resources'][key] for row in eligible), maximum=max(row['resources'][key] for row in eligible))
                    if complete else None for key in METRICS}))
    return dict(status='shared_host_resource_panel_snapshot', execution_scope=executor.SHARED_SCOPE,
        planned_attempts=27, reviewed_attempts=len(by_index), resource_reviewed_attempts=len(reviews),
        eligible_attempts=sum(row['comparative_timing_eligible'] is True for row in rows),
        excluded_attempts=[row['index'] for row in rows if row['comparative_timing_eligible'] is False],
        pre_native_aborted_indices=[row['index'] for row in rows if row['status'] == 'reviewed_pre_native_abort'],
        all_planned_attempts_reviewed=len(by_index)==27,
        all_cells_have_three_eligible_repeats=all(cell['eligible_repeats']==3 for cell in cells),
        primary_scopes=SCOPES, runs=rows, cells=cells, uncontended_timing=False,
        scientific_timings_admitted=False, publication_ready=False,
        limitations=[*previous.summarize(runs, [])['limitations'],
            'Pre-native aborts are reviewed failed attempts with null endpoints, not zero timings or unattempted identities.',
            'The resource-reviewed count excludes pre-native aborts; exclusions retain both missing and measured failures.'])


def collect():
    results = ROOT / 'benchmark_tools/results'
    work = ROOT / 'benchmarks/work/threadripper_shared_execution_20261003'
    plan_path = results / 'threadripper_private_commands_20260928.json'
    plan_ref = executor.record(plan_path)
    plan = executor.read_frozen(plan_path, executor.PRIVATE_PLAN_SHA)
    paths = [results / 'threadripper_shared_attempt_22396.json']
    paths += [work / f'review_run_{index:02d}/summary.json' for index in range(1,27) if index not in (17, 20)]
    present = [path for path in paths if path.exists()]
    refs, reviews, sessions = [plan_ref], [], []
    for path in present:
        pin = executor.record(path)
        review = executor.read(pin)
        session_ref = review['session']
        session = executor.read(session_ref)
        controller_ref = session['controller']
        previous.terminal_matches(review, session, executor.read(controller_ref), plan_ref)
        reviews.append(review)
        refs += [pin, session_ref, controller_ref]
        sessions.append(session_ref)
        for category in previous.REVIEW_CATEGORIES:
            category_ref = review['reviews'][category]
            value = executor.read(category_ref)
            executor.expect(value, dict(schema='threadripper_panel_review_v1', index=review['index'], job_id=review['job_id'],
                plan_sha256=plan_ref['sha256'], category=category,
                decision='failed' if category=='environment' and not review['shared_host_resources_reviewed'] else 'passed',
                execution_scope=executor.SHARED_SCOPE, uncontended_timing=False))
            for supporting in [value['source'], *value['evidence']]:
                executor.check(supporting)
            refs += [category_ref, value['source']]
            if category == 'resources':
                matches = [pin for pin in value['evidence'] if Path(pin['path']).name == 'resources.json']
                if len(matches) != 1:
                    raise ValueError('Missing independently derived endpoints')
                resource = executor.read(matches[0])
                if resource['primary'] != review['resources'] or resource['primary_scopes'] != review['resource_scopes']:
                    raise ValueError('Resource summary differs from raw derivation')
                refs += matches
            if category == 'environment':
                matches = [pin for pin in value['evidence'] if Path(pin['path']).name == 'environment_replay.json']
                if len(matches) != 1:
                    raise ValueError('Missing independent environment replay')
                replay = executor.read(matches[0])
                previous.environment_matches(review, replay, executor.read(replay['preflight']))
                refs += [matches[0], replay['preflight']]
    aborts, abort_refs = [], []
    for index, job in ((17, 22413), (20, 22416)):
        abort_ref = executor.record(results / f'threadripper_shared_prenative_failure_{job}.json')
        resolved_ref = executor.record(work / f'resolved_session_{index:02d}_preparation_sync_20261004.json')
        if executor.read(resolved_ref)['pre_native_audit'] != abort_ref:
            raise ValueError('Resolved history borrows another abort')
        aborts.append(executor.read(abort_ref))
        abort_refs.extend([abort_ref, resolved_ref])
        sessions.insert(index, resolved_ref)
    # Index 0 retains its explicit cadence resolution rather than its original failed session.
    sessions[0] = executor.record(work / 'resolved_session_00.json')
    bound = executor.bind(plan_ref, sessions, command=str(ROOT / 'benchmark_tools/run_threadripper_shared_scaling.sh'),
        cwd=str(ROOT), time_limit=executor.TIME_LIMIT, allocation_mode='shared')
    if bound['progress']['status'] not in {'next_identity_requires_preflight','all_attempts_reviewed'}:
        raise ValueError('Report contains an unresolved history prefix')
    report = summarize(plan['runs'], reviews, aborts)
    refs += [*abort_refs, *bound['evidence']]
    for pin in refs:
        executor.check(pin)
    if present != [path for path in paths if path.exists()]:
        raise ValueError('Review inventory changed during snapshot')
    report.update(plan=plan_ref, evidence=refs, source=executor.record(__file__),
        prior_reporter=executor.record(previous.__file__))
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = collect()
    outputs = previous.export(args.output.resolve(), report)
    print(json.dumps(dict(reviewed=report['reviewed_attempts'], resource_reviewed=report['resource_reviewed_attempts'],
        eligible=report['eligible_attempts'], excluded=report['excluded_attempts'], outputs=outputs), indent=2))
