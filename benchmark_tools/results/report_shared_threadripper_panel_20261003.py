"""Export reviewed shared-host resources, preserving excluded and missing repeats."""
import argparse
import csv
import json
import math
from pathlib import Path
from statistics import median
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from benchmark_tools.prepare_scaling_inputs import METHODS, SIZES, planned_runs
from benchmark_tools.derive_threadripper_resources import SCOPES
from benchmark_tools import run_threadripper_scaling as executor
from benchmark_tools.verify_threadripper_controller import validate

METRICS = ('wall_seconds', 'cpu_seconds', 'peak_memory_bytes')
IDENTITY = ('index', 'method', 'proteomes', 'repeat')
REVIEW_CATEGORIES = ('runtime', 'environment', 'resources', 'outputs_or_failure')


def terminal_matches(review, session, controller, plan_ref):
    executor.expect(session, dict(schema='threadripper_panel_session_v1',
        index=review['index'], job_id=review['job_id'], plan_sha256=plan_ref['sha256'],
        phase='terminal', native_outcome='exited_zero', reviews=review['reviews']))
    executor.expect(controller, dict(command=['scontrol', 'show', 'job', str(review['job_id']), '--oneliner'],
        returncode=0))
    allocation = validate(controller['stdout'], review['job_id'], 'terminal',
        command=str(ROOT / 'benchmark_tools/run_threadripper_shared_scaling.sh'),
        cwd=str(ROOT), time_limit=executor.TIME_LIMIT, allocation_mode='shared')
    if (allocation['scheduler_state'], allocation['scheduler_exit_code']) != (
            review['scheduler_state'], review['scheduler_exit_code']):
        raise ValueError('Published terminal outcome differs from retained controller')


def environment_matches(review, replay, preflight):
    executor.expect(replay, dict(execution_scope=executor.SHARED_SCOPE,
        uncontended_timing=False, scientific_timings_admitted=False))
    executor.expect(preflight, dict(index=review['index'], job_id=review['job_id'],
        execution_scope=executor.SHARED_SCOPE, decision='passed', uncontended_timing=False))
    processes, pressure = replay['processes'], replay['pressure']
    decisions = (processes['sampled_process_policy_satisfied'], pressure['sampled_pressure_evidence_satisfied'])
    if any(type(value) is not bool for value in decisions) or all(decisions) is not review['shared_host_resources_reviewed']:
        raise ValueError('Published eligibility differs from independent environment replay')
    if (preflight['observed_foreign_average_cores'] != review['preflight_foreign_average_cores']
            or processes['maximum_observed_foreign_average_cores'] != review['whole_run_maximum_foreign_average_cores']):
        raise ValueError('Published contention differs from retained observations')


def summarize(runs, reviews):
    identities = [{key: run[key] for key in IDENTITY} for run in runs]
    expected = planned_runs()
    if identities != expected or any(type(row[key]) is not type(want[key])
            for row, want in zip(identities, expected) for key in IDENTITY):
        raise ValueError('Changed frozen panel identities or order')
    by_index, jobs = {}, set()
    for review in reviews:
        index, job = review['index'], review['job_id']
        if type(index) is not int or not 0 <= index < 27 or index in by_index or type(job) is not int or job <= 0 or job in jobs:
            raise ValueError('Duplicate or invalid reviewed identity/job')
        if any(type(review.get(key)) is not type(runs[index][key]) or review[key] != runs[index][key] for key in IDENTITY):
            raise ValueError('Review belongs to a different frozen identity')
        if (review['execution_scope'] != 'shared_host_matched_resources'
                or review['uncontended_timing'] is not False or review['scientific_timings_admitted'] is not False
                or review['resource_scopes'] != SCOPES or review['primary_resources_replayed'] is not True):
            raise ValueError('Changed resource scopes, isolation or admission claim')
        eligible = review['shared_host_resources_reviewed']
        if type(eligible) is not bool or review['original_environment_protocol_passed'] is not eligible:
            raise ValueError('Inconsistent eligibility decision')
        status = 'shared_attempt_independently_reviewed' if eligible else 'shared_attempt_reviewed_with_environment_failure'
        if review['status'] != status:
            raise ValueError('Review status and resource eligibility disagree')
        if eligible and (review['scheduler_state'], review['scheduler_exit_code']) != ('COMPLETED', '0:0'):
            raise ValueError('Eligible resource review lacks successful terminal scheduler evidence')
        if not eligible and (review['scheduler_state'], review['scheduler_exit_code']) != ('FAILED', '1:0'):
            raise ValueError('Excluded cadence outcome differs from preserved failure')
        resources = review['resources']
        if set(resources) != set(METRICS):
            raise ValueError('Missing or added resource endpoint')
        for key, value in resources.items():
            if type(value) not in (int, float) or not math.isfinite(value) or value <= 0:
                raise ValueError('Invalid resource measurement')
        if type(resources['peak_memory_bytes']) is not int or resources['peak_memory_bytes'] > 128 * 1024**3:
            raise ValueError('Native peak has wrong units or exceeds the enforced ceiling')
        if resources['wall_seconds'] > 85830:
            raise ValueError('Native wall exceeds the frozen timeout/cleanup boundary')
        by_index[index] = review
        jobs.add(job)
    if set(by_index) != set(range(len(by_index))):
        raise ValueError('Reviewed prefix has a gap; do not skip unresolved attempts')
    rows = []
    for identity in identities:
        review = by_index.get(identity['index'])
        row = dict(identity, status='not_yet_reviewed', job_id=None,
            comparative_timing_eligible=None, resources=None,
            preflight_foreign_average_cores=None, whole_run_maximum_foreign_average_cores=None)
        if review is not None:
            row.update(status='reviewed_shared_observation' if review['shared_host_resources_reviewed'] else 'reviewed_excluded_attempt',
                job_id=review['job_id'], comparative_timing_eligible=review['shared_host_resources_reviewed'],
                resources=review['resources'], preflight_foreign_average_cores=review['preflight_foreign_average_cores'],
                whole_run_maximum_foreign_average_cores=review['whole_run_maximum_foreign_average_cores'])
        rows.append(row)
    cells = []
    for method in METHODS:
        for size in SIZES:
            selected = [row for row in rows if row['method'] == method and row['proteomes'] == size]
            eligible = [row for row in selected if row['comparative_timing_eligible'] is True]
            excluded = [row['index'] for row in selected if row['comparative_timing_eligible'] is False]
            reviewed = len(eligible) + len(excluded)
            complete = len(eligible) == 3
            cells.append(dict(method=method, proteomes=size, planned_repeats=3, reviewed_repeats=reviewed,
                eligible_repeats=len(eligible), excluded_indices=excluded,
                pending_indices=[row['index'] for row in selected if row['status'] == 'not_yet_reviewed'],
                summary_status='three_eligible_repeats' if complete else 'incomplete_eligible_repeats',
                resources={key: dict(median=median([row['resources'][key] for row in eligible]),
                    minimum=min(row['resources'][key] for row in eligible), maximum=max(row['resources'][key] for row in eligible))
                    if complete else None for key in METRICS}))
    return dict(status='shared_host_resource_panel_snapshot', execution_scope='shared_host_matched_resources',
        planned_attempts=27, reviewed_attempts=len(reviews), eligible_attempts=sum(row['shared_host_resources_reviewed'] for row in reviews),
        excluded_attempts=[row['index'] for row in reviews if not row['shared_host_resources_reviewed']],
        all_planned_attempts_reviewed=len(reviews) == 27,
        all_cells_have_three_eligible_repeats=all(cell['eligible_repeats'] == 3 for cell in cells),
        primary_scopes=SCOPES, runs=rows, cells=cells, uncontended_timing=False,
        scientific_timings_admitted=False, publication_ready=False,
        limitations=['Matched-resource, shared-host observations: contention distortion is unknown and potentially method dependent.',
            'Not-yet-reviewed rows can include live or unrun attempts; this offline report is not a scheduler observation.',
            'Excluded raw measurements are retained in the per-attempt table, not pooled into eligible summaries.',
            'Cell medians/ranges require all three eligible repeats; blanks are missing coverage, never zeros.',
            'Observed ranges are descriptive, not confidence intervals; no fastest-repeat selection, overhead subtraction or speedup ranking.',
            'Native CPU includes wrapper bracket work; step peak includes launcher memory, not pure algorithm RSS.',
            'One nested taxon series; proteome count and taxon composition co-vary.',
            'Historical DGX/ARM/shared-host times are not pooled; completed resources do not prove full publication readiness.'])


def collect():
    plan_path = ROOT / 'benchmark_tools/results/threadripper_private_commands_20260928.json'
    plan_ref = executor.record(plan_path)
    plan = executor.read_frozen(plan_path, executor.PRIVATE_PLAN_SHA)
    work = ROOT / 'benchmarks/work/threadripper_shared_execution_20261003'
    paths = [ROOT / 'benchmark_tools/results/threadripper_shared_attempt_22396.json']
    paths += [work / f'review_run_{index:02d}/summary.json' for index in range(1, 27)]
    present = [path for path in paths if path.exists()]
    refs, reviews = [plan_ref], []
    for path in present:
        pin = executor.record(path)
        review = executor.read(pin)
        reviews.append(review)
        refs.append(pin)
        session_ref = review['session']
        session = executor.read(session_ref)
        controller_ref = session['controller']
        terminal_matches(review, session, executor.read(controller_ref), plan_ref)
        refs.extend([session_ref, controller_ref])
        for category in REVIEW_CATEGORIES:
            category_ref = review['reviews'][category]
            value = executor.read(category_ref)
            executor.expect(value, dict(schema='threadripper_panel_review_v1', index=review['index'],
                job_id=review['job_id'], plan_sha256=plan_ref['sha256'], category=category,
                decision='failed' if category == 'environment' and not review['shared_host_resources_reviewed'] else 'passed',
                execution_scope=executor.SHARED_SCOPE, uncontended_timing=False))
            executor.check(value['source'])
            refs.append(value['source'])
            for supporting in value['evidence']:
                executor.check(supporting)
            refs.append(category_ref)
            if category == 'resources':
                resource_refs = [item for item in value['evidence'] if Path(item['path']).name == 'resources.json']
                if len(resource_refs) != 1:
                    raise ValueError('Review omits the independently derived resource record')
                resource = executor.read(resource_refs[0])
                if resource['primary'] != review['resources'] or resource['primary_scopes'] != review['resource_scopes']:
                    raise ValueError('Published summary differs from derived endpoints')
                refs.append(resource_refs[0])
            if category == 'environment':
                replay_refs = [item for item in value['evidence'] if Path(item['path']).name == 'environment_replay.json']
                if len(replay_refs) != 1:
                    raise ValueError('Review omits the independent environment replay')
                replay = executor.read(replay_refs[0])
                environment_matches(review, replay, executor.read(replay['preflight']))
                refs.extend([replay_refs[0], replay['preflight']])
    result = summarize(plan['runs'], reviews)
    for pin in refs:
        executor.check(pin)
    if present != [path for path in paths if path.exists()]:
        raise ValueError('Review inventory changed during collection; preserve this as a snapshot, not a mixed panel')
    result.update(plan=plan_ref, evidence=refs, source=executor.record(__file__))
    return result


def export(output, report):
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)
    with (output / 'panel.json').open('x') as stream:
        json.dump(report, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')
    with (output / 'attempts.tsv').open('x', newline='') as stream:
        columns = [*IDENTITY, 'job_id', 'status', 'comparative_timing_eligible', *METRICS,
            'preflight_foreign_average_cores', 'whole_run_maximum_foreign_average_cores']
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter='\t')
        writer.writeheader()
        for row in report['runs']:
            writer.writerow({**{key: row[key] for key in columns if key not in METRICS},
                **dict.fromkeys(METRICS, None), **(row['resources'] or {})})
    with (output / 'cells.tsv').open('x', newline='') as stream:
        columns = ['method', 'proteomes', 'planned_repeats', 'reviewed_repeats', 'eligible_repeats', 'summary_status',
            *[metric + '_' + statistic for metric in METRICS for statistic in ('median', 'minimum', 'maximum')]]
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter='\t')
        writer.writeheader()
        for cell in report['cells']:
            row = {key: cell[key] for key in columns[:6]}
            for metric in METRICS:
                for statistic in ('median', 'minimum', 'maximum'):
                    row[metric + '_' + statistic] = None if cell['resources'][metric] is None else cell['resources'][metric][statistic]
            writer.writerow(row)
    return [executor.record(output / name) for name in ('panel.json', 'attempts.tsv', 'cells.tsv')]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = collect()
    outputs = export(args.output.resolve(), report)
    print(json.dumps(dict(status=report['status'], reviewed=report['reviewed_attempts'],
        eligible=report['eligible_attempts'], excluded=report['excluded_attempts'], outputs=outputs), indent=2))
