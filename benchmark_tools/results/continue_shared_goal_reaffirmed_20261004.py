"""Rebind only the authorized goal reaffirmation, preserving the prepared panel."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from benchmark_tools.results import continue_shared_prepared_panel_20261004 as previous

executor, WORK, RESULTS = previous.executor, previous.WORK, previous.RESULTS
GOAL = ROOT / 'benchmark_tools/PUBLICATION_GOAL_20261003.txt'
FILES = {key: RESULTS / f'threadripper_goal_reaffirmed_{key}_20261004.json'
    for key in ('policy', 'readiness', 'preparation')}
ADDITION = ('  > Shared-host timing authorization reaffirmed 2026-10-04: on resumption, '
    'proceed on the current Threadripper despite competing analyses, provided safe capacity '
    'and valid accounting are maintained. Do not wait for a quiet window or require the DGX. '
    'Record contention for each run and disclose in timing tables, figures, Methods and '
    'limitations that its effect is unknown and may differ between tools; matching resource '
    'limits does not establish isolated performance. Preserve all other scientific requirements, '
    'retained results and failure-handling rules. This update authorizes continuation under '
    'that scope; it does not request restarting completed runs or disrupting unrelated work.\n  >\n')


def reaffirmed(old, new):
    first, rest = old.split(b'  >\n', 1)
    if new != first + b'  >\n' + ADDITION.encode('ascii') + rest:
        raise ValueError('Only the authorized goal reaffirmation may change')


def replace_pin(refs, old, new):
    if sum(pin == old for pin in refs) != 1:
        raise ValueError('Require exactly one original documentary binding')
    return [new if pin == old else pin for pin in refs]


def prepare():
    if any(path.exists() or path.is_symlink() for path in FILES.values()):
        raise FileExistsError('Goal reaffirmation preparation already exists')
    previous.load()
    parent_ref = executor.record(previous.PREPARATION)
    parent = executor.read(parent_ref)
    for key in ('lookup', 'binding', 'plan', 'recipe', 'policy', 'readiness', 'resources'):
        executor.check(parent[key])
    policy = executor.read(parent['policy'])
    ready = executor.read(parent['readiness'])
    goals = [pin for pin in policy['evidence'] if pin['path'] == str(GOAL)]
    if len(goals) != 1:
        raise ValueError('Missing original goal binding')
    old_goal = goals[0]
    old_bytes = subprocess.check_output(['git', 'show',
        parent['source_commit'] + ':' + str(GOAL.relative_to(ROOT))], cwd=ROOT)
    if (len(old_bytes), hashlib.sha256(old_bytes).hexdigest()) != (old_goal['bytes'], old_goal['sha256']):
        raise ValueError('Retained goal does not match preparation commit')
    goal_ref = executor.record(GOAL)
    reaffirmed(old_bytes, GOAL.read_bytes())
    source_ref = executor.record(__file__)
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    for pin in (source_ref, goal_ref):
        path = Path(pin['path'])
        if path.read_bytes() != subprocess.check_output(['git', 'show',
                commit + ':' + str(path.relative_to(ROOT))], cwd=ROOT):
            raise ValueError('Reaffirmation source/goal must be committed before preparation')
    policy_new = dict(policy, evidence=replace_pin(policy['evidence'], old_goal, goal_ref))
    for pin in policy_new['evidence']:
        executor.check(pin)
    policy_ref = previous.save(FILES['policy'], policy_new)
    ready_new = dict(ready, environment_policy=policy_ref,
        evidence=replace_pin(replace_pin(ready['evidence'], old_goal, goal_ref), parent['policy'], policy_ref))
    for pin in ready_new['evidence']:
        executor.check(pin)
    ready_ref = previous.save(FILES['readiness'], ready_new)
    result = dict(parent, policy=policy_ref, readiness=ready_ref,
        documentary_reaffirmation=dict(parent=parent_ref, old_goal=old_goal, goal=goal_ref,
            source=source_ref, source_commit=commit, historical_goal_verified_from_git=True,
            only_documentary_evidence_changed=True, native_runs_started=False))
    return previous.save(FILES['preparation'], result)


def launch(index):
    _, _, plan, plan_ref = previous.load()
    prefix, _ = previous.history(index, plan_ref)
    if not 24 <= index < 27:
        raise ValueError('Reaffirmed continuation covers only remaining indices 24-26')
    preparation = executor.read(executor.record(FILES['preparation']))
    for key in ('lookup', 'binding', 'plan', 'recipe', 'policy', 'readiness', 'resources'):
        executor.check(preparation[key])
    executor.review(preparation['policy'], dict(decision='reviewed'))
    executor.review(preparation['readiness'], dict(decision='passed'))
    run = plan['runs'][index]
    request_path = WORK / f'request_run_{index:02d}.json'
    session = Path(run['measurement_directory']).parent.parent / 'sessions' / f'run_{index:02d}'
    if request_path.exists() or session.exists() or Path(run['measurement_directory']).parent.exists():
        raise FileExistsError('Existing attempt must be observed/reviewed, never retried')
    suffix = f'{index:02d}'
    raw = previous.first.command('submission_' + suffix, ['sbatch', '--hold', '--parsable', '--partition=gpu',
        '--chdir=' + str(ROOT), '--output=' + str(WORK / f'slurm_run_{suffix}-%j.log'),
        '--error=' + str(WORK / f'slurm_run_{suffix}-%j.log'), str(previous.SCRIPT),
        '--request', str(request_path), '--request-sha256', 'scheduler-comment'])
    match = re.fullmatch(r'([0-9]+)(?:;[^\s]+)?\n?', raw)
    if match is None:
        raise ValueError('Cannot identify held job; preserve submission and inspect scheduler')
    job = int(match[1])
    held = previous.first.command('held_controller_' + suffix, ['scontrol', 'show', 'job', str(job), '--oneliner'])
    pairs = re.findall(r'(?<!\S)([A-Za-z][^\s=]*)=([^\s]+)', held)
    fields = dict(pairs)
    required = dict(JobId=str(job), JobState='PENDING', Reason='JobHeldUser', NumCPUs='64',
        NumTasks='1', MinMemoryNode='128G', TimeLimit=executor.TIME_LIMIT, Requeue='0',
        Command=str(previous.SCRIPT), WorkDir=str(ROOT))
    required['CPUs/Task'] = '64'
    if len(fields) != len(pairs) or any(fields.get(key) != value for key, value in required.items()):
        raise ValueError('Held resources/identity differ; do not release')
    request = dict(schema='threadripper_execution_request_v1', job_id=job,
        execution_authorized=True, deployment='private_v2_20260928', execution_scope=executor.SHARED_SCOPE,
        runtime_lookup=preparation['lookup'], plan_sha256=plan_ref['sha256'], lookup_sha256=previous.LOOKUP_SHA,
        index=index, allocation_cwd=str(ROOT), scheduler_command=str(previous.SCRIPT),
        recipe=preparation['recipe'], history=prefix, readiness_review=preparation['readiness'],
        environment_preflight_path=str(session / 'environment_preflight.json'))
    request_ref = previous.save(request_path, request)
    chosen, _, selected_history, sources = executor.select(request, ROOT, job)
    if chosen != run or selected_history['progress']['index'] != index:
        raise ValueError('Actual held-job selection differs from frozen next identity')
    for pin in [request_ref, *sources]:
        executor.check(pin)
    previous.first.command('request_comment_' + suffix, ['scontrol', 'update', 'JobId=' + str(job),
        'Comment=' + request_ref['sha256']])
    pinned = previous.first.command('pinned_controller_' + suffix, ['scontrol', 'show', 'job', str(job), '--oneliner'])
    if ('Comment=' + request_ref['sha256']) not in pinned.split():
        raise ValueError('Request digest not bound in scheduler; leave held')
    previous.first.command('release_' + suffix, ['scontrol', 'release', str(job)])
    result = dict(status='repaired_shared_job_bound_and_released', job_id=job, index=index,
        method=run['method'], proteomes=run['proteomes'], repeat=run['repeat'], source=executor.record(__file__),
        request=request_ref, preparation=executor.record(FILES['preparation']), execution_scope=executor.SHARED_SCOPE,
        automatic_retry=False, uncontended_timing=False, scientific_timings_admitted=False,
        native_phase='requires_live_observation')
    previous.save(WORK / f'launch_{suffix}.json', result)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['prepare', 'launch'])
    parser.add_argument('--index', type=int)
    args = parser.parse_args()
    if args.action == 'launch' and args.index is None:
        parser.error('launch requires --index')
    print(json.dumps(prepare() if args.action == 'prepare' else launch(args.index), indent=2))
