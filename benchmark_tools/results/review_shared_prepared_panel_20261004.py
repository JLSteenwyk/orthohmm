"""Review a terminal shared-host attempt without submitting or retrying work."""
import argparse
import csv
import io
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from benchmark_tools.results import continue_shared_prepared_panel_20261004 as launch
from benchmark_tools import run_threadripper_scaling as executor
from benchmark_tools.audit_threadripper_native_outcome import audit
from benchmark_tools.derive_threadripper_resources import endpoints, protocol
from benchmark_tools.review_threadripper_process_stream import evaluate as process_evaluate
from benchmark_tools.review_threadripper_pressure_stream import evaluate as pressure_evaluate
from benchmark_tools.slurm_resource_snapshot import scoped_path
from benchmark_tools.verify_threadripper_controller import validate


def require(value, expected):
    executor.expect(value, expected)


def review(index, revision=1, prior_native=None, terminal_controller=None):
    lookup, binding, plan, plan_ref = launch.load()
    if type(index) is not int or not 21 <= index < len(plan['runs']):
        raise ValueError('Invalid frozen index')
    run = plan['runs'][index]
    request_ref = executor.record(launch.WORK / f'request_run_{index:02d}.json')
    request = executor.read(request_ref)
    job = request['job_id']
    argv = ['scontrol', 'show', 'job', str(job), '--oneliner']
    accounting = None
    if terminal_controller is None:
        started = time.time_ns()
        query = subprocess.run(argv, text=True, capture_output=True, timeout=5, check=True)
        observed = dict(command=argv, returncode=query.returncode, stdout=query.stdout,
            stderr=query.stderr, started_unix_ns=started, finished_unix_ns=time.time_ns())
    else:
        observed = executor.read(executor.record(terminal_controller))
        if observed['command'] != argv or observed['returncode'] != 0:
            raise ValueError('Retained terminal controller queried another job or failed')
        accounting_argv = ['sacct', '-j', str(job), '--noheader', '--parsable2',
            '--format=JobID,State,ExitCode,Start,End']
        started = time.time_ns()
        query = subprocess.run(accounting_argv, text=True, capture_output=True, timeout=5, check=True)
        accounting = dict(command=accounting_argv, returncode=query.returncode,
            stdout=query.stdout, stderr=query.stderr, started_unix_ns=started,
            finished_unix_ns=time.time_ns(), retained_controller=executor.record(terminal_controller))
    allocation = validate(observed['stdout'], job, 'terminal', command=str(launch.SCRIPT),
        cwd=str(ROOT), time_limit=executor.TIME_LIMIT, allocation_mode='shared')
    if accounting is not None:
        rows = [row for row in csv.reader(io.StringIO(query.stdout), delimiter='|') if row and row[0] == str(job)]
        fields = allocation['fields']
        if rows != [[str(job), fields['JobState'], fields['ExitCode'], fields['StartTime'], fields['EndTime']]]:
            raise ValueError('Fresh accounting does not corroborate the retained terminal controller')
    if type(revision) is not int or revision < 1:
        raise ValueError('Invalid review revision')
    directory = launch.WORK / (f'review_run_{index:02d}' if revision == 1 else f'review_run_{index:02d}_revision_{revision}')
    directory.mkdir(exist_ok=False)
    controller_ref = launch.save(directory / 'controller.json', observed)
    accounting_ref = launch.save(directory / 'accounting.json', accounting) if accounting is not None else None
    try:
        if allocation['fields'].get('Comment') != request_ref['sha256']:
            raise ValueError('Terminal scheduler request digest differs')
        scheduler_succeeded = allocation['scheduler_state'] == 'COMPLETED' and allocation['scheduler_exit_code'] == '0:0'
        selected, selected_lookup, _, sources = executor.select(request, ROOT, job)
        if selected != run or selected_lookup != lookup:
            raise ValueError('Terminal attempt differs from frozen selection')
        measurement = Path(run['measurement_directory'])
        session = measurement.parent.parent / 'sessions' / f'run_{index:02d}'
        verification_ref = executor.record(measurement.parent / 'verification.json')
        verification = executor.read(verification_ref)
        result_ref = executor.record(session / 'result.json')
        result = executor.read(result_ref)
        require(result, dict(index=index, job_id=job, execution_scope=executor.SHARED_SCOPE,
            uncontended_timing=False, next_submission_authorized=False))
        if not scheduler_succeeded and (allocation['scheduler_state'] != 'FAILED'
                or result.get('status') != 'executor_failed'
                or result.get('error') != 'Whole-run sampled environment policy was not satisfied'):
            raise ValueError('Unclassified infrastructure failure requires separate review')
        if scheduler_succeeded and result.get('status') != 'measurement_returned_pending_independent_review':
            raise ValueError('Successful scheduler and executor disagree')
        if result['request'] != request_ref or result['wrapper'] != verification:
            raise ValueError('Executor and independently retained wrapper differ')
        require(verification, dict(status='command_exited_zero', scientific_results_admitted=False))
        if verification['source_sha256'] != executor.record(ROOT / 'benchmark_tools/run_verified_slurm_measurement.py')['sha256']:
            raise ValueError('Runtime verification wrapper source differs')
        runtime_refs = []
        for number, side in enumerate(('before', 'after'), start=1):
            check_ref = executor.record(session / f'lookup_checks/checked_{number:02d}.json')
            checked = executor.read(check_ref)
            require(checked, dict(status='runtime_and_lookup_checked', scientific_execution_authorized=False))
            if verification[side]['runtime'] != checked:
                raise ValueError('Wrapper and lookup check disagree')
            expected_runtime = []
            for path, sha in binding['runtime_specs']:
                manifest = executor.read_frozen(Path(path), sha)
                expected_runtime.append(dict(path=path, sha256=sha,
                    status='runtime_tree_identity_matches', records=len(manifest['records']),
                    scientific_execution_authorized=False))
            if checked['runtime'] != expected_runtime:
                raise ValueError('Runtime trees or retained record counts differ')
            for name in ('orthohmm', 'orthofinder'):
                require(checked['lookup'][name], dict(status='native_python_lookup_matches'))
                executor.check(checked['lookup'][name]['report'])
                runtime_refs.append(checked['lookup'][name]['report'])
            if verification[side]['original_inputs'] != run['dataset']['inputs']:
                raise ValueError('Runtime input verification differs')
            runtime_refs.append(check_ref)
        baseline = executor.read(lookup['baseline'])
        if prior_native is None:
            native = audit(run, baseline, job)
            native_ref = launch.save(directory / 'native_audit.json', native)
        else:
            native_ref = executor.record(prior_native)
            native = executor.read(native_ref)
            if native['source'] != executor.record(ROOT / 'benchmark_tools/audit_threadripper_native_outcome.py'):
                raise ValueError('Retained native audit source differs')
            if native['replay']['measured'] != verification['measurement']:
                raise ValueError('Retained native audit belongs to another measurement')
            for ref in [*native['evidence'], *native['replay']['evidence'], *native['outputs']['native']['checked_files']]:
                executor.check(ref)
        require(native, dict(status='native_success_outputs_verified', native_outcome='exited_zero', job_id=job))
        prepared = executor.read(executor.record(launch.PREPARATION))
        resource_protocol = prepared['resources']
        ready_policy = executor.read(executor.read(request['readiness_review'])['environment_policy'])
        if resource_protocol not in ready_policy['evidence']:
            raise ValueError('Resource protocol is not bound by the actual launch policy')
        protocol_ref, _ = protocol(Path(resource_protocol['path']), resource_protocol['sha256'])
        replayed = native['replay']
        resource = endpoints(replayed, replayed['measured']['native'], job)
        resource.update(execution_scope=executor.SHARED_SCOPE, protocol=protocol_ref,
            uncontended_timing=False, contention_distortion='unknown_potentially_method_dependent')
        resources_ref = launch.save(directory / 'resources.json', resource)
        preflight_ref = executor.record(session / 'environment_preflight.json')
        preflight = executor.read(preflight_ref)
        require(preflight, dict(job_id=job, index=index, execution_scope=executor.SHARED_SCOPE,
            decision='passed', background_competition_recorded=True, foreign_cpu_used_for_eligibility=False,
            preflight_pressure_limits_used=False, uncontended_timing=False))
        readiness = executor.read(request['readiness_review'])
        policy_ref = readiness['environment_policy']
        policy = executor.read(policy_ref)
        if not executor.shared_environment(policy) or preflight['environment_policy'] != policy_ref:
            raise ValueError('Shared preflight policy differs')
        if any(value < policy['minimum_available_memory_bytes'] for value in preflight['available_memory_bytes']):
            raise ValueError('Launch capacity did not meet the prospective minimum')
        retained_ref = result['process_stream_review']
        retained = executor.read(retained_ref)
        require(retained, dict(job_id=job, index=index, execution_scope=executor.SHARED_SCOPE,
            background_cpu_used_for_eligibility=False,
            pressure_thresholds_used_for_eligibility=False, uncontended_timing=False))
        for ref in [*preflight['evidence'], *retained['evidence']]:
            executor.check(ref)
        handoff = launch.observer_handoff(session, measurement, request_ref, policy_ref, index, job, terminal=True)
        handoff_ref = launch.save(directory / 'prepared_handoff.json', handoff)
        ready = executor.read(executor.record(measurement / 'ready.json'))
        scope = scoped_path(ready['cgroup'], job)
        job_scope = next(path for path in scope.parents if path.name == f'job_{job}')
        done = replayed['measured']['native']
        with (measurement / 'host_processes.jsonl').open() as stream:
            observer = json.loads(stream.readline())['observer_pid']
            stream.seek(0)
            processes = process_evaluate(stream, executor.read(policy['process_policy']),
                boot_id=preflight['boot_id'], job_scope=str(job_scope), observer_pid=observer,
                launch=done['started_ns']/1e9, end=done['finished_ns']/1e9,
                maximum_foreign_average_cores=policy['maximum_foreign_average_cores'],
                maximum_sample_period_s=policy['maximum_sample_period_s'])
        for key, value in processes.items():
            if key not in {'schema', 'limitations'} and retained.get(key) != value:
                raise ValueError('Independent process replay differs: ' + key)
        points = sorted(measurement.glob('point_*.json'))
        if points != [measurement / f'point_{i:06d}.json' for i in range(len(points))]:
            raise ValueError('Pressure point inventory differs')
        pressure = pressure_evaluate((json.loads(path.read_text()) for path in points),
            boot_id=preflight['boot_id'], job_scope=str(job_scope), launch_ns=done['started_ns'],
            end_ns=done['finished_ns'], limits=policy['maximum_pressure_percent'],
            maximum_period_s=policy['maximum_pressure_sample_period_s'], pressure_role='diagnostic_only')
        environment_passed = processes['sampled_process_policy_satisfied'] and pressure['sampled_pressure_evidence_satisfied']
        if retained['pressure_review'] != pressure or retained['sampled_environment_policy_satisfied'] != environment_passed:
            raise ValueError('Independent shared environment replay differs')
        if scheduler_succeeded != environment_passed:
            raise ValueError('Scheduler outcome does not match the reviewed environment gate')
        env_ref = launch.save(directory / 'environment_replay.json',
            dict(processes=processes, pressure=pressure, preflight=preflight_ref,
                retained_review=retained_ref, execution_scope=executor.SHARED_SCOPE,
                uncontended_timing=False, scientific_timings_admitted=False))
        evidence = dict(runtime=[request_ref, verification_ref, result_ref, *runtime_refs,
            request['recipe'], request['runtime_lookup'], lookup['binding'], lookup['baseline']],
            resources=[native_ref, resources_ref, protocol_ref],
            environment=[preflight_ref, retained_ref, env_ref, policy_ref, handoff_ref, *handoff['evidence']],
            outputs_or_failure=[native_ref])
        for ref in [request_ref, controller_ref, *sources, *native['evidence'],
                    *retained['evidence'], *preflight['evidence'], *runtime_refs]:
            executor.check(ref)
        reviews = {category: launch.save(directory / (category + '_review.json'),
            dict(schema='threadripper_panel_review_v1', index=index, job_id=job,
                plan_sha256=plan_ref['sha256'], category=category,
                decision='failed' if category == 'environment' and not environment_passed else 'passed',
                evidence=refs, source=executor.record(__file__),
                execution_scope=executor.SHARED_SCOPE, uncontended_timing=False,
                review_reference='Independent terminal identity, raw replay and native-output review under the prospective shared-host policy'))
            for category, refs in evidence.items()}
        session_ref = launch.save(directory / 'session.json',
            dict(schema='threadripper_panel_session_v1', index=index, job_id=job,
                plan_sha256=plan_ref['sha256'], phase='terminal', controller=controller_ref,
                native_outcome='exited_zero', native_audit=native_ref, reviews=reviews))
        summary = dict(status='shared_attempt_independently_reviewed' if environment_passed else 'shared_attempt_reviewed_with_environment_failure',
            index=index, job_id=job, scheduler_state=allocation['scheduler_state'],
            scheduler_exit_code=allocation['scheduler_exit_code'],
            method=run['method'], proteomes=run['proteomes'], repeat=run['repeat'],
            session=session_ref, reviews=reviews, resources=resource['primary'],
            terminal_accounting=accounting_ref, prepared_handoff=handoff_ref,
            resource_scopes=resource['primary_scopes'], execution_scope=executor.SHARED_SCOPE,
            native_outputs=native['outputs']['native'],
            preflight_foreign_average_cores=preflight['observed_foreign_average_cores'],
            whole_run_maximum_foreign_average_cores=processes['maximum_observed_foreign_average_cores'],
            uncontended_timing=False, scientific_timings_admitted=False,
            shared_host_resources_reviewed=environment_passed, primary_resources_replayed=True,
            original_environment_protocol_passed=environment_passed,
            original_environment_failures=dict(processes=processes['failures'], pressure=pressure['failures']),
            next_submission_authorized=False,
            limitations=['Shared-host observations only; contention distortion is unknown and potentially method dependent.',
                'Periodic observations are not continuous containment or isolation certification.',
                'CPU includes wrapper bracket work; memory is native-step lifetime peak, not pure algorithm RSS.',
                'All-prefix history and a fresh launch preflight remain required for the next frozen identity.'])
        summary_ref = launch.save(directory / 'summary.json', summary)
        return {**summary, 'summary': summary_ref}
    except BaseException as error:
        launch.save(directory / 'failure.json', dict(status='shared_attempt_review_failed',
            index=index, job_id=job, error_type=type(error).__name__, error=str(error),
            controller=controller_ref, automatic_retry=False, next_submission_authorized=False))
        raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--index', type=int, required=True)
    parser.add_argument('--revision', type=int, default=1)
    parser.add_argument('--prior-native-audit', type=Path)
    parser.add_argument('--terminal-controller', type=Path)
    args = parser.parse_args()
    result = review(args.index, args.revision, args.prior_native_audit, args.terminal_controller)
    print(json.dumps({key: result[key] for key in ('status', 'job_id', 'index', 'resources', 'summary')}, indent=2))
