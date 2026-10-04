"""Own one same-allocation preflight child and retain its terminal disposition."""

import json
import os
from pathlib import Path
import subprocess
import sys
import time

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.slurm_resource_snapshot import scoped_path


class EnvironmentWorker:
    def __init__(self, session, request_ref, policy_ref, root, *, popen=subprocess.Popen):
        self.session, self.root = Path(session), Path(root)
        self.request_ref, self.policy_ref, self.popen = request_ref, policy_ref, popen
        self.process = self.log = None
        self.entered = self.joined = False
        self.report = dict(status="not_started", automatic_retry=False, scientific_timings_admitted=False)

    def __enter__(self):
        if self.entered:
            raise RuntimeError("An environmental worker owner cannot be reused")
        self.entered = True
        cache = self.session / "environment_worker_python_cache"
        if cache.exists() or cache.is_symlink():
            raise FileExistsError(cache)
        env = os.environ.copy()
        env.update(PYTHONPATH=str(self.root), PYTHONNOUSERSITE="1", PYTHONHASHSEED="0",
                   PYTHONDONTWRITEBYTECODE="1", PYTHONPYCACHEPREFIX=str(cache))
        for name in ("PYTHONHOME", "LD_PRELOAD", "LD_LIBRARY_PATH", "LD_AUDIT"):
            env.pop(name, None)
        command = [sys.executable, "-B", "-m", "benchmark_tools.threadripper_environment_worker",
                   "--request", self.request_ref["path"], "--request-sha256", self.request_ref["sha256"],
                   "--policy", self.policy_ref["path"], "--policy-sha256", self.policy_ref["sha256"]]
        self.report.update(command=command, cwd=str(self.root), request=self.request_ref,
                           policy=self.policy_ref, started_unix_ns=time.time_ns(),
                           slurm_job_id=env.get("SLURM_JOB_ID"), cache_path=str(cache))
        try:
            self.log = (self.session / "environment_worker.log").open("x")
            self.process = self.popen(command, cwd=self.root, env=env, stdin=subprocess.DEVNULL,
                                      stdout=self.log, stderr=subprocess.STDOUT)
            self.report.update(status="started", pid=self.process.pid)
            save(self.session / "environment_worker_started.json", self.report)
        except BaseException as error:
            self.__exit__(type(error), error, error.__traceback__)
            raise
        return self

    def wait_response(self, path, seconds=20):
        deadline = time.monotonic() + seconds
        while not Path(path).exists():
            if self.process.poll() is not None and not Path(path).exists():
                raise RuntimeError("Environmental worker exited without a response")
            if time.monotonic() >= deadline:
                raise TimeoutError("Environmental worker response deadline exceeded")
            time.sleep(.02)
        return json.loads(Path(path).read_text())

    def wait_prepared(self, seconds=4000):
        if self.process is None or self.report.get("prepared") is not None:
            raise RuntimeError("Preparation requires one newly owned environmental worker")
        path = self.session / "environment_worker_prepared.json"
        prepared = self.wait_response(path, seconds)
        ref = record(path)
        check(self.request_ref)
        check(self.policy_ref)
        request = json.loads(Path(self.request_ref["path"]).read_text())
        expected = dict(schema="threadripper_environment_worker_prepared_v1",
            status="prepared_waiting_for_release_request", request=self.request_ref,
            policy=self.policy_ref, job_id=request["job_id"], index=request["index"],
            pid=self.process.pid, native_release_authorized=False, scientific_timings_admitted=False)
        if any(type(prepared.get(key)) is not type(value) or prepared[key] != value
               for key, value in expected.items()):
            raise ValueError("Environmental preparation belongs to another worker or request")
        stamp = prepared.get("prepared_unix_ns")
        if (type(stamp) is not int or not self.report["started_unix_ns"] <= stamp <= time.time_ns()
                or prepared.get("boot_id") != Path("/proc/sys/kernel/random/boot_id").read_text().strip()):
            raise ValueError("Environmental preparation has stale or invalid identity")
        scope = scoped_path(Path(f"/proc/{self.process.pid}/cgroup").read_text(), request["job_id"])
        job_scope = str(Path(*scope.parts[:scope.parts.index(f"job_{request['job_id']}") + 1]))
        if prepared.get("job_scope") != job_scope or self.process.poll() is not None:
            raise ValueError("Prepared environmental worker is not live in this allocation")
        if json.loads(path.read_text()) != prepared:
            raise ValueError("Environmental preparation changed during validation")
        for pin in (self.request_ref, self.policy_ref, ref):
            check(pin)
        self.report["prepared"] = ref
        return ref

    def finish(self):
        code = self.process.wait(timeout=5)
        self.report["exit_code"] = code
        if code != 0:
            raise RuntimeError(f"Environmental worker exited with status {code}")
        self.joined = True

    def __exit__(self, kind, error, traceback):
        cleanup_error = None
        try:
            if self.process is not None:
                if self.process.poll() is None:
                    self.report["cleanup_requested"] = True
                    self.process.terminate()
                    try:
                        self.process.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        self.report["kill_requested"] = True
                        self.process.kill()
                        self.process.wait(timeout=10)
                self.report["exit_code"] = self.process.returncode
                self.report["terminal"] = self.process.poll() is not None
            else:
                self.report["terminal"] = True
        except BaseException as failure:
            cleanup_error = failure
            self.report.update(cleanup_error_type=type(failure).__name__, cleanup_error=str(failure),
                               terminal=self.process is not None and self.process.poll() is not None)
        finally:
            if self.log is not None:
                self.log.close()
                self.report["log"] = record(self.session / "environment_worker.log")
            self.report.update(status="completed" if kind is None and cleanup_error is None
                               and self.joined
                               and not self.report.get("cleanup_requested")
                               and self.report.get("exit_code") == 0 else "failed_or_cancelled",
                               finished_unix_ns=time.time_ns(),
                               parent_error_type=None if kind is None else kind.__name__)
            save(self.session / "environment_worker_lifecycle.json", self.report)
        if cleanup_error is not None and kind is None:
            raise cleanup_error
        if kind is None and not self.joined:
            raise RuntimeError("Environmental worker was not successfully joined at release")
        return False
