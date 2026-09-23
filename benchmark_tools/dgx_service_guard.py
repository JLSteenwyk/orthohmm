"""Guard the approved Samwise suppression within one held DGX session.

No automatic context-manager cleanup: a lost observation is not permission to
restore a competing service while the bound Slurm allocation may still run.
This guard does not submit jobs or establish whole-host timing eligibility.
"""

import hashlib
from datetime import datetime
import os
from pathlib import Path
import re
import subprocess
import time

from benchmark_tools.probe_dgx_step_separation import save

UNIT = "samwise-daemon-samwise.service"
SOURCE = Path("/home/jlsteenwyk/.config/systemd/user") / UNIT
ENABLED = SOURCE.parent / "default.target.wants" / UNIT
MASK = Path("/run/user/1000/systemd/user.control") / UNIT
LOWER_MASK = Path("/run/user/1000/systemd/user") / UNIT
SOURCE_SHA = "cd2e76207b56adc59dbdaba796995f59ea7dea1320abc97c8fd2ecf5015cb84a"
TERMINAL = {"COMPLETED", "FAILED", "CANCELLED", "TIMEOUT", "NODE_FAIL",
            "OUT_OF_MEMORY", "BOOT_FAIL", "DEADLINE", "PREEMPTED", "REVOKED"}


class ServiceGuard:
    def __init__(self, evidence):
        self.evidence = Path(evidence)
        self.sequence = 0
        self.prior = None
        self.stopped = False
        self.owned = False
        self.job = None
        self.submission_attempted = False
        self.restored = False

    def command(self, argv):
        self.sequence += 1
        started = time.time_ns()
        receipt = {"command": argv, "started_unix_ns": started}
        try:
            result = subprocess.run(argv, capture_output=True, text=True, timeout=15)
            receipt.update(returncode=result.returncode, stdout=result.stdout, stderr=result.stderr)
        except (OSError, subprocess.TimeoutExpired) as error:
            receipt.update(error_type=type(error).__name__, error=str(error))
            raise
        finally:
            receipt["finished_unix_ns"] = time.time_ns()
            save(self.evidence / f"command_{self.sequence:06d}.json", receipt)
        if result.returncode:
            raise RuntimeError(f"Service guard command failed: {argv[0]}")
        return result.stdout

    def ctl(self, *args):
        return self.command(["systemctl", "--user", *args])

    def persistent_identity(self):
        if (hashlib.sha256(SOURCE.read_bytes()).hexdigest() != SOURCE_SHA
                or not ENABLED.is_symlink() or os.readlink(ENABLED) != str(SOURCE)):
            raise ValueError("Persistent Samwise configuration changed")

    def state(self):
        fields = ("LoadState", "ActiveState", "SubState", "MainPID", "FragmentPath")
        raw = self.ctl("show", UNIT, *[f"--property={key}" for key in fields])
        rows = [line.split("=", 1) for line in raw.splitlines()]
        if (any(len(row) != 2 for row in rows) or len(rows) != len(fields)
                or {row[0] for row in rows} != set(fields)):
            raise ValueError("Incomplete or duplicate service state")
        return dict(rows)

    def begin(self):
        if (os.environ.get("ORTHOHMM_APPROVED_SAMWISE_STOP") != "20260920"
                or os.uname().nodename != "spark-7ff0" or os.getuid() != 1000):
            raise ValueError("Require approved DGX user session")
        self.evidence.mkdir(parents=True, exist_ok=False)
        self.persistent_identity()
        if any(p.exists() or p.is_symlink() for p in (MASK, LOWER_MASK)):
            raise ValueError("Refuse preexisting runtime masks")
        state = self.state()
        if state["LoadState"] != "loaded" or state["ActiveState"] not in {
                "active", "activating", "failed", "inactive"}:
            raise ValueError("Unsupported prior service state")
        self.prior = state["ActiveState"]
        save(self.evidence / "prior.json", {"state": state, "source_sha256": SOURCE_SHA})
        # Mark before stop: even a timed-out observation may have stopped it.
        self.stopped = True
        self.ctl("stop", UNIT)
        MASK.parent.mkdir(parents=True, exist_ok=True)
        MASK.symlink_to("/dev/null")
        self.owned = True
        self.ctl("daemon-reload")
        self.ctl("reset-failed", UNIT)
        self.check()

    def check(self):
        if not self.owned or self.restored:
            raise ValueError("No active owned service mask")
        self.persistent_identity()
        if not MASK.is_symlink() or os.readlink(MASK) != "/dev/null":
            raise ValueError("Owned service mask changed or disappeared")
        if self.state() != {"LoadState": "masked", "ActiveState": "inactive",
                "SubState": "dead", "MainPID": "0", "FragmentPath": str(MASK)}:
            raise ValueError("Service suppression not verified")

    def before_submission(self):
        self.check()
        if self.submission_attempted:
            raise ValueError("No automatic repeat submission")
        self.submission_attempted = True
        save(self.evidence / "submission_attempted.json", {"automatic_retry": False})

    def bind_job(self, job):
        if type(job) is not int or job <= 0 or self.job is not None:
            raise ValueError("Require one positive job identity")
        # Bind before observation: a failed check must not erase a live job.
        self.submission_attempted = True
        self.job = job
        save(self.evidence / "job.json", {"job_id": job})
        self.check()

    def restore(self):
        if not self.stopped or self.restored:
            raise ValueError("No pending restoration")
        if self.submission_attempted and self.job is None:
            raise ValueError("Submission identity unresolved; retain suppression")
        if self.job is not None:
            raw = self.command(["scontrol", "show", "job", str(self.job), "--oneliner"])
            lines = [line for line in raw.splitlines() if line.strip()]
            pairs = re.findall(r"(?<!\S)([A-Za-z][^\s=]*)=([^\s]+)", raw)
            fields = dict(pairs)
            if (len(lines) != 1 or len(fields) != len(pairs)
                    or fields.get("JobId") != str(self.job)
                    or fields.get("JobState") not in TERMINAL
                    or any(key in fields for key in ("ArrayJobId", "ArrayTaskId", "HetJobId"))
                    or fields.get("EndTime") in {None, "", "Unknown", "None"}):
                raise ValueError("Bound job is not authoritatively terminal; retain suppression")
            datetime.strptime(fields["EndTime"], "%Y-%m-%dT%H:%M:%S")
        self.persistent_identity()
        if self.owned:
            if not MASK.is_symlink() or os.readlink(MASK) != "/dev/null":
                raise ValueError("Mask changed or disappeared; refuse automatic restoration")
            MASK.unlink()
            self.owned = False
        self.ctl("daemon-reload")
        if self.prior != "inactive":
            self.ctl("start", UNIT)
        state = self.state()
        self.persistent_identity()
        if state["LoadState"] != "loaded" or MASK.exists() or MASK.is_symlink():
            raise ValueError("Restored service configuration not verified")
        save(self.evidence / "restored.json", {"prior_active_state": self.prior,
            "state": state, "start_requested": self.prior != "inactive",
            "application_health_verified": False, "scientific_timings_admitted": False})
        self.restored = True
