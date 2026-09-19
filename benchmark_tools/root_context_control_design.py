"""Fixed root-context engineering conditions and finite owned-service command."""

from pathlib import PurePosixPath

ORDER = (("idle", "steady", "churn", "user-contended"),
         ("steady", "user-contended", "idle", "churn"),
         ("churn", "idle", "user-contended", "steady"))
MODES = {mode for block in ORDER for mode in block}


def unit_name(job, index):
    if type(job) is not int or job <= 0 or type(index) is not int or not 0 <= index < 12:
        raise ValueError("Invalid owned-control identity")
    return f"orthohmm-root-{job}-{index}.service"


def check_service_scope(membership, manager, unit):
    if not membership.startswith("0::/") or len(membership.strip().splitlines()) != 1:
        raise ValueError("Require unified service membership")
    scope = PurePosixPath(membership.strip()[3:])
    base = PurePosixPath(manager)
    if (not manager.startswith("/user.slice/") or ".." in base.parts or str(base) != manager
            or base not in scope.parents or scope.name != unit or ".." in scope.parts):
        raise ValueError("Service outside expected user manager or wrong owned unit")
    return str(scope)


def service_command(python, source, directory, unit, cpu):
    if type(cpu) is not int or cpu < 0:
        raise ValueError("Invalid service CPU")
    return ["systemd-run", "--user", "--quiet", "--wait", "--pipe", "--collect", "--unit=" + unit,
        "--working-directory=" + str(source.parent.parent), "--property=RuntimeMaxSec=45s",
        "--property=MemoryMax=256M", "--property=TasksMax=8", "taskset", "-c", str(cpu),
        python, "-B", str(source), "--directory", str(directory), "--competitor", str(cpu)]
