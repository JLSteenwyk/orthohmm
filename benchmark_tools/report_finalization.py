"""Retain reporting-stage duration and a later job-cgroup memory observation."""

from contextlib import contextmanager
import time

from benchmark_tools.probe_dgx_step_separation import save


@contextmanager
def observe(directory, job_id, scope, read_memory, *, clock=time.monotonic_ns):
    receipt = dict(schema="threadripper_report_finalization_v1", job_id=job_id,
        scope=str(scope), started_ns=clock(), status="reporting_failed",
        scientific_timings_admitted=False,
        limitations=["Reporting duration excludes this receipt and its final memory read.",
            "Job peak includes preparation, inference and observer/reporting through the read; not complete-job peak.",
            "Job and native-step peaks overlap; no subtraction or addition is valid.",
            "Abrupt termination can prevent this receipt; absence is incomplete evidence."])
    failure = None
    try:
        yield
        receipt["status"] = "reporting_completed"
    except BaseException as error:
        failure = error
        receipt.update(error_type=type(error).__name__, error=str(error))
        raise
    finally:
        receipt["finished_ns"] = clock()
        receipt["reporting_wall_s"] = (receipt["finished_ns"] - receipt["started_ns"]) / 1e9
        try:
            receipt["job_memory"] = read_memory(scope)
        except Exception as error:
            receipt.update(status="reporting_memory_observation_failed",
                memory_error_type=type(error).__name__, memory_error=str(error))
            save(directory / "report_finalization.json", receipt)
            if failure is None:
                raise
        else:
            save(directory / "report_finalization.json", receipt)


def validate(receipt, job_id, after_native):
    """Validate phase identity/time; the caller separately validates raw memory."""
    if (receipt.get("schema") != "threadripper_report_finalization_v1"
            or receipt.get("status") != "reporting_completed"
            or type(receipt.get("job_id")) is not int or receipt["job_id"] != job_id
            or receipt.get("scientific_timings_admitted") is not False
            or receipt.get("scope") != after_native["scope"]):
        raise ValueError("Incomplete or mismatched report finalization")
    start, end = receipt["started_ns"], receipt["finished_ns"]
    memory = receipt["job_memory"]
    if (type(start) is not int or type(end) is not int
            or not after_native["finished_ns"] < start <= end < memory["started_ns"]
            or type(receipt["reporting_wall_s"]) not in (int, float)
            or receipt["reporting_wall_s"] != (end - start) / 1e9):
        raise ValueError("Invalid report finalization timing")
    return memory
