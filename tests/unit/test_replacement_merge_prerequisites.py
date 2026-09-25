import pytest

from benchmark_tools.run_blast_recovery_merge import require_completed


def accounting(replacement):
    rows = ["JobID|State|ExitCode|NodeList|AllocCPUS|Elapsed", "22148|COMPLETED|0:0|bizon|2|00:01:00"]
    for i in range(20):
        rows.append(f"22103_{i}|{'TIMEOUT' if replacement and i == 14 else 'COMPLETED'}|0:0|bizon|180|00:01:00")
        rows.append(f"22105_{i}|{'PENDING' if replacement and i == 14 else 'COMPLETED'}|0:0|bizon|2|00:01:00")
    if replacement:
        rows.extend(["22160_14|COMPLETED|0:0|bizon|180|00:01:00", "22161|COMPLETED|0:0|bizon|2|00:01:00"])
    return "\n".join(rows) + "\n"


def test_original_panel_unchanged():
    require_completed(accounting(False))


def test_explicit_replacement_required():
    require_completed(accounting(True), True)
    with pytest.raises(ValueError):
        require_completed(accounting(True))
    with pytest.raises(ValueError):
        require_completed(accounting(False), True)


@pytest.mark.parametrize("job", ["22160_14", "22161", "22148", "22103_13", "22105_19"])
@pytest.mark.parametrize("problem", ["missing", "duplicate", "running", "exit", "cpu"])
def test_each_required_job_must_pass(job, problem):
    rows = accounting(True).splitlines()
    selected = next(row for row in rows if row.startswith(job + "|"))
    if problem == "missing":
        rows.remove(selected)
    elif problem == "duplicate":
        rows.append(selected)
    else:
        fields = selected.split("|")
        index, value = {"running": (1, "RUNNING"), "exit": (2, "1:0"), "cpu": (4, "999")}[problem]
        fields[index] = value
        rows[rows.index(selected)] = "|".join(fields)
    with pytest.raises(ValueError, match="prerequisite"):
        require_completed("\n".join(rows), True)
