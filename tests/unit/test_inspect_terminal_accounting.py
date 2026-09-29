import pytest

from benchmark_tools.inspect_terminal_accounting import FIELDS, assess


def text(memory="", cpu="00:00:00"):
    return FIELDS.replace(",", "|") + f"\n123|COMPLETED|0:0|10|64|{cpu}|{memory}||start|end\n"


@pytest.mark.parametrize("plugin", [None, "(null)", "jobacct_gather/none"])
def test_disabled_is_not_zero_usage(plugin):
    result = assess(text(), [123], plugin)
    assert result["unavailable_cpu_rows"] == result["missing_memory_rows"] == 1
    assert result["complete_job_resource_accounting_verified"] is False


def test_available_fields_never_automatically_admit_resources():
    result = assess(text("123K", "00:00:03"), [123], "jobacct_gather/cgroup")
    assert result["unavailable_cpu_rows"] == result["missing_memory_rows"] == 0
    assert result["complete_job_resource_accounting_verified"] is False


def test_missing_memory_is_unavailable_with_enabled_collector():
    assert assess(text(), [123], "jobacct_gather/cgroup")["missing_memory_rows"] == 1


@pytest.mark.parametrize("bad,jobs", [("bad", [123]), (text(), [124]), (text(), [123, 124]),
    (text() + text().splitlines()[1] + "\n", [123])])
def test_bad_inventories_rejected(bad, jobs):
    with pytest.raises(ValueError):
        assess(bad, jobs, "(null)")
