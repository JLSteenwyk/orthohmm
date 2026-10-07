"""CPU selection uses synthetic topology; no affinity or scheduler mutation."""

from copy import deepcopy

import pytest

from benchmark_tools import native_factorial_cpu_selection as module


@pytest.fixture
def case():
    values = list(range(40, 72)) + list(range(136, 168))
    topology = [dict(cpu=cpu, package=0, core=(cpu - 40) % 96, numa_node=1) for cpu in values]
    return values, topology


def test_chooses_allocated_physical_cores_without_binding(case):
    result = module.choose(*case)
    assert result["native_cpu_ids"] == list(range(40, 72))
    assert result["native_physical_core_count"] == 32
    assert result["allocated_logical_cpu_count"] == 64
    assert result["historical_fixed_mask_compatible"] is False
    assert int(result["native_os_cpu_mask"], 16).bit_count() == 32
    for key in ("affinity_changed", "scientific_execution_authorized", "scientific_timings_admitted",
                "publication_ready"):
        assert result[key] is False


def test_topology_order_does_not_change_selection(case):
    values, topology = case
    assert module.choose(values, topology) == module.choose(values, list(reversed(topology)))


def test_original_mask_is_retained_when_allocated():
    values = list(range(32)) + list(range(96, 128))
    topology = [dict(cpu=cpu, package=0, core=cpu % 96) for cpu in values]
    result = module.choose(values, topology)
    assert result["historical_fixed_mask_compatible"] is True
    assert result["native_cpu_ids"] == list(range(32))
    assert result["native_os_cpu_mask"] == "0xffffffff"


def test_first_32_logical_ids_are_not_assumed_to_be_physical_cores():
    values = list(range(64))
    topology = [dict(cpu=cpu, package=0, core=cpu // 2) for cpu in values]
    result = module.choose(values, topology)
    assert result["native_cpu_ids"] == list(range(0, 64, 2))
    assert result["historical_fixed_mask_compatible"] is False


@pytest.mark.parametrize("change", ["fewer", "more", "duplicate", "unsorted", "negative", "outside",
    "bool", "tuple", "missing_topology", "duplicate_topology", "unallocated", "package", "core",
    "extra_field", "missing_field", "numa", "numa_bool", "sibling_numa", "fewer_physical",
    "more_physical", "non_sibling"])
def test_refuses_unmatched_or_invalid_allocations(case, change):
    values, topology = deepcopy(case)
    if change == "fewer": values.pop()
    elif change == "more": values.append(168)
    elif change == "duplicate": values[-1] = values[0]
    elif change == "unsorted": values.reverse()
    elif change == "negative": values[0] = -1
    elif change == "outside": values[-1] = 192
    elif change == "bool": values[0] = True
    elif change == "tuple": values = tuple(values)
    elif change == "missing_topology": topology.pop()
    elif change == "duplicate_topology": topology[-1] = dict(topology[0])
    elif change == "unallocated": topology[-1]["cpu"] = 100
    elif change == "package": topology[-1]["package"] = -1
    elif change == "core": topology[-1]["core"] = True
    elif change == "extra_field": topology[-1]["ignored"] = 1
    elif change == "missing_field": topology[-1].pop("core")
    elif change == "numa": topology[-1]["numa_node"] = -1
    elif change == "numa_bool": topology[-1]["numa_node"] = True
    elif change == "sibling_numa": topology[-1]["numa_node"] = 2
    elif change == "fewer_physical":
        for row in topology: row["core"] %= 16
    elif change == "more_physical":
        for row in topology: row["core"] = row["cpu"]
    else: topology[-1]["package"] = 1
    with pytest.raises(ValueError):
        module.choose(values, topology)


def write_topology(root, values):
    for cpu in values:
        directory = root / f"cpu{cpu}"
        (directory / "topology").mkdir(parents=True)
        (directory / "topology/physical_package_id").write_text("0\n")
        (directory / "topology/core_id").write_text(str((cpu - 40) % 96))
        (directory / "online").write_text("1\n")
        (directory / "node1").mkdir()


def test_host_topology_reads_recorded_fields(tmp_path, case):
    values, topology = case
    write_topology(tmp_path, values)
    assert module.host_topology(values, tmp_path) == topology


@pytest.mark.parametrize("change", ["offline", "ambiguous", "missing", "nonnumeric"])
def test_refuses_invalid_host_topology(tmp_path, case, change):
    values, _ = case
    write_topology(tmp_path, values)
    directory = tmp_path / "cpu40"
    if change == "offline": (directory / "online").write_text("0\n")
    elif change == "ambiguous": (directory / "node2").mkdir()
    elif change == "missing": (directory / "topology/core_id").unlink()
    else: (directory / "topology/core_id").write_text("unknown")
    with pytest.raises((ValueError, FileNotFoundError)):
        module.host_topology(values, tmp_path)
