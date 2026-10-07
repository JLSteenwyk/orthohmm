"""Choose 32 physical cores within a measured 64-slot step, without binding."""

from pathlib import Path


def cpu_ids(values):
    if (type(values) is not list or len(values) != 64
            or any(type(cpu) is not int or not 0 <= cpu < 192 for cpu in values)
            or len(set(values)) != 64 or values != sorted(values)):
        raise ValueError("Require sorted distinct OS IDs from an actual 64-slot step")
    return values


def host_topology(values, cpu_root=Path("/sys/devices/system/cpu")):
    records = []
    for cpu in cpu_ids(values):
        root = Path(cpu_root) / f"cpu{cpu}"
        online = root / "online"
        if online.exists() and online.read_text().strip() != "1":
            raise ValueError("Allocated CPU is offline")
        nodes = list(root.glob("node[0-9]*"))
        if len(nodes) > 1:
            raise ValueError("Ambiguous CPU NUMA node")
        records.append(dict(cpu=cpu,
            package=int((root / "topology/physical_package_id").read_text()),
            core=int((root / "topology/core_id").read_text()),
            numa_node=int(nodes[0].name[4:]) if nodes else None))
    return records


def choose(values, topology):
    values = cpu_ids(values)
    if type(topology) is not list or len(topology) != 64:
        raise ValueError("Require topology for every allocated OS CPU")
    by_cpu = {}
    groups = {}
    for row in topology:
        if (not isinstance(row, dict) or not {"cpu", "package", "core"} <= set(row)
                or not set(row) <= {"cpu", "package", "core", "numa_node"}
                or any(type(row[k]) is not int or row[k] < 0 for k in ("cpu", "package", "core"))
                or row["cpu"] not in values or row["cpu"] in by_cpu
                or row.get("numa_node") is not None
                and (type(row["numa_node"]) is not int or row["numa_node"] < 0)):
            raise ValueError("Invalid, duplicate or unallocated topology record")
        by_cpu[row["cpu"]] = row
        groups.setdefault((row["package"], row["core"]), []).append(row["cpu"])
    if set(by_cpu) != set(values) or len(groups) != 32 or any(len(group) != 2 for group in groups.values()):
        raise ValueError("Require exactly 32 physical cores with two allocated SMT slots each")
    records = []
    for (package, core), group in sorted(groups.items()):
        group = sorted(group)
        nodes = {by_cpu[cpu].get("numa_node") for cpu in group}
        if len(nodes) != 1:
            raise ValueError("SMT siblings disagree on NUMA node")
        records.append(dict(package=package, core=core, allocated_os_cpu_ids=group,
                            selected_os_cpu_id=group[0], numa_node=nodes.pop()))
    selected = sorted(row["selected_os_cpu_id"] for row in records)
    return dict(schema="native_factorial_cpu_selection_v1",
        strategy="lowest_os_cpu_id_per_allocated_physical_core",
        allocated_os_cpu_ids=values, allocated_logical_cpu_count=64, physical_cores=records,
        native_cpu_ids=selected, native_physical_core_count=32,
        native_os_cpu_mask=hex(sum(1 << cpu for cpu in selected)),
        historical_fixed_mask_compatible=selected == list(range(32)),
        affinity_changed=False, scientific_execution_authorized=False,
        scientific_timings_admitted=False, publication_ready=False,
        limitations=["Selection only; caller must obtain and bind the actual running step allocation.",
            "No scheduler submission, CPU binding, cgroup limit check, inference or timing is performed.",
            "A fresh step-placement, runtime, capacity and accounting check remains mandatory.",
            "Different CPU/NUMA placement must be recorded; matched counts do not establish isolated performance."])
