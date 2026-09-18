import json
from pathlib import Path

import pytest

from benchmark_tools.audit_fastoma_resource_probe import IMAGE_DIGEST, validate_probe, validate_config


def fixture():
    return ["name\tstatus\texit\nresource_probe\tCOMPLETED\t0\n",
            json.dumps({"cpu_max": "100000 100000", "memory_max": "268435456"}),
            f"docker run -i --cpus 1.0 --memory 256m --network none {IMAGE_DIGEST} /bin/bash task.sh\n", "0\n"]


def test_observed_limits_and_wrapper_agree():
    result = validate_probe(*fixture())
    assert result["observed_limits"]["memory_max"] == "268435456"
    assert IMAGE_DIGEST in result["docker_argv"]


@pytest.mark.parametrize("mutation", ["live", "cached", "failed", "extra_task", "unlimited_cpu",
    "two_cpus", "zero_period", "memory", "mutable_image", "host_network", "duplicate_cpu", "privileged"])
def test_invalid_probe_refused(mutation):
    args = fixture()
    if mutation == "live":
        args[0] = args[0].replace("COMPLETED", "RUNNING")
    elif mutation == "cached":
        args[0] = args[0].replace("COMPLETED", "CACHED")
    elif mutation == "failed":
        args[3] = "1"
    elif mutation == "extra_task":
        args[0] += "resource_probe\tCOMPLETED\t0\n"
    elif mutation in ("unlimited_cpu", "two_cpus", "zero_period", "memory"):
        output = json.loads(args[1])
        if mutation == "memory":
            output["memory_max"] = "max"
        else:
            output["cpu_max"] = {"unlimited_cpu": "max 100000", "two_cpus": "200000 100000",
                                 "zero_period": "0 0"}[mutation]
        args[1] = json.dumps(output)
    elif mutation == "mutable_image":
        args[2] = args[2].replace(IMAGE_DIGEST, "dessimozlab/fastoma:0.3.5")
    elif mutation == "host_network":
        args[2] = args[2].replace("--network none", "--network host")
    elif mutation == "duplicate_cpu":
        args[2] = args[2].replace("--cpus 1.0", "--cpus 1.0 --cpus 1")
    else:
        args[2] = args[2].replace("docker run", "docker run --privileged")
    with pytest.raises(ValueError):
        validate_probe(*args)


def effective_config():
    path = Path(__file__).resolve().parents[2] / "benchmark_tools/results/fastoma_corrected_resource_probe_20260918.json"
    return json.loads(path.read_text())["effective_configuration"]


def test_real_effective_configuration():
    result = validate_config(effective_config())
    assert result["process.container"] == repr(IMAGE_DIGEST)
    assert result["params.nr_repr_per_hog"] == "5"


@pytest.mark.parametrize("key", ["process.container", "executor.cpus", "executor.memory", "params.max_cpus",
    "params.max_memory", "params.fasta_header_id_transformer", "params.force_pairwise_ortholog_generation",
    "params.filter_method", "params.filter_gap_ratio_row", "params.filter_gap_ratio_col", "params.nr_repr_per_hog",
    "params.min_sequence_length", "docker.enabled", "docker.runOptions", "trace.enabled",
    "process.'withName:omamer_run'.memory", "process.'withName:collect_subhogs'.memory",
    "process.'withName:extract_pairwise_ortholog_relations'.memory"])
def test_changed_configuration_refused(key):
    lines = effective_config().splitlines()
    text = "\n".join(key + " = 'changed'" if s.startswith(key + " = ") else s for s in lines)
    with pytest.raises(ValueError, match="configuration differs"):
        validate_config(text)


def test_duplicate_configuration_key_refused():
    with pytest.raises(ValueError, match="Duplicate"):
        validate_config(effective_config() + "\nexecutor.cpus = 180\n")
