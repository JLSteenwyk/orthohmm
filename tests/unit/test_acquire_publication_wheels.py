import io
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmark_tools import acquire_publication_wheels as module


@pytest.fixture
def preparation(tmp_path, monkeypatch):
    rows, payloads, urls = [], {}, {}
    for name in sorted(set.union(*module.ROLES.values())):
        filename = name.replace("-", "_") + "-1.0.0-py3-none-any.whl"
        path = tmp_path / filename
        payload = ("Synthetic wheel: " + name).encode()
        path.write_bytes(payload)
        row = dict(name=name, version="1.0.0", wheel=module.record(path))
        rows.append(row)
        wheel_url = "https://files.pythonhosted.org/packages/fixed/" + filename
        url = f"https://pypi.org/pypi/{name}/1.0.0/json"
        metadata = dict(info=dict(name=name, version="1.0.0"), urls=[dict(filename=filename,
            size=len(payload), digests=dict(sha256=row["wheel"]["sha256"]), yanked=False,
            packagetype="bdist_wheel", url=wheel_url)])
        payloads[url], payloads[wheel_url] = json.dumps(metadata).encode(), payload
        urls[name] = url
    inventory = tmp_path / "inventory.json"
    inventory.write_text(json.dumps(dict(status="selected_wheel_elf_inventory", wheels=rows)))
    monkeypatch.setattr(module, "INVENTORY_SHA", module.record(inventory)["sha256"])
    locks = {}
    for role in module.ROLES:
        path = tmp_path / (role + ".txt")
        path.write_text("# Synthetic exact lock, never interpreted or installed\n")
        locks[role] = path
    monkeypatch.setattr(module, "LOCK_SHAS", {role: module.record(p)["sha256"] for role, p in locks.items()})
    supplied = {name: Path(next(r["wheel"]["path"] for r in rows if r["name"] == name)) for name in ("pip", "orthohmm")}
    args = SimpleNamespace(inventory=inventory, inference_lock=locks["inference"], reader_lock=locks["reader"],
        orthohmm_wheel=supplied["orthohmm"], pip_wheel=supplied["pip"], output=tmp_path / "fresh output",
        acknowledge_historical_runtime=True, timeout=60)
    calls = []

    class Response(io.BytesIO):
        def __init__(self, payload, url):
            super().__init__(payload)
            self.headers, self.url = {"Content-Length": str(len(payload))}, url
        def geturl(self):
            return self.url

    class Opener:
        def open(self, request, timeout):
            calls.append(request.full_url)
            assert timeout == 60
            return Response(payloads[request.full_url], request.full_url)

    monkeypatch.setattr(module.network, "network_opener", lambda: Opener())
    return args, calls, payloads, urls


def test_two_exact_wheel_sets_without_installation(preparation):
    args, calls, _, urls = preparation
    result = module.acquire(args)
    assert result["status"] == "historical_inference_reader_wheels_acquired"
    assert len(calls) == 20 and len(set(calls)) == 20
    assert urls["pip"] not in calls and urls["orthohmm"] not in calls
    assert (result["union_wheels"], result["public_wheels"], result["supplied_wheels"]) == (12, 10, 2)
    assert {k: len(v["wheels"]) for k, v in result["environments"].items()} == dict(inference=11, reader=5)
    for role in module.ROLES:
        assert (args.output / (role + "_requirements.txt")).read_bytes() == getattr(args, role + "_lock").read_bytes()
    assert result == json.loads((args.output / "complete.json").read_bytes())
    assert all(result[k] is False for k in ("installation_performed", "native_inference_executed",
        "publication_ready", "security_clearance", "redistribution_clearance", "retry"))


@pytest.mark.parametrize("defect", ["inventory", "lock", "pip", "project", "filename", "symlink", "ack", "timeout", "existing"])
def test_preflight_before_output_and_network(preparation, tmp_path, defect):
    args, calls, _, _ = preparation
    if defect in {"inventory", "lock", "pip", "project"}:
        path = {"inventory": args.inventory, "lock": args.reader_lock, "pip": args.pip_wheel,
                "project": args.orthohmm_wheel}[defect]
        path.write_text("changed")
    elif defect == "filename":
        alias = tmp_path / "wrong.whl"
        alias.write_bytes(args.pip_wheel.read_bytes())
        args.pip_wheel = alias
    elif defect == "symlink":
        alias = tmp_path / "wheel-link"
        alias.symlink_to(args.pip_wheel)
        args.pip_wheel = alias
    elif defect == "ack":
        args.acknowledge_historical_runtime = False
    elif defect == "timeout":
        args.timeout = True
    else:
        args.output.mkdir()
    with pytest.raises((ValueError, FileExistsError)):
        module.acquire(args)
    assert not calls
    if defect != "existing":
        assert not args.output.exists()


@pytest.mark.parametrize("defect", ["version", "name", "duplicate", "size", "hash", "yanked", "host", "basename", "kind"])
def test_metadata_must_identify_exact_wheel(preparation, defect):
    args, _, payloads, urls = preparation
    wheels = module.wheel_rows(json.loads(args.inventory.read_bytes()))
    wheel = wheels["numpy"]
    metadata = json.loads(payloads[urls["numpy"]])
    row = metadata["urls"][0]
    if defect in {"version", "name"}:
        metadata["info"][defect] = "different"
    elif defect == "duplicate":
        metadata["urls"].append(dict(row))
    elif defect == "size":
        row["size"] += 1
    elif defect == "hash":
        row["digests"]["sha256"] = "0" * 64
    elif defect == "yanked":
        row["yanked"] = True
    elif defect == "host":
        row["url"] = "https://pypi.org/" + wheel["filename"]
    elif defect == "basename":
        row["url"] = "https://files.pythonhosted.org/different.whl"
    else:
        row["packagetype"] = "sdist"
    with pytest.raises(ValueError):
        module.public_selection(metadata, wheel)


@pytest.mark.parametrize("defect", ["network", "metadata", "wheel", "supplied_changed", "copy"])
def test_failure_retained_without_retry(preparation, monkeypatch, defect):
    args, calls, payloads, urls = preparation
    if defect == "metadata":
        payloads[urls["biopython"]] = b"{}"
    elif defect == "wheel":
        url = json.loads(payloads[urls["biopython"]])["urls"][0]["url"]
        payloads[url] = b"x" * len(payloads[url])
    elif defect == "network":
        monkeypatch.setattr(module.network, "network_opener", lambda: (_ for _ in ()).throw(OSError("network failure")))
    elif defect == "supplied_changed":
        copy = module.copy_wheel
        def change(source, target, expected):
            result = copy(source, target, expected)
            if Path(source) == args.pip_wheel:
                args.pip_wheel.write_text("changed during copy")
            return result
        monkeypatch.setattr(module, "copy_wheel", change)
    else:
        monkeypatch.setattr(module, "copy_wheel", lambda *args: (_ for _ in ()).throw(OSError("copy failure")))
    with pytest.raises((OSError, KeyError, ValueError)):
        module.acquire(args)
    failed = json.loads((args.output / "failed.json").read_bytes())
    assert failed["retry"] is False and failed["attempts"] == 1
    assert not (args.output / "complete.json").exists()
    assert len(set(calls)) == len(calls)
    with pytest.raises(FileExistsError):
        module.acquire(args)


@pytest.mark.parametrize("defect", ["duplicate", "missing", "size", "hash", "filename", "name"])
def test_retained_inventory_contract(preparation, defect):
    args, _, _, _ = preparation
    data = json.loads(args.inventory.read_bytes())
    row = data["wheels"][0]
    if defect == "duplicate":
        data["wheels"].append(dict(row))
    elif defect == "missing":
        data["wheels"].pop()
    elif defect == "size":
        row["wheel"]["bytes"] = True
    elif defect == "hash":
        row["wheel"]["sha256"] = "bad"
    elif defect == "filename":
        row["wheel"]["path"] = "wrong file.whl"
    else:
        row["name"] = "unknown"
    with pytest.raises(ValueError):
        module.wheel_rows(data)
