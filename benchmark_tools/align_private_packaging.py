"""Align the private candidate to hash-proven historical packaging source."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import zipfile

from benchmark_tools.build_private_timing_environment import record, write, wheel_metadata


def compare_payload(wheel, records, old_root):
    expected = {str(Path(r["path"]).relative_to(old_root)): r for r in records
                if r["kind"] == "file" and Path(r["path"]).is_relative_to(old_root)}
    if not expected:
        raise ValueError("Missing historical packaging payload")
    with zipfile.ZipFile(wheel) as archive:
        names = [n for n in archive.namelist() if n.startswith("packaging/") and not n.endswith("/")]
        if len(names) != len(set(names)) or {n.removeprefix("packaging/") for n in names} != set(expected):
            raise ValueError("Packaging wheel member inventory differs")
        checked = []
        for name in names:
            raw = archive.read(name)
            prior = expected[name.removeprefix("packaging/")]
            digest = hashlib.sha256(raw).hexdigest()
            if digest != prior["sha256"] or len(raw) != prior["bytes"]:
                raise ValueError("Packaging wheel does not match historical source")
            checked.append(dict(member=name, bytes=len(raw), sha256=digest, prior_path=prior["path"]))
    return checked


def run(candidate, wheel, historical, expected_sha, output):
    if output.exists():
        raise FileExistsError(output)
    raw = historical.read_bytes()
    if hashlib.sha256(raw).hexdigest() != expected_sha:
        raise ValueError("Historical runtime inventory changed")
    if wheel_metadata(wheel) != ("packaging", "26.1"):
        raise ValueError("Require packaging 26.1 wheel")
    checked = compare_payload(wheel, json.loads(raw)["records"],
        Path("/home/bizon/anaconda3/lib/python3.10/site-packages/packaging"))
    installation = candidate / "result.json"
    prior = json.loads(installation.read_text())
    if prior["status"] != "private_timing_environment_candidate_installed" or prior["selected"]["packaging"] != "26.0":
        raise ValueError("Wrong private candidate")
    output.mkdir(parents=True)
    python = candidate / "venv/bin/python"
    env = json.loads((candidate / "imports_started.json").read_text())["environment"]
    env.update(NUMBA_CACHE_DIR=str(output / "numba_cache"))
    Path(env["NUMBA_CACHE_DIR"]).mkdir()
    wheel_record = record(wheel)
    lock = output / "packaging.txt"
    lock.write_text(f"packaging==26.1 --hash=sha256:{wheel_record['sha256']}\n")
    write(output / "started.json", dict(source=record(Path(__file__)), previous=record(installation),
        historical=record(historical), wheel=wheel_record, checked_payload=checked))
    stages = []
    def execute(name, command):
        write(output / f"{name}_started.json", dict(command=command, environment=env))
        with (output / f"{name}.log").open("x") as log:
            result = subprocess.run(command, cwd=output, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=180)
        item = dict(command=command, exit_code=result.returncode, log=record(output / f"{name}.log"))
        write(output / f"{name}_finished.json", item)
        stages.append(item)
        if result.returncode:
            raise RuntimeError(f"{name} failed; preserve attempt")
    execute("install", [str(python), "-I", "-B", "-m", "pip", "--isolated", "install",
        "--no-index", "--no-deps", "--only-binary=:all:", "--require-hashes", "--find-links", str(wheel.parent), "-r", str(lock)])
    execute("dependencies", [str(python), "-I", "-B", "-m", "pip", "check"])
    command = json.loads((candidate / "imports_started.json").read_text())["command"]
    execute("imports", command)
    observed = json.loads((output / "imports.log").read_text())
    selected = dict(prior["selected"], packaging="26.1")
    if observed["packages"].get("packaging") != "26.1":
        raise ValueError("Private package version not aligned")
    installed_root = Path(observed["modules"]["packaging"]).parent
    installed = []
    for row in checked:
        item = record(installed_root / row["member"].removeprefix("packaging/"))
        if item["sha256"] != row["sha256"] or item["bytes"] != row["bytes"]:
            raise ValueError("Installed packaging payload differs")
        installed.append(item)
    result = dict(status="private_packaging_historical_payload_aligned", selected=selected,
        previous=record(installation), baseline=prior["baseline"], source=record(Path(__file__)),
        historical=record(historical), wheel=wheel_record, lock=record(lock), stages=stages,
        checked_payload=checked, installed=installed, import_report=record(output / "imports.log"),
        scientific_execution_authorized=False,
        limitations=["Only the private candidate changed; previous installation and comparison receipts remain historical.",
                     "Historical module bytes, not stale distribution metadata, establish packaging 26.1 selection.",
                     "PyYAML build difference and native integration remain unvalidated; no timing admission."])
    write(output / "result.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("candidate", "wheel", "historical", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--sha256", required=True)
    args = parser.parse_args()
    result = run(args.candidate.resolve(), args.wheel.resolve(), args.historical.resolve(), args.sha256, args.output.resolve())
    print(json.dumps(dict(status=result["status"], files=len(result["installed"])) ))
