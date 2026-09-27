"""Validate an isolated reader dependency upgrade against retained readback."""

import argparse
import json
from pathlib import Path
import shlex
import subprocess

from packaging.requirements import Requirement
from packaging.utils import canonicalize_name

from benchmark_tools.audit_dependency_lock import evaluate
from benchmark_tools.audit_frozen_overlay_install import local_install_wheels
from benchmark_tools.audit_recovery_advisories import verify_inventory
from benchmark_tools.audit_recovery_install import installed_payload
from benchmark_tools.export_publication_readers import verify
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_publication_pipeline import save


def verify_lock(text, wheels):
    expected = []
    for line in text.splitlines():
        parts = shlex.split(line, comments=True)
        if not parts:
            continue
        if len(parts) != 2 or not parts[1].startswith("--hash=sha256:"):
            raise ValueError("Expected one exact requirement and SHA256 per line")
        requirement = Requirement(parts[0])
        specs = list(requirement.specifier)
        if (len(specs) != 1 or specs[0].operator != "==" or "*" in specs[0].version
                or requirement.url or requirement.marker or requirement.extras):
            raise ValueError("Requirement is not an unconditional exact pin")
        expected.append((canonicalize_name(requirement.name), specs[0].version, parts[1][14:]))
    actual = [(canonicalize_name(w["name"]), w["version"], w["wheel"]["sha256"]) for w in wheels]
    if not expected or len({r[0] for r in expected}) != len(expected) or sorted(expected) != sorted(actual):
        raise ValueError("Installed wheels differ from reader lock")


def compare_reports(previous, current):
    rows = []
    for name in ("structure", "sequences", "events", "hierarchy", "result"):
        old, new = previous / (name + ".json"), current / (name + ".json")
        left, right = json.loads(old.read_text()), json.loads(new.read_text())
        versions = {}
        if name == "sequences":
            versions = dict(previous=left.pop("biopython_version"), current=right.pop("biopython_version"))
        for report, root in ((left, previous), (right, current)):
            for item in report.get("checked_records", []):
                check(item)
                path = Path(item["path"])
                if path.parent == root:
                    if path.name not in {"structure.json", "events.json"}:
                        raise ValueError("Unexpected internal report reference")
                    # Each referenced report is independently compared in this loop.
                    item.clear()
                    item["compared_report"] = path.name
        # Only the aggregate's pointers to newly written reports may differ.
        if name == "result":
            for report, root in ((left, previous), (right, current)):
                if set(report["reports"]) != {"structure", "sequences", "events", "hierarchy"}:
                    raise ValueError("Unexpected reader report inventory")
                for key, item in report["reports"].items():
                    if item != record(root / (key + ".json")):
                        raise ValueError("Changed report reference")
                del report["reports"]
        if left != right:
            raise ValueError("Reader report changed: " + name)
        rows.append(dict(name=name, previous=record(old), current=record(new), equal=True,
                         biopython_versions=versions))
    return rows


def validate(runtime, readers, native, previous, snapshot, lock, output_directory):
    readback = output_directory / "readback"
    output_directory.mkdir(parents=True, exist_ok=False)
    if readback.exists():
        raise FileExistsError(readback)
    manifest = verify(readers)
    watched = [record(p) for p in (snapshot, lock, readers / "manifest.json", runtime / "install_report.json")]
    alerts = json.loads(snapshot.read_text())
    if (alerts["repository"] != "JLSteenwyk/orthohmm" or alerts["state_filter"] != "open"
            or not alerts["alerts"] or any(a["state"] != "open" for a in alerts["alerts"])):
        raise ValueError("Unexpected advisory snapshot")
    installation = json.loads((runtime / "install_report.json").read_text())
    wheels = local_install_wheels(installation, runtime / "wheels")
    verify_lock(lock.read_text(), wheels)
    python = str(runtime / "venv/bin/python")
    probe = "import json,sysconfig,importlib.metadata as m;print(json.dumps(dict(site=sysconfig.get_path('purelib'),packages=[dict(name=d.metadata['Name'],version=d.version) for d in m.distributions()])))"
    observed = json.loads(subprocess.check_output([python, "-I", "-c", probe], text=True, timeout=60))
    inventory = verify_inventory(installation, observed["packages"])
    site = Path(observed["site"])
    if not site.resolve().is_relative_to((runtime / "venv").resolve()):
        raise ValueError("Site-packages escaped reader environment")
    pip_check = subprocess.check_output([python, "-I", "-m", "pip", "check"], text=True, timeout=60)
    before = [installed_payload(Path(w["wheel"]["path"]), site) for w in wheels]
    save(output_directory / "payload_before.json", before)
    trace = output_directory / "readers.strace"
    code = ("import sys;from pathlib import Path;sys.path.insert(0," + repr(str(readers)) + ");"
            "from benchmark_tools.audit_publication_pipeline import audit;"
            "print(audit(Path(" + repr(str(native)) + "),Path(" + repr(str(readback)) + "))['status'])")
    command = ["/usr/bin/strace", "-f", "-s", "4096", "-e", "trace=%file", "-o", str(trace),
               python, "-I", "-B", "-c", code]
    environment = dict(HOME="/tmp", PATH="/usr/bin:/bin", LANG="C.UTF-8", OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
    output = subprocess.check_output(command, env=environment, cwd="/tmp", text=True, timeout=300)
    if output.strip() != "canonical_full_pipeline_scientific_readback_verified":
        raise ValueError("Readback did not succeed")
    comparisons = compare_reports(previous, readback)
    after = [installed_payload(Path(w["wheel"]["path"]), site) for w in wheels]
    save(output_directory / "payload_after.json", after)
    if before != after or verify(readers) != manifest:
        raise ValueError("Payload or reader source changed during execution")
    for item in watched:
        check(item)
    advisory_rows = evaluate(inventory, alerts["alerts"])
    return dict(status="reader_upgrade_fixture_equivalence_verified", command=command, environment=environment,
        inventory=inventory, pip_check=pip_check, wheels=wheels, comparisons=comparisons,
        checked_records=watched, source=record(__file__), trace=record(trace),
        payload_before=record(output_directory / "payload_before.json"), payload_after=record(output_directory / "payload_after.json"),
        matched_payload_files=sum(r["matched_files"] for r in before), advisory_comparisons=advisory_rows,
        affected_alerts=sum(r["status"] == "affected" for r in advisory_rows),
        inference_rerun=False, full_dataset_equivalence=False, comprehensive_security_clearance=False,
        publication_ready=False, limitations=["Same-host 16-gene fixture; no satellite merges or full-dataset equivalence.",
            "Version-range check covers only the supplied repository advisory snapshot, not exploitability or all dependencies.",
            "Payload checks exclude generated RECORD, bytecode and relocated non-site wheel data.",
            "The recorded file trace is not a sandbox; base Python and OS remain shared.",
            "Historical installations and lockfiles remain unchanged."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("runtime", "readers", "native", "previous", "snapshot", "lock", "output-directory", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = validate(*(getattr(args, n).resolve() for n in
                        ("runtime", "readers", "native", "previous", "snapshot", "lock", "output_directory")))
    save(args.output, result)
    print(json.dumps({k: result[k] for k in ("status", "matched_payload_files", "affected_alerts")}))
