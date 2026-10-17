"""Verify supporting protocol coverage, explicit checks and retained artifact identities."""

import collections
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import xml.etree.ElementTree as ET


def record(path):
    path = Path(path).resolve()
    before = path.stat()
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1048576), b""):
            digest.update(block)
    after = path.stat()
    assert (before.st_ino, before.st_size, before.st_mtime_ns) == (after.st_ino, after.st_size, after.st_mtime_ns)
    return {"path": str(path), "bytes": after.st_size, "sha256": digest.hexdigest()}


def identity(pin):
    return str(Path(pin["path"]).resolve()), pin["bytes"], pin["sha256"]


def follow(value, path):
    for key in path:
        value = value[key]
    return value


def selected_module(classname, modules):
    matches = [module for module in modules if classname == module or classname.startswith(module + ".")]
    if len(matches) != 1:
        raise ValueError("JUnit case must belong to exactly one selected test module: " + classname)
    return matches[0]


def collect(value, found):
    if isinstance(value, dict):
        if ({"path", "bytes", "sha256"} <= value.keys() and isinstance(value["path"], str)
                and type(value["bytes"]) is int and isinstance(value["sha256"], str) and len(value["sha256"]) == 64):
            found[identity(value)] = {key: value[key] for key in ("path", "bytes", "sha256")}
        for key, child in value.items():
            if key != "retained_details":
                collect(child, found)
    elif isinstance(value, list):
        for child in value:
            collect(child, found)


def historical_source(pin, repo):
    path = Path(pin["path"]).resolve()
    try:
        relative = path.relative_to(repo)
    except ValueError:
        return None
    if not (relative.parent == Path("benchmark_tools") or str(relative).startswith("orthohmm/")):
        return None
    if path.suffix not in (".py", ".c", ".json"):
        return None
    commits = subprocess.check_output(["git", "log", "--format=%H", "--", str(relative)], text=True).splitlines()
    for commit in commits:
        observed = subprocess.run(["git", "show", commit + ":" + str(relative)], capture_output=True)
        if (observed.returncode == 0 and len(observed.stdout) == pin["bytes"]
                and hashlib.sha256(observed.stdout).hexdigest() == pin["sha256"]):
            return {"revision": commit, "path_in_git": str(relative),
                "scope": "Retained historical source/protocol bytes in Git, not a match to the current mutable path or retrospective execution attestation"}
    return None


def main(root, destination):
    if destination.exists():
        raise FileExistsError("Preserve previous evidence reviews")
    repo = Path.cwd().resolve()
    spec_path = root / "supporting_review_specification.json"
    spec = json.loads(spec_path.read_text())
    inventory_pin = spec["inventory"]
    assert identity(record(inventory_pin["path"])) == identity(inventory_pin)
    inventory = json.loads(Path(inventory_pin["path"]).read_text())
    assert len(inventory["protocol_inventory"]) == 67
    assigned = [i for group in spec["groups"] for i in group["protocol_indices"]]
    assert sorted(assigned) == list(range(67)), "Every protocol must be reviewed exactly once"
    documents, pins, groups = {}, {}, []
    for group in spec["groups"]:
        items, receipts = [], []
        for name in group["receipts"]:
            path = (Path(spec["receipt_root"]) / name).resolve()
            if str(path) not in documents:
                pin = record(path)
                document = json.loads(path.read_text())
                assert record(path) == pin
                documents[str(path)] = (pin, document)
                collect(document, pins)
                pins[identity(pin)] = pin
            pin, document = documents[str(path)]
            items.append(document)
            receipts.append({"source": pin, "status": document.get("status")})
        for check in group.get("checks", []):
            observed = follow(items[check["receipt"]], check["path"])
            assert type(observed) is type(check["equals"]) and observed == check["equals"], (group["name"], check)
        protocols = []
        for index in group["protocol_indices"]:
            protocol = inventory["protocol_inventory"][index]
            assert identity(record(protocol["source"]["path"])) == identity(protocol["source"])
            protocols.append({"index": index, "source": protocol["source"],
                "first_committed_revision": protocol["first_committed_revision"]})
        groups.append({**group, "protocols": protocols, "receipts": receipts,
            "disposition": "reviewed_with_stated_scope_and_conditional_stops"})
    regression_path = repo / "benchmarks/results/hillclimb_campaign_regression_v2_20261016/selection.json"
    regression_pin = record(regression_path)
    regression = json.loads(regression_path.read_text())
    samples = regression["primary_samples_before"]
    for pin in samples["file_records"]:
        pins[identity(pin)] = pin
    original_files = {str(Path(p["path"]).resolve()) for p in samples["file_records"]}
    assert {str(p.resolve()) for p in (repo / "tests/samples").rglob("*") if p.is_file()} == original_files
    fingerprints = {}
    for key, command in (("status_sha256", ["git", "status", "--short", "--", "tests/samples"]),
            ("tracked_diff_sha256", ["git", "diff", "--binary", "--", "tests/samples"])):
        value = subprocess.check_output(command)
        fingerprints[key] = hashlib.sha256(value).hexdigest()
        assert fingerprints[key] == samples[key]
    assert subprocess.check_output(["git", "diff", "de89fd84..HEAD", "--", "orthohmm"]) == b""
    primary_xml = repo / "benchmarks/results/hillclimb_campaign_regression_v2_20261016/regression.xml"
    suites = list(ET.parse(primary_xml).getroot().iter("testsuite"))
    primary_counts = {key: sum(int(s.attrib.get(key, 0)) for s in suites)
                      for key in ("tests", "failures", "errors", "skipped")}
    assert primary_counts == dict(tests=2349, failures=0, errors=0, skipped=0)
    modules = {path[:-3].replace("/", ".") for path in regression["selected_paths"]}
    executed = collections.Counter()
    for case in ET.parse(primary_xml).getroot().iter("testcase"):
        module = selected_module(case.attrib["classname"], modules)
        if case.find("skipped") is None:
            executed[module] += 1
    assert len(executed) == 99 and set(executed) == modules and all(executed.values())
    reader_xml = repo / "benchmarks/results/hillclimb_supporting_reader_tests_20261016.xml"
    reader_suites = list(ET.parse(reader_xml).getroot().iter("testsuite"))
    reader_counts = {key: sum(int(s.attrib.get(key, 0)) for s in reader_suites)
                     for key in ("tests", "failures", "errors", "skipped")}
    assert reader_counts == dict(tests=78, failures=0, errors=0, skipped=0)
    current, dispositions, unresolved, xml = {}, [], [], []
    by_digest = collections.defaultdict(list)
    for pin in pins.values():
        by_digest[(pin["bytes"], pin["sha256"])].append(pin)
    for item in pins.values():
        path = str(Path(item["path"]).resolve())
        if path not in current:
            current[path] = record(path) if Path(path).is_file() else None
        actual = current[path]
        if actual is not None and identity(actual) == identity(item):
            dispositions.append({"source": item, "verification": "current_bytes_match"})
            if Path(path).suffix == ".xml":
                suites = list(ET.parse(path).getroot().iter("testsuite"))
                xml.append({"source": item, "counts": {key: sum(int(s.attrib.get(key, 0)) for s in suites)
                    for key in ("tests", "failures", "errors", "skipped")}})
            continue
        historical = historical_source(item, repo)
        if historical is not None:
            dispositions.append({"source": item, "verification": "historical_source_git_binding", **historical})
            continue
        archive = None
        for candidate in by_digest[(item["bytes"], item["sha256"])]:
            archived_path = str(Path(candidate["path"]).resolve())
            if archived_path == path:
                continue
            if archived_path not in current:
                current[archived_path] = record(archived_path) if Path(archived_path).is_file() else None
            if current[archived_path] is not None and identity(current[archived_path]) == identity(candidate):
                archive = candidate
                break
        if archive is not None:
            dispositions.append({"source": item, "verification": "retained_exact_archive_binding", "archive": archive})
        else:
            unresolved.append({"expected": item, "actual": actual})
    for pin, _ in documents.values():
        assert record(pin["path"]) == pin
    assert record(regression_path) == regression_pin
    report = dict(status="supporting_protocol_evidence_review_verified" if not unresolved else "supporting_protocol_evidence_bindings_pending",
        specification=record(spec_path), recipe=record(__file__), inventory=inventory_pin, groups=groups,
        protocol_count=67, group_count=len(groups), unique_receipts=len(documents),
        artifact_dispositions=dispositions, unresolved=unresolved, parsed_test_reports=xml,
        primary_regression={"source": record(primary_xml), "counts": primary_counts, "executed_modules": len(executed)},
        reader_tests={"source": record(reader_xml), "counts": reader_counts,
            "new_regression_test": record(repo / "tests/unit/test_supporting_campaign_readback.py")},
        primary_worktree_preserved={"selection": regression_pin, "sample_files": len(original_files), "fingerprints": fingerprints},
        failed_attempt={"job": "24387 FAILED 1:0", "recipe_revision": "8c84427b",
            "log": record("benchmarks/results/hillclimb_supporting_audit_20261016-24387.out"),
            "reason": "Reader counted 107 JUnit classes as modules; all cases map uniquely to the exact 99 selected modules. Original XML still has 2349 passes and zero failures/errors/skips. Corrected module preflight now precedes artifact hashing; original failed result namespace has no output."},
        actual_files_hashed=len(current), actual_bytes_hashed=sum(p["bytes"] for p in current.values() if p),
        limitations=["Supporting receipt gates, direct artifact identities and historical source/archive bindings are explicit; Git bindings do not retroactively prove historical execution or current-file equality.",
            "The already audited 24 retained comparator cells and 480 registered identities are inherited rather than rehashing unrelated OrthoMCL BPO files through retained_details. Their historical input-consumption limits remain explicit.",
            "Per-protocol findings and untriggered conditional requirements retain the original stopping rules. No raw scoring, inference, new accuracy result or isolated speed ranking is created.",
            "Referenced historical failed/skipped test reports are retained with their counts; they are not silently counted as passes or summed into independent replicate totals."],
        scientific_inference_rerun=False, defaults_promoted=False, goal_completed=False)
    with destination.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"status": report["status"], "protocols": 67, "groups": len(groups),
        "receipts": len(documents), "dispositions": len(dispositions), "unresolved": len(unresolved),
        "actual_files_hashed": len(current), "actual_GiB_hashed": report["actual_bytes_hashed"] / 2**30,
        "historical_bindings": [d["source"]["path"] for d in dispositions if d["verification"] != "current_bytes_match"]}))


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
