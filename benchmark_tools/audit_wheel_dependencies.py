"""Check declared wheel dependencies for a hash-pinned lock on this interpreter."""

import argparse
from email.parser import BytesParser
import json
from pathlib import Path, PurePosixPath
import platform
import zipfile

from packaging.markers import default_environment
from packaging.requirements import Requirement
from packaging.specifiers import SpecifierSet
from packaging.utils import canonicalize_name, parse_wheel_filename
from packaging.version import Version

from benchmark_tools.audit_pypi_releases import read_pins
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


def metadata(path, pin):
    ref = record(path.resolve())
    name, version, _, _ = parse_wheel_filename(path.name)
    if (name != pin["name"] or version != Version(pin["version"])
            or ref["sha256"] not in pin["hashes"]):
        raise ValueError("Wheel filename or hash differs from lock")
    with zipfile.ZipFile(path) as archive:
        candidates = [n for n in archive.namelist() if n.endswith(".dist-info/METADATA")
                      and len(PurePosixPath(n).parts) == 2]
        if len(candidates) != 1:
            raise ValueError("Require exactly one wheel METADATA")
        message = BytesParser().parsebytes(archive.read(candidates[0]))
    for key in ("Name", "Version", "Requires-Python"):
        if len(message.get_all(key, [])) > 1:
            raise ValueError("Duplicate singleton metadata field")
    if (canonicalize_name(message["Name"]) != pin["name"]
            or Version(message["Version"]) != version):
        raise ValueError("Wheel metadata identity differs from lock")
    return dict(name=pin["name"], version=str(version), wheel=ref,
                requires_python=message.get("Requires-Python"),
                requires_dist=message.get_all("Requires-Dist", []))


def evaluate(packages, environment):
    if len({p["name"] for p in packages}) != len(packages):
        raise ValueError("Duplicate package metadata")
    versions = {p["name"]: p["version"] for p in packages}
    rows, failures = [], []
    for package in packages:
        spec = package["requires_python"]
        if spec and not SpecifierSet(spec).contains(environment["python_full_version"], prereleases=True):
            failures.append(dict(package=package["name"], reason="requires_python", requirement=spec))
        for declaration in package["requires_dist"]:
            req = Requirement(declaration)
            active = req.marker is None or req.marker.evaluate(dict(environment, extra=""))
            target = canonicalize_name(req.name)
            installed = versions.get(target)
            status = "inactive_marker"
            if active:
                status = ("unsupported_url_or_extras" if req.url or req.extras else
                          "missing" if installed is None else
                          "satisfied" if req.specifier.contains(installed, prereleases=True) else "version_mismatch")
                if status != "satisfied":
                    failures.append(dict(package=package["name"], reason=status, requirement=declaration))
            rows.append(dict(package=package["name"], requirement=declaration,
                             target=target, locked_version=installed, active=active, status=status))
    return dict(dependency_rows=rows, failures=failures, declared_runtime_dependencies_satisfied=not failures)


def audit(lock, wheelhouse, output):
    if output.exists():
        raise FileExistsError(output)
    lock_ref = record(lock.resolve())
    pins = read_pins(lock.read_text())
    packages = []
    wheels = list(wheelhouse.glob("*.whl"))
    for pin in pins:
        candidates = [p for p in wheels if parse_wheel_filename(p.name)[:2] == (pin["name"], Version(pin["version"]))]
        if len(candidates) != 1:
            raise ValueError("Require exactly one wheel per pinned package: " + pin["name"])
        packages.append(metadata(candidates[0], pin))
    environment = default_environment()
    result = dict(status="declared_wheel_dependency_audit", lock=lock_ref, packages=packages,
        environment=environment, python=platform.python_version(), source=record(Path(__file__).resolve()),
        **evaluate(packages, environment), publication_ready=False,
        limitations=["Only declared wheel runtime dependencies with no optional extras on the recorded interpreter.",
            "Active URL or extra-bearing requirements are unresolved, never silently ignored.",
            "No installation, imports, wheel-tag/native-library compatibility, OS closure or vulnerability assessment.",
            "Missing or incorrect upstream dependency declarations cannot be detected here."])
    for ref in [lock_ref, *[p["wheel"] for p in packages]]:
        check(ref)
    with output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("lock", "wheelhouse", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.lock, args.wheelhouse, args.output)
    print(json.dumps(dict(packages=len(result["packages"]), edges=len(result["dependency_rows"]), failures=result["failures"])))
