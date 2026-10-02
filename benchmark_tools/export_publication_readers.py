"""Export the static analysis-module closure of the independent pipeline readers."""

import argparse
import ast
import hashlib
import json
from pathlib import Path
import subprocess

REVISION = "487d759"
ENTRYPOINT = "benchmark_tools.audit_publication_pipeline"


def identity(data):
    return dict(bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def dependencies(payload):
    names = set()
    for node in ast.walk(ast.parse(payload)):
        if isinstance(node, ast.Import):
            names.update(a.name for a in node.names if a.name.startswith("benchmark_tools."))
        elif isinstance(node, ast.ImportFrom):
            if node.level:
                raise ValueError("Relative imports require an explicit export rule")
            if node.module == "benchmark_tools":
                names.update("benchmark_tools." + a.name for a in node.names)
            elif node.module and node.module.startswith("benchmark_tools."):
                names.add(node.module)
    return sorted(names)


def closure(reader):
    pending, payloads = [ENTRYPOINT], {}
    while pending:
        module = pending.pop()
        if module in payloads:
            continue
        payload = reader(module.replace(".", "/") + ".py")
        payloads[module] = payload
        pending.extend(dependencies(payload))
    return {m.replace(".", "/") + ".py": p for m, p in sorted(payloads.items())}


def verify(directory):
    manifest = json.loads((directory / "manifest.json").read_text())
    expected = manifest["files"]
    observed = {str(p.relative_to(directory)) for p in directory.rglob("*") if p.is_file()}
    if observed != set(expected) | {"manifest.json"}:
        raise ValueError("Changed export inventory")
    for name, item in expected.items():
        relative = Path(name)
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError("Invalid export path")
        path = directory / relative
        if path.is_symlink() or identity(path.read_bytes()) != item:
            raise ValueError("Changed exported source")
    return manifest


def write_export(payloads, revision, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    for name in payloads:
        relative = Path(name)
        if relative.is_absolute() or ".." in relative.parts or relative.as_posix() != name:
            raise ValueError("Invalid export path")
    output.mkdir(parents=True)
    for name, payload in payloads.items():
        path = output / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(payload)
        path.chmod(0o644)
    manifest = dict(revision=revision, entrypoint=ENTRYPOINT,
        files={name: identity(payload) for name, payload in payloads.items()},
        publication_ready=False, scope="Static benchmark_tools import closure; reader entrypoint only",
        limitations=["External Python dependencies and scientific artifacts are not included",
            "Dynamic imports and arbitrary helper entrypoints are not guaranteed by static closure",
            "Manifest hashes detect changes, not authenticity; pin the manifest independently",
            "No native execution, independent biological accuracy or runtime relocation implied"])
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    (output / "manifest.json").chmod(0o644)
    return verify(output)


def export(repo, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    revision = subprocess.check_output(["git", "rev-parse", REVISION + "^{commit}"], cwd=repo, text=True).strip()
    payloads = closure(lambda name: subprocess.check_output(["git", "show", revision + ":" + name], cwd=repo))
    return write_export(payloads, revision, output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(export(args.repo.resolve(), args.output.absolute()), indent=2))
