"""Inspect actual native interpreter lookup and exercised import identities."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.run_simulation_methods import read_frozen, execution_environment
from benchmark_tools.isolated_numba_cache import fresh_cache
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.probe_dgx_step_separation import save

PROBE = r'''
import importlib, json, os, pathlib, sys
requested = json.loads(sys.argv[1])
for name in requested:
    importlib.import_module(name)
def label(value):
    value = value if isinstance(value, type) or hasattr(value, '__qualname__') else type(value)
    return value.__module__ + '.' + value.__qualname__
paths = []
for entry in sys.path:
    path = pathlib.Path(entry or os.getcwd()).resolve()
    names = sorted(p.name for p in path.iterdir() if p.name not in {'__pycache__', '.git'}) if path.is_dir() else None
    paths.append(dict(entry=entry, resolved=str(path), is_directory=path.is_dir(), exists=path.exists(), names=names))
modules = {name: str(pathlib.Path(module.__file__).resolve()) for name, module in list(sys.modules.items())
           if getattr(module, '__file__', None) and not module.__file__.startswith('<')}
editable = {name: {key: getattr(module, key) for key in ('MAPPING', 'NAMESPACES') if hasattr(module, key)}
            for name, module in list(sys.modules.items()) if name.startswith('__editable__')}
mapped = set()
ipc_mappings = []
for line in pathlib.Path('/proc/self/maps').read_text().splitlines():
    fields = line.split(maxsplit=5)
    if len(fields) == 6 and fields[5].startswith('/'):
        path = pathlib.Path(fields[5])
        if not path.is_file():
            if fields[5].startswith(('/dev/shm/sem.', '/dev/shm/pym-')) and 'x' not in fields[1]:
                ipc_mappings.append(dict(path=fields[5], permissions=fields[1]))
                continue
            raise ValueError('Missing/deleted mapped runtime file: ' + str(path))
        mapped.add(str(path.resolve()))
print(json.dumps(dict(python=sys.version, executable=str(pathlib.Path(sys.executable).resolve()),
    cwd=os.getcwd(), requested=requested, paths=paths, modules=modules, mapped_files=sorted(mapped),
    meta_path=[label(x) for x in sys.meta_path], path_hooks=[label(x) for x in sys.path_hooks],
    editable=editable, ipc_mappings=ipc_mappings,
    dont_write_bytecode=sys.dont_write_bytecode, pycache_prefix=sys.pycache_prefix)))
'''


def scientific_origins(report, package, expected_root):
    root = Path(expected_root)
    selected = {name: path for name, path in report["modules"].items()
                if name == package or name.startswith(package + ".")}
    if not selected or any(not Path(path).is_relative_to(root) for path in selected.values()):
        raise ValueError("Scientific package resolved outside frozen source: " + package)
    return dict(package=package, root=str(root), checked_modules=len(selected))


def lookup_signature(report):
    return {key: report[key] for key in ("python", "executable", "cwd", "requested", "paths",
        "modules", "mapped_files", "meta_path", "path_hooks", "editable", "dont_write_bytecode", "files")}


def compare_lookup(expected, observed):
    if (expected["coverage"] != dict(missing=[], changed=[], all_covered=True)
            or observed["coverage"] != expected["coverage"]
            or not expected["dont_write_bytecode"] or not observed["dont_write_bytecode"]
            or expected["scientific_origin"] != observed["scientific_origin"]
            or lookup_signature(expected) != lookup_signature(observed)):
        raise ValueError("Native Python lookup or imported file identity differs")
    return dict(status="native_python_lookup_matches", scientific_execution_authorized=False,
                limitations=["Declared import-chain and startup lookup match, not continuous import enforcement."])


def covered_files(files, records):
    index = {}
    for row in records:
        if row["kind"] == "file":
            index[row["path"]] = (row["bytes"], row["sha256"])
        elif row["kind"] == "symlink" and "target_sha256" in row:
            index[row["resolved"]] = (row["target_bytes"], row["target_sha256"])
    missing, changed = [], []
    for row in files:
        if row["path"] not in index:
            missing.append(row["path"])
        elif index[row["path"]] != (row["bytes"], row["sha256"]):
            changed.append(row["path"])
    return dict(missing=missing, changed=changed, all_covered=not missing and not changed)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("baseline", "binding", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--baseline-sha256", required=True)
    parser.add_argument("--binding-sha256", required=True)
    args = parser.parse_args()
    baseline = read_frozen(args.baseline, args.baseline_sha256)
    binding = read_frozen(args.binding, args.binding_sha256)
    args.output.mkdir(parents=True, exist_ok=False)
    records = [r for path, sha in binding["runtime_specs"]
               for r in read_frozen(Path(path), sha)["records"]]
    env, _ = execution_environment(baseline)
    for key in ("LD_PRELOAD", "LD_LIBRARY_PATH", "LD_AUDIT"):
        env.pop(key, None)
    env.update(PYTHONDONTWRITEBYTECODE="1", PYTHONPYCACHEPREFIX=str(args.output / "python_cache"))
    interpreters = dict(orthohmm=baseline["tool_entrypoints"]["orthohmm_python"]["absolute_path"],
        orthofinder=str(Path(baseline["tool_entrypoints"]["orthofinder"]["absolute_path"]).parent / "python"))
    requested = dict(orthohmm=["orthohmm.orthohmm", "orthohmm.search.engine", "igraph", "leidenalg", "dendropy"],
                     orthofinder=["orthofinder.run.main"])
    of_roots = {str(Path(r["absolute_path"]).parent) for r in baseline["orthofinder_distribution"]
                if r["absolute_path"].endswith("/orthofinder/__init__.py")}
    if len(of_roots) != 1:
        raise ValueError("Ambiguous frozen OrthoFinder package root")
    package_roots = dict(orthohmm=str(Path(baseline["core_root"]) / "orthohmm"), orthofinder=of_roots.pop())
    reports = {}
    for name, executable in interpreters.items():
        with fresh_cache(args.output / (name + "_numba_cache")):
            child_env = dict(env, NUMBA_CACHE_DIR=os.environ["NUMBA_CACHE_DIR"],
                             NUMBA_CACHE_LOCATOR_CLASSES=os.environ["NUMBA_CACHE_LOCATOR_CLASSES"])
            completed = subprocess.run([executable, "-B", "-c", PROBE, json.dumps(requested[name])],
                cwd=baseline["core_root"], env=child_env, capture_output=True, text=True, timeout=90)
        save(args.output / (name + "_process.json"), dict(exit_code=completed.returncode,
                                                        stdout=completed.stdout, stderr=completed.stderr))
        if completed.returncode:
            raise RuntimeError("Lookup probe failed: " + name)
        report = json.loads(completed.stdout)
        paths = sorted(set(report["modules"].values()) | set(report["mapped_files"]) | {report["executable"]})
        report["files"] = [record(path) for path in paths]
        report["coverage"] = covered_files(report["files"], records)
        report["scientific_origin"] = scientific_origins(report, name, package_roots[name])
        report["lookup_signature"] = lookup_signature(report)
        save(args.output / (name + ".json"), report)
        reports[name] = dict(source=record(args.output / (name + ".json")),
                             modules=len(report["modules"]), coverage=report["coverage"])
    save(args.output / "result.json", dict(status="native_lookup_inspected", reports=reports,
        baseline=record(args.baseline), binding=record(args.binding), source=record(__file__),
        scientific_execution_authorized=False, limitations=[
            "Exercises declared imports only, not every native workload branch.",
            "Directory names are inspected nonrecursively; unrelated project contents are not inventoried.",
            "Checks observed file hashes against pinned snapshots; not a full runtime-tree revalidation.",
            "Lookup identities are observed, not a sandbox against mid-run mutation or new imports."]))


if __name__ == "__main__":
    main()
