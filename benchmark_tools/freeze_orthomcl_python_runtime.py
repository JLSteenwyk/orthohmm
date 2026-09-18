"""Inspect and inventory the dedicated Python environment used for BPO preparation."""

import argparse
import importlib.metadata
import json
from pathlib import Path
import site
import sys
import sysconfig

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.snapshot_runtime_trees import inventory, verify
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def require_packages(packages):
    if packages != {"biopython": "1.86", "numpy": "2.2.6"}:
        raise ValueError("Unexpected dedicated Python package inventory")


def minimal_roots(paths):
    roots = sorted(set(Path(p).resolve() for p in paths), key=lambda p: (len(p.parts), str(p)))
    selected = []
    for path in roots:
        if not any(path == root or path.is_relative_to(root) for root in selected):
            selected.append(path)
    return selected


def inspect():
    if (not sys.flags.isolated or not sys.dont_write_bytecode or site.ENABLE_USER_SITE
            or sys.prefix == sys.base_prefix or sys.version_info[:2] != (3, 10)
            or not sys.pycache_prefix or Path(sys.pycache_prefix).exists()):
        raise ValueError("Require isolated Python 3.10 venv, disabled bytecode writes and absent cache prefix")
    prefix = Path(sys.prefix).resolve()
    purelib = Path(sysconfig.get_path("purelib")).resolve()
    if not purelib.is_relative_to(prefix):
        raise ValueError("Package directory outside dedicated environment")
    packages = {d.metadata["Name"].lower(): d.version
                for d in importlib.metadata.distributions(path=[str(purelib)])}
    require_packages(packages)
    import benchmark_tools.prepare_qfo_corrected_bpo  # Exercise the production import chain.
    modules = {name: str(Path(mod.__file__).resolve()) for name, mod in list(sys.modules.items())
               if getattr(mod, "__file__", None) and not str(mod.__file__).startswith("<")}
    site_modules = {name: path for name, path in modules.items() if "/site-packages/" in path}
    if any(not Path(path).is_relative_to(purelib) or name.split(".")[0] not in {"Bio", "numpy"}
           for name, path in site_modules.items()):
        raise ValueError("Unexpected site module or editable startup hook")
    mapped = set()
    for line in Path("/proc/self/maps").read_text().splitlines():
        fields = line.split(maxsplit=5)
        if len(fields) == 6 and fields[5].startswith("/"):
            path = Path(fields[5])
            if not path.is_file():
                raise ValueError("Missing/deleted mapped runtime file: " + str(path))
            mapped.add(path.resolve())
    stdlib = Path(sysconfig.get_path("stdlib")).resolve()
    roots = minimal_roots([prefix, Path(sys.executable).resolve(), *mapped,
        *[p for p in stdlib.iterdir() if p.name not in {"site-packages", "__pycache__"}
          and p.suffix not in {".pyc", ".pyo"}]])
    return {"python_version": sys.version, "prefix": str(prefix), "base_prefix": sys.base_prefix,
            "executable": record(sys.executable), "packages": packages,
            "site_modules": site_modules, "mapped_files": [record(p) for p in sorted(mapped)],
            "imported_sources": [record(p) for p in sorted(set(modules.values()))],
            "roots": list(map(str, roots))}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    observed = inspect()
    runtime = inventory(observed["roots"])
    verification = verify(runtime)
    result = {"status": "dedicated_bpo_python_runtime_inventoried", "source": record(__file__),
              "inspection": observed, "runtime": runtime, "verification": verification,
              "execution_authorized": False, "publication_ready": False,
              "limitations": [
                  "Import-chain and mapped-file inventory plus full environment and base stdlib; not an OS-wide hermetic image.",
                  "Base standard library and interpreter are shared read-only dependencies and require before/after verification.",
                  "Actual imported helper sources are recorded separately; a frozen execution worktree remains required.",
                  "Wheel lock targets CPython 3.10 Linux x86_64; other platforms need separately verified artifacts."]}
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps({"status": result["status"], "runtime_records": len(runtime["records"])}))
