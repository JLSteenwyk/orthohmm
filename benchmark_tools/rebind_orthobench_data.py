"""Rebind frozen OrthoBench input identities to a separately acquired local tree."""

import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath


ROLES = {"fasta": "BENCHMARKS/Input", "references": "BENCHMARKS/RefOGs",
         "uncertain": "BENCHMARKS/RefOGs/low_certainty_assignments"}


def record(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return dict(path=str(path.absolute()), bytes=path.stat().st_size, sha256=digest.hexdigest())


def rebind(manifest, sha256, acquisition, output):
    manifest, acquisition, output = Path(manifest), Path(acquisition), Path(output)
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    original = record(manifest)
    if original["sha256"] != sha256:
        raise ValueError("Frozen data manifest checksum differs")
    data = json.loads(manifest.read_text())
    if (set(data) != {"dataset", "genes", "upstream_revision", *ROLES}
            or data["dataset"] != "orthobench" or data["genes"] != 251378
            or [len(data[r]) for r in ROLES] != [12, 70, 11]):
        raise ValueError("Require the full OrthoBench manifest schema and dimensions")
    acquisition = acquisition.resolve(strict=True)
    result, checked = dict(data), []
    original_roots = set()
    for role, directory in ROLES.items():
        rows, names = [], set()
        for old in data[role]:
            if set(old) != {"path", "bytes", "sha256"}:
                raise ValueError("Unexpected input identity fields")
            path = PurePosixPath(old["path"])
            suffix = PurePosixPath(directory)
            if (not path.is_absolute() or ".." in path.parts or "\\" in old["path"]
                    or path.parts[-len(suffix.parts)-1:-1] != suffix.parts
                    or path.name in names):
                raise ValueError("Invalid original role path or duplicate basename")
            original_roots.add(str(path.parents[len(suffix.parts)]))
            names.add(path.name)
            local = acquisition / directory / path.name
            if local.is_symlink() or not local.resolve().is_relative_to(acquisition):
                raise ValueError("Input symlink or path escapes acquisition root")
            item = record(local)
            if any(item[k] != old[k] for k in ("bytes", "sha256")):
                raise ValueError("Acquired input differs from frozen bytes: " + path.name)
            rows.append(item)
            checked.append(item)
        result[role] = rows
    if len(original_roots) != 1:
        raise ValueError("Original inputs do not share one acquisition root")
    if record(manifest) != original or any(record(Path(r["path"])) != r for r in checked):
        raise ValueError("Input changed during rebinding")
    output.mkdir(parents=True, exist_ok=False)
    target = output / "data.json"
    with target.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    receipt = dict(status="orthobench_data_paths_rebound", original_manifest=original,
                   acquisition_root=str(acquisition), manifest=record(target), files=len(checked),
                   source=record(Path(__file__).resolve()), inference_rerun=False,
                   publication_ready=False, limitations=[
                       "Verifies frozen file bytes and order, not acquisition history or data rights.",
                       "Upstream revision is retained metadata, not a fresh Git repository audit.",
                       "Caller supplies the manifest digest independently and retains the new digest."])
    with (output / "rebind.json").open("x") as handle:
        json.dump(receipt, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("manifest", "acquisition", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--sha256", required=True)
    args = parser.parse_args()
    print(json.dumps(rebind(args.manifest, args.sha256, args.acquisition, args.output)))
