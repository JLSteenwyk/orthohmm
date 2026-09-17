"""Create fresh nested dataset directories from the frozen transfer bundle."""

import argparse
import hashlib
import json
from pathlib import Path
import shutil

TRANSFER_SHA = "fcff891fc3d787df52585eb3788a3dbe6d400e8699c4f0f80ca2b66a4e6b0eb0"


def record(path):
    return {"path": str(path.resolve()), "bytes": path.stat().st_size,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def materialize(bundle, output, expected_sha=TRANSFER_SHA):
    bundle, output = Path(bundle).resolve(), Path(output).absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    manifest_path = bundle / "manifest.json"
    manifest_record = record(manifest_path)
    if manifest_record["sha256"] != expected_sha:
        raise ValueError("Wrong frozen transfer manifest")
    manifest = json.loads(manifest_path.read_text())
    files = {row["path"]: row for row in manifest["inputs"]}
    if len(files) != 12 or len(manifest["inputs"]) != 12:
        raise ValueError("Expected twelve unique proteomes")
    names = []
    for relative, row in files.items():
        path = bundle / relative
        if Path(relative).parts != ("inputs", Path(relative).name) or path.is_symlink():
            raise ValueError("Invalid transfer input path")
        actual = record(path)
        if any(actual[key] != row[key] for key in ("sha256", "bytes")):
            raise ValueError("Transfer input changed")
        names.append(path.name)
    if sorted(p.name for p in (bundle / "inputs").iterdir()) != sorted(names):
        raise ValueError("Unexpected transfer input inventory")
    if [d["proteomes"] for d in manifest["datasets"]] != [4, 8, 12]:
        raise ValueError("Changed nested dataset sizes")
    for dataset in manifest["datasets"]:
        if dataset["inputs"] != list(files)[:dataset["proteomes"]]:
            raise ValueError("Changed nested dataset membership")
    output.mkdir(parents=True)
    datasets = []
    for dataset in manifest["datasets"]:
        target = output / str(dataset["proteomes"])
        target.mkdir()
        copied = []
        for relative in dataset["inputs"]:
            destination = target / Path(relative).name
            shutil.copy2(bundle / relative, destination)
            actual = record(destination)
            if any(actual[key] != files[relative][key] for key in ("sha256", "bytes")):
                raise ValueError("Dataset copy differs from frozen input")
            copied.append(actual)
        datasets.append({**dataset, "input_directory": str(target.resolve()), "inputs": copied})
    if record(manifest_path) != manifest_record:
        raise ValueError("Transfer manifest changed while copying")
    report = {"status": "nested_inputs_materialized_unrun", "source": record(Path(__file__)),
              "transfer_manifest": manifest_record, "datasets": datasets,
              "planned_runs": manifest["planned_runs"], "inference_started": False}
    (output / "manifest.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    materialize(args.bundle, args.output)
