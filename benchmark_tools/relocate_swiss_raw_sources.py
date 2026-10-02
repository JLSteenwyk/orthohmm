"""Stage exact retained raw dependencies without altering historical manifests."""

import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import shutil

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


def frozen_json(path, digest):
    data = Path(path).read_bytes()
    if hashlib.sha256(data).hexdigest() != digest:
        raise ValueError("Frozen relocation manifest changed")
    return json.loads(data)


def validate_records(items):
    if not isinstance(items, list) or not items:
        raise ValueError("Require nonempty retained source records")
    seen = {}
    for item in items:
        if (not isinstance(item, dict) or set(item) != {"path", "bytes", "sha256"}
                or not isinstance(item["path"], str) or not Path(item["path"]).is_absolute()
                or type(item["bytes"]) is not int or item["bytes"] < 0
                or not isinstance(item["sha256"], str) or len(item["sha256"]) != 64
                or any(c not in "0123456789abcdef" for c in item["sha256"])):
            raise ValueError("Malformed retained source record")
        prior = seen.setdefault(item["path"], item)
        if prior != item:
            raise ValueError("Conflicting repeated source identity")


def restore_inputs(items, source, bindings=None, bindings_sha=None):
    """Return hash-checked relocated records, retaining original identities separately."""
    if bindings is None and bindings_sha is None:
        return items, None
    if bindings is None or bindings_sha is None:
        raise ValueError("Relocation requires both manifest and independent digest")
    validate_records(items)
    bindings = Path(bindings).resolve()
    document = frozen_json(bindings, bindings_sha)
    if (not isinstance(document, dict)
            or set(document) != {"schema_version", "source", "entries", "redistribution_authorized"}
            or type(document["schema_version"]) is not int or document["schema_version"] != 1
            or document["redistribution_authorized"] is not False):
        raise ValueError("Unexpected relocation schema or rights claim")
    validate_records([document["source"]])
    if any(document["source"][k] != source[k] for k in ("bytes", "sha256")):
        raise ValueError("Relocation belongs to a different frozen source manifest")
    entries = document["entries"]
    if not isinstance(entries, list) or len(entries) != len(items):
        raise ValueError("Relocation must preserve every source occurrence")
    root, restored = bindings.parent, []
    for original, entry in zip(items, entries):
        if not isinstance(entry, dict) or set(entry) != {"original", "artifact"}:
            raise ValueError("Malformed relocation entry")
        if entry["original"] != original or entry["artifact"] != "inputs/" + original["sha256"]:
            raise ValueError("Changed source identity, order or content address")
        relative = PurePosixPath(entry["artifact"])
        path = root.joinpath(*relative.parts)
        if path.resolve() != path or not path.is_file():
            raise ValueError("Relocated input must be a regular file without symlink parents")
        item = {**original, "path": str(path)}
        check(item)
        restored.append(item)
    binding_record = record(bindings)
    if binding_record["sha256"] != bindings_sha:
        raise ValueError("Relocation manifest changed during validation")
    return restored, dict(binding=binding_record, helper=record(__file__),
                         original_source=document["source"], entries=entries,
                         redistribution_authorized=False)


def stage(source, source_sha, key, output):
    """Copy pinned inputs once by content address; retain partial attempts on failure."""
    source, output = Path(source).resolve(), Path(output).absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    if output.resolve() != output:
        raise ValueError("Staging destination has symlink parents")
    if key not in ("checked_inputs", "records"):
        raise ValueError("Unexpected raw record collection")
    data = frozen_json(source, source_sha)
    if not isinstance(data, dict):
        raise ValueError("Require a frozen source manifest object")
    source_record = record(source)
    if source_record["sha256"] != source_sha:
        raise ValueError("Frozen source manifest changed during validation")
    items = data[key]
    validate_records(items)
    for item in items:
        check(item)
    output.mkdir(parents=True)
    (output / "inputs").mkdir()
    entries = []
    for item in items:
        artifact = "inputs/" + item["sha256"]
        target = output / artifact
        if not target.exists():
            with Path(item["path"]).open("rb") as src, target.open("xb") as dst:
                shutil.copyfileobj(src, dst, length=1024 * 1024)
        check({**item, "path": str(target)})
        entries.append(dict(original=item, artifact=artifact))
    for item in items:
        check(item)
    check(source_record)
    document = dict(schema_version=1, source=source_record, entries=entries,
                    redistribution_authorized=False)
    path = output / "bindings.json"
    with path.open("x") as stream:
        json.dump(document, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    restore_inputs(items, source_record, path, record(path)["sha256"])
    return record(path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--source-sha256", required=True)
    parser.add_argument("--records-key", choices=("checked_inputs", "records"), required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(stage(args.source, args.source_sha256, args.records_key, args.output), sort_keys=True))
