"""Copy only admitted, input-identical raw gene-tree caches into a fresh arm."""

import json
from pathlib import Path
import re
import shutil

from benchmark_tools.run_qfo_order_replay import record


def check(item):
    path = Path(item["path"])
    if path.is_symlink() or record(path) != item:
        raise ValueError("Changed or linked cache artifact: " + str(path))


def cache_paths(family):
    if re.fullmatch(r"Family[0-9]{7,}", family) is None:
        raise ValueError("Unsafe family identifier")
    return (f"checkpoints/{family}.json", f"gene_trees/{family}.raw.nwk",
            f"candidate_fastas/{family}.faa", f"alignments/{family}.faa")


def select_cache(source, families, admitted_outputs, expected_hash, tool_records):
    """expected_hash must recompute the native tree-input hash on target sequences.

    Caller must first admit the completed native run and scientific readback.
    This function checks artifact eligibility, not that upstream admission.
    """
    source = Path(source).resolve()
    for item in tool_records:
        check(item)
    if not tool_records:
        raise ValueError("External tool identities are required")
    admitted = {r["path"]: r for r in admitted_outputs}
    if len(admitted) != len(admitted_outputs):
        raise ValueError("Duplicate admission paths")
    expected_paths = {Path(r["path"]) for r in admitted_outputs
                      if Path(r["path"]).parent == source / "checkpoints"}
    if set((source / "checkpoints").glob("*.json")) != expected_paths:
        raise ValueError("Checkpoint inventory differs from admission")
    files, included, excluded = [], [], []
    for path in sorted(expected_paths):
        check(admitted[str(path)])
        checkpoint = json.loads(path.read_text())
        family = checkpoint["family_id"]
        relatives = cache_paths(family)
        genes = checkpoint["genes"]
        if (path.name != family + ".json" or checkpoint["status"] != "complete"
                or checkpoint["schema_version"] != 2 or not genes
                or not all(isinstance(g, str) and g for g in genes) or len(set(genes)) != len(genes)):
            raise ValueError("Invalid admitted checkpoint")
        target = families.get(family)
        if target is None or tuple(sorted(genes)) != tuple(sorted(target)):
            excluded.append(dict(family=family, reason="changed_or_absent_membership"))
            continue
        if checkpoint["tree_input_sha256"] != expected_hash(family, tuple(sorted(target))):
            excluded.append(dict(family=family, reason="changed_tree_inputs"))
            continue
        for relative in relatives:
            key = str(source / relative)
            if key not in admitted:
                raise ValueError("Unadmitted cache artifact: " + key)
            item = admitted[key]
            check(item)
            if relative.endswith(".raw.nwk") and item["sha256"] != checkpoint["raw_tree_sha256"]:
                raise ValueError("Raw tree differs from checkpoint")
            files.append(dict(relative=relative, source=item))
        included.append(family)
    for item in tool_records:
        check(item)
    return dict(source_directory=str(source), included_families=included,
                excluded_families=excluded, files=files, tool_records=tool_records,
                species_tree_cache_copied=False, reconciliation_outputs_copied=False)


def copy_cache(manifest, destination):
    """Use independent copies, never links to mutable previous-run artifacts."""
    destination = Path(destination)
    if destination.exists():
        raise FileExistsError(destination)
    allowed = {name for family in manifest["included_families"] for name in cache_paths(family)}
    names = [r["relative"] for r in manifest["files"]]
    if len(names) != len(set(names)) or set(names) != allowed:
        raise ValueError("Cache file inventory is incomplete or unsafe")
    for row in manifest["files"]:
        if row["source"]["path"] != str(Path(manifest["source_directory"]) / row["relative"]):
            raise ValueError("Cache source escapes its admitted directory")
        check(row["source"])
    for item in manifest["tool_records"]:
        check(item)
    destination.mkdir(parents=True)
    copied = []
    for row in manifest["files"]:
        path = destination / row["relative"]
        path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(row["source"]["path"], path)
        result = record(path)
        if any(result[k] != row["source"][k] for k in ("bytes", "sha256")):
            raise ValueError("Copied cache differs")
        check(row["source"])
        copied.append(result)
    return copied
