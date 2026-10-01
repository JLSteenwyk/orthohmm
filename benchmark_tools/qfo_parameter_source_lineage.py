"""Explicit historical binding for one superseded, non-scientific exporter."""

import hashlib
from pathlib import Path
import subprocess

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

SOURCE = "benchmark_tools/export_qfo_corrected_comparison.py"
ARCHIVE = "benchmarks/work/publication_qfo_cpm_candidates_admission_v3/" + SOURCE
OLD_COMMIT = "5efb206b23a44f85386d8cb7e90f3e815a3d162a"
CURRENT_COMMIT = "01ac4f66b2ac4b3565e5617b7abfe789e7ffe9aa"
OLD_BYTES = 8652
OLD_SHA = "4239f8d263295b313c5cc27b866bf7c6494feac44a8fbfe18e160bd92ff69bfb"
CURRENT_BYTES = 9572
CURRENT_SHA = "8fd211aeaaac5765b6f5ffa01812d8e61ac4134b8c1ceb1bada4f497df7d7aa0"


def direct_record(path):
    path = Path(path)
    if not path.is_absolute() or path.resolve() != path or path.is_symlink():
        raise ValueError("Require direct historical-source path")
    return record(path)


def identity(root, name, size, sha):
    return dict(path=str(root / name), bytes=size, sha256=sha)


def verify_blob(root, revision, ref):
    blob = subprocess.check_output(["git", "-C", str(root), "show", revision + ":" + SOURCE])
    if len(blob) != ref["bytes"] or hashlib.sha256(blob).hexdigest() != ref["sha256"]:
        raise ValueError("Historical/current exporter differs from retained Git blob")


def check_records(records, root, historical_exporter_binding=False):
    """Check unchanged evidence, optionally binding the exact old exporter to Git."""
    root = Path(root)
    original = identity(root, SOURCE, OLD_BYTES, OLD_SHA)
    archived = identity(root, ARCHIVE, OLD_BYTES, OLD_SHA)
    current = identity(root, SOURCE, CURRENT_BYTES, CURRENT_SHA)
    lineage = None
    if historical_exporter_binding:
        if original not in records:
            raise ValueError("Historical exporter binding requested without original source evidence")
        if direct_record(root / SOURCE) != current or direct_record(root / ARCHIVE) != archived:
            raise ValueError("Known historical/current exporter identity changed")
        verify_blob(root, OLD_COMMIT, archived)
        verify_blob(root, CURRENT_COMMIT, current)
        lineage = dict(original_record=original, historical_copy=archived,
            historical_revision=OLD_COMMIT, current_record=current, current_revision=CURRENT_COMMIT,
            role="transitively captured comparator-report exporter; not parameter scoring or bootstrap",
            source=record(__file__), scientific_source_substitution=False)
    effective = {}
    for ref in records:
        bound = archived if lineage is not None and ref == original else ref
        prior = effective.get(bound["path"])
        if prior is not None and prior != bound:
            raise ValueError("Conflicting evidence identities during historical-source binding")
        effective[bound["path"]] = bound
    if lineage is not None:
        for ref in (current, lineage["source"]):
            prior = effective.get(ref["path"])
            if prior is not None and prior != ref:
                raise ValueError("Unknown exporter/source identity cannot use historical binding")
            effective[ref["path"]] = ref
    for ref in effective.values():
        check(ref)
    return list(effective.values()), lineage
