"""Invented relocation fixtures; no selected family scores or inference rerun."""

import csv
import io
import json
from pathlib import Path
import shutil
import statistics
import tarfile

import pytest

from benchmark_tools import readback_swiss_model_divergence as edge
from benchmark_tools import readback_swiss_model_divergence_strata as rational
from benchmark_tools import reproduce_swiss_model_divergence as replay


def index_fixture(tmp_path, monkeypatch):
    directory = tmp_path / "component"
    directory.mkdir()
    rows, pins = [], {}
    for name in sorted(set(replay.PINS) | set(replay.SUPPORT)):
        data = ("invented " + name).encode()
        (directory / name).write_bytes(data)
        ref = replay.identity(data)
        rows.append(dict(path=name, **ref))
        if name in replay.PINS:
            pins[name] = ("invented/" + name, ref["sha256"])
    monkeypatch.setattr(replay, "PINS", pins)
    index = dict(schema=replay.SCHEMA, files=rows, publication_ready=False,
                 native_inference_reproduced=False, raw_counts_readmitted=False,
                 new_bootstrap_draws=0, redistribution_clearance=False)
    replay.save(directory / replay.INDEX, index)
    return directory, index


def index_sha(directory):
    return replay.identity((directory / replay.INDEX).read_bytes())["sha256"]


def rewrite_index(directory, value):
    (directory / replay.INDEX).write_text(json.dumps(value))


def test_checked_flat_component(tmp_path, monkeypatch):
    directory, index = index_fixture(tmp_path, monkeypatch)
    assert replay.checked_component(directory, index_sha(directory)) == (directory, index)


@pytest.mark.parametrize("corruption", ["anchor", "payload", "extra", "duplicate", "symlink", "scope", "pin"])
def test_component_rejects_corruption_before_loading_copied_code(tmp_path, monkeypatch, corruption):
    directory, index = index_fixture(tmp_path, monkeypatch)
    anchor = index_sha(directory)
    if corruption == "anchor":
        anchor = "0" * 64
    elif corruption == "payload":
        (directory / "counts.json").write_text("changed")
    elif corruption == "extra":
        (directory / "unexpected").mkdir()
    elif corruption == "duplicate":
        index["files"][-1] = index["files"][0]
        rewrite_index(directory, index)
        anchor = index_sha(directory)
    elif corruption == "symlink":
        path = directory / "counts.json"
        data = path.read_bytes()
        path.unlink()
        target = tmp_path / "outside"
        target.write_bytes(data)
        path.symlink_to(target)
    elif corruption == "scope":
        index["native_inference_reproduced"] = True
        rewrite_index(directory, index)
        anchor = index_sha(directory)
    elif corruption == "pin":
        replay.PINS["counts.json"] = ("invented/counts.json", "0" * 64)
    with pytest.raises(ValueError):
        replay.checked_component(directory, anchor)


@pytest.mark.parametrize("kind", ["directory", "dangling_symlink"])
def test_refuses_occupied_output_before_git_access(tmp_path, kind):
    target = tmp_path / "occupied"
    if kind == "directory":
        target.mkdir()
    else:
        target.symlink_to(tmp_path / "missing")
    with pytest.raises(FileExistsError):
        replay.build(tmp_path, "unused", target)


def archive_fixture(tmp_path, corruption=None):
    archive = tmp_path / "evidence.tar.gz"
    refs = []
    names = [f"invented/{i:03d}.txt" for i in range(183)]
    if corruption == "unsafe":
        names[0] = "../escape"
    for name in names:
        refs.append(dict(path="/old/benchmarks/results/" + name,
                         **replay.identity(name.encode())))
    with tarfile.open(archive, "w:gz") as container:
        for i, name in enumerate(names):
            if corruption == "missing" and i == 182:
                continue
            member = tarfile.TarInfo(names[0] if corruption == "duplicate" and i == 182 else name)
            if corruption == "symlink" and i == 0:
                member.type, member.linkname = tarfile.SYMTYPE, "/outside"
                container.addfile(member)
                continue
            raw = b"changed" if corruption == "checksum" and i == 0 else name.encode()
            member.size = len(raw)
            container.addfile(member, io.BytesIO(raw))
    return archive, dict(checked_inputs=refs)


def test_all_payloads_checked_and_restored_to_new_namespace(tmp_path):
    archive, features = archive_fixture(tmp_path)
    output = tmp_path / "restored"
    assert replay.restore_payloads(archive, features, output) == 183
    assert len(list(output.rglob("*.txt"))) == 183
    assert (output / "invented/000.txt").read_bytes() == b"invented/000.txt"


@pytest.mark.parametrize("corruption", ["unsafe", "missing", "duplicate", "symlink", "checksum"])
def test_bad_evidence_is_not_extracted(tmp_path, corruption):
    archive, features = archive_fixture(tmp_path, corruption)
    output = tmp_path / "restored"
    with pytest.raises(ValueError):
        replay.restore_payloads(archive, features, output)
    assert not output.exists() and not (tmp_path / "escape").exists()


def test_historical_mapping_requires_one_exact_namespace():
    assert replay.archive_name("/old/benchmarks/results/family/tree") == "family/tree"
    for path in ("/old/other/results/file", "/benchmarks/results/benchmarks/results/file"):
        with pytest.raises(ValueError):
            replay.archive_name(path)


def numerical_fixture(tmp_path):
    directory = tmp_path / "component"
    directory.mkdir()
    shutil.copyfile(edge.__file__, directory / "edge_reader.py")
    shutil.copyfile(rational.__file__, directory / "rational_reader.py")
    evidence = tmp_path / "restored"
    run = evidence / "swiss_model_divergence_20261007_v1"
    memberships, descriptors, pair_rows = {}, {}, []
    for family, lengths in (("InventedA", (0.1, 0.2, 0.3)), ("InventedB", (1.0, 1.5, 2.0))):
        members = [family + str(i) for i in range(3)]
        memberships[family] = members
        path = run / family / "inference.treefile"
        path.parent.mkdir(parents=True)
        path.write_text("(" + ",".join(f"{g}:{v}" for g, v in zip(members, lengths)) + ");\n")
        descriptors[family], pairs = edge.edge_distances(path, members)
        pair_rows.extend(dict(family=family, gene_a=a, gene_b=b, distance=d) for (a, b), d in pairs.items())
    cutoff = statistics.median(v["median_pair_distance"] for v in descriptors.values())
    bins = dict(all=sorted(memberships), lower_or_equal_median=["InventedA"], higher_than_median=["InventedB"])
    feature_report = dict(memberships=memberships, strata=bins, median_family_distance=cutoff)
    (run / "report.json").write_text(json.dumps(feature_report))
    with (run / "pairs.tsv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, delimiter="\t", fieldnames=("family", "gene_a", "gene_b", "distance"))
        writer.writeheader()
        writer.writerows(pair_rows)
    family_rows = []
    points = {}
    for i, cell in enumerate(rational.CELLS):
        for j, family in enumerate(sorted(memberships)):
            counts = dict(TP=2 + i + j, FP=i + 1, FN=j + 2, TN=7)
            points[cell, family] = rational.point(counts)
            family_rows.append(dict(cell=cell, family=family, counts_without_prior=counts,
                                    **{k: float(v) for k, v in points[cell, family].items()}))
    rows, differences = [], []
    for stratum, families in bins.items():
        values = {c: rational.macro([points[c, f] for f in families]) for c in rational.CELLS}
        for cell in rational.CELLS:
            rows.append(dict(cell=cell, stratum=stratum, families=len(families), family_members=families,
                status="descriptive", prediction_semantics="native_pair" if cell == rational.CELLS[1] else "group_clique",
                **{m: float(v) for m, v in values[cell].items()}))
        for contrast, cell in (("R_at_P0_C0", rational.CELLS[1]), ("C_at_P0_R0", rational.CELLS[2])):
            differences.append(dict(contrast=contrast, candidate=cell, reference=rational.CELLS[0],
                stratum=stratum, families=len(families), family_members=families, status="descriptive",
                **{m + "_pp": float(100 * (values[cell][m] - values[rational.CELLS[0]][m])) for m in rational.METRICS}))
    for name, records, fields in (
        ("scores.tsv", rows, ("stratum", "cell", "families", "status", *rational.METRICS, "prediction_semantics")),
        ("differences.tsv", differences, ("stratum", "contrast", "families", "status", *(m + "_pp" for m in rational.METRICS)))):
        with (directory / name).open("w", newline="") as stream:
            writer = csv.DictWriter(stream, delimiter="\t", fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(records)
    lines = ["| " + " | ".join([r["stratum"], str(r["families"]), r["cell"],
        *[f"{100 * r[m]:.3f}" for m in rational.METRICS]]) + " |" for r in rows]
    lines += ["| " + " | ".join([r["stratum"], str(r["families"]), r["contrast"],
        *[f"{r[m + '_pp']:+.3f}" for m in rational.METRICS]]) + " |" for r in differences]
    (directory / "TABLE.md").write_text("\n".join(lines) + "\n")
    features = dict(features=descriptors, strata=bins, median_family_distance=cutoff)
    counts = dict(memberships=memberships, family_rows=family_rows)
    projection = dict(memberships=memberships, family_rows=family_rows, bins=bins,
                      median_family_distance=cutoff, rows=rows, differences=differences)
    return directory, evidence, features, counts, projection


def test_toy_relocated_arithmetic_uses_only_passed_inputs(tmp_path):
    result = replay.replay_values(*numerical_fixture(tmp_path))
    assert tuple(result[k] for k in ("families", "proteins", "pairs", "family_rows", "score_rows", "differences")) == (2, 6, 6, 6, 9, 6)
    assert result["bins"]["lower_or_equal_median"] == ["InventedA"]


@pytest.mark.parametrize("corruption", ["cutoff", "score", "pair", "table"])
def test_replayed_values_reject_changed_features_points_and_presentation(tmp_path, corruption):
    args = numerical_fixture(tmp_path)
    directory, evidence, features, _, projection = args
    if corruption == "cutoff":
        features["median_family_distance"] += 1
    elif corruption == "score":
        projection["rows"][1]["F1"] += 0.01
    elif corruption == "pair":
        path = evidence / "swiss_model_divergence_20261007_v1/pairs.tsv"
        path.write_text(path.read_text().replace("0.30000000000000004", "0.9"))
    else:
        (directory / "TABLE.md").write_text("not a table")
    with pytest.raises(ValueError):
        replay.replay_values(*args)


def test_created_attempt_retains_failure_without_retry_or_success_receipt(tmp_path, monkeypatch):
    import Bio
    import numpy
    monkeypatch.setattr(Bio, "__version__", "1.87")
    monkeypatch.setattr(numpy, "__version__", "2.2.6")
    directory = tmp_path / "component"
    directory.mkdir()
    for name, value in {
        "features.json": {"status": "features_verified"},
        "counts.json": {}, "counts_reader.json": {"family_rows_checked": 54},
        "projection_reader.json": {"status": "projection_verified"},
        "projection.json": {"cells": [{}, {"timing_eligible": False, "timing_admitted": False}]},
    }.items():
        (directory / name).write_text(json.dumps(value))
    monkeypatch.setattr(replay, "checked_component", lambda *args: (directory, {}))

    def fail(*args):
        raise ValueError("invented corrupted payload")

    monkeypatch.setattr(replay, "restore_payloads", fail)
    output = tmp_path / "attempt"
    with pytest.raises(ValueError, match="invented corrupted"):
        replay.replay(directory, "invented", output)
    assert not (output / "replay.json").exists()
    failure = json.loads((output / "failed_replay.json").read_text())
    assert failure["status"] == "portable_replay_failed" and failure["retry"] is False
    assert failure["publication_ready"] is False
