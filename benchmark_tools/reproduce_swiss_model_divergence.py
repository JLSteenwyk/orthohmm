"""Replay retained distance/table arithmetic from an externally anchored component."""

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path, PurePosixPath
import statistics
import subprocess
import tarfile

PINS = {
    "evidence.tar.gz": ("benchmark_tools/results/swiss_model_divergence_evidence_23932_v1.tar.gz", "eeed6e6786f5720ee8c19fb6123fad64c38fdb25c33afe57f5e61bcae4a6113d"),
    "features.json": ("benchmark_tools/results/swiss_model_divergence_readback_23932_v1.json", "45cb12b3d201daf56e3fee5c5833153a69dd1732a2fee1c28d59261130eeb312"),
    "counts.json": ("benchmark_tools/results/native_qfo_three_cell_strata_20261007_v1/report.json", "55ed72ce771c756424877cfa412da329ff104efb7dca4b35113b0d80703a59e5"),
    "counts_reader.json": ("benchmark_tools/results/native_qfo_three_cell_strata_readback_20261007_v2.json", "6f3eef63c0e96e1d85805fcfb70ed9b66fbb1db76a08ca2a105483422a0f6969"),
    "projection.json": ("benchmark_tools/results/swiss_model_divergence_strata_20261007_v1/report.json", "ddfa27757c4c256df9c881079ac352c37d7e17646f75fbd417b624e7980be78f"),
    "projection_reader.json": ("benchmark_tools/results/swiss_model_divergence_strata_readback_20261007_v1.json", "186c0151137fd825513ef36c27df595b1a0f1cd4b2b4bf68de4b9665ad328d0e"),
    "scores.tsv": ("benchmark_tools/results/swiss_model_divergence_strata_20261007_v1/scores.tsv", "12d12819107fd9fd5436eb458110eea32f9870bccda23694811abc24c8b6a1c3"),
    "differences.tsv": ("benchmark_tools/results/swiss_model_divergence_strata_20261007_v1/differences.tsv", "3a120223d8f02d3eb2b1802c98829c7c4708543eb62ffac60d7ad76d15ebfa34"),
    "TABLE.md": ("benchmark_tools/results/swiss_model_divergence_strata_20261007_v1/TABLE.md", "46be38ef7ddd144e3b6f499ef1e3fae908281932012c0fb9d40cd80a2e37d739"),
    "protocol.md": ("benchmark_tools/results/SWISS_MODEL_DIVERGENCE_PROTOCOL_20261007.md", "a596371c443471585b9a277cbdd1466bb3f9870200f3187d20286461d3699bfd"),
    "edge_reader.py": ("benchmark_tools/readback_swiss_model_divergence.py", "eb453e0274f6730ec3c8765e3a99260f1eeab46e9f873321486c9c4f51b1f573"),
    "rational_reader.py": ("benchmark_tools/readback_swiss_model_divergence_strata.py", "c4d3a0e6fcaea31eda2374ae48b3adbd5009a898d5667430d785b85cb57f3bbd"),
}
SUPPORT = {"replay.py": "benchmark_tools/reproduce_swiss_model_divergence.py",
           "README.md": "benchmark_tools/SWISS_MODEL_DIVERGENCE_REPLAY.md",
           "requirements.txt": "benchmark_tools/swiss_model_divergence_replay_requirements.txt"}
INDEX = "REPLAY_INDEX.json"
SCHEMA = "swiss_model_divergence_portable_v1"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def identity(content):
    return dict(bytes=len(content), sha256=hashlib.sha256(content).hexdigest())


def fresh(path):
    path = Path(path).absolute()
    if path.exists() or path.is_symlink():
        raise FileExistsError(path)
    return path


def save(path, value):
    with Path(path).open("x", encoding="ascii") as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")


def build(repo, revision, output):
    repo, output = Path(repo).resolve(strict=True), fresh(output)
    commit = subprocess.check_output(["git", "-C", str(repo), "rev-parse", "--verify",
                                      revision + "^{commit}"], text=True).strip()
    contents, rows = {}, []
    for name, path in {**{n: p for n, (p, _) in PINS.items()}, **SUPPORT}.items():
        raw = subprocess.check_output(["git", "-C", str(repo), "show", commit + ":" + path])
        if name in PINS:
            require(identity(raw)["sha256"] == PINS[name][1], "Changed frozen input: " + name)
        contents[name] = raw
        rows.append(dict(path=name, git_path=path, git_revision=commit, **identity(raw)))
    require(contents["replay.py"] == Path(__file__).read_bytes(), "Build source not precommitted")
    output.mkdir(parents=True, exist_ok=False)
    for name, raw in contents.items():
        with (output / name).open("xb") as stream:
            stream.write(raw)
    index = dict(schema=SCHEMA, files=rows, publication_ready=False,
                 native_inference_reproduced=False, raw_counts_readmitted=False,
                 new_bootstrap_draws=0, redistribution_clearance=False)
    save(output / INDEX, index)
    return dict(status="portable_component_built", source_commit=commit, files=len(rows),
                payload_bytes=sum(r["bytes"] for r in rows),
                index=identity((output / INDEX).read_bytes()), publication_ready=False)


def checked_component(directory, index_sha):
    directory = Path(directory).resolve(strict=True)
    index_path = directory / INDEX
    require(not index_path.is_symlink() and identity(index_path.read_bytes())["sha256"] == index_sha,
            "Index differs from external anchor")
    index = json.loads(index_path.read_text())
    require(index["schema"] == SCHEMA and all(index[k] is False for k in (
        "publication_ready", "native_inference_reproduced", "raw_counts_readmitted",
        "redistribution_clearance")) and type(index["new_bootstrap_draws"]) is int
        and index["new_bootstrap_draws"] == 0, "Inflated component scope")
    expected = set(PINS) | set(SUPPORT)
    require(len(index["files"]) == len(expected)
            and {r["path"] for r in index["files"]} == expected, "Incorrect component inventory")
    require({p.name for p in directory.iterdir()} == expected | {INDEX}, "Extra component entries")
    for row in index["files"]:
        path = directory / row["path"]
        require(not path.is_symlink() and path.is_file(), "Nonregular component file")
        require(identity(path.read_bytes()) == {k: row[k] for k in ("bytes", "sha256")},
                "Changed component file: " + row["path"])
        if row["path"] in PINS:
            require(row["sha256"] == PINS[row["path"]][1], "Changed scientific pin")
    return directory, index


def archive_name(path):
    parts = PurePosixPath(path).parts
    matches = [i for i in range(len(parts) - 1) if parts[i:i + 2] == ("benchmarks", "results")]
    require(len(matches) == 1, "Ambiguous historical archive mapping")
    return PurePosixPath(*parts[matches[0] + 2:]).as_posix()


def restore_payloads(path, features, output):
    expected = {archive_name(r["path"]): r for r in features["checked_inputs"]
                if ("benchmarks", "results") in list(zip(PurePosixPath(r["path"]).parts,
                                                        PurePosixPath(r["path"]).parts[1:]))}
    require(len(expected) == 183, "Incomplete retained payload inventory")
    with tarfile.open(path, "r:gz") as archive:
        members = archive.getmembers()
        require(len(members) == len(expected) and {m.name for m in members} == set(expected)
                and all(m.isfile() for m in members), "Incorrect archive inventory or types")
        for member in members:
            relative = PurePosixPath(member.name)
            require(not relative.is_absolute() and ".." not in relative.parts
                    and str(relative) == member.name and "\\" not in member.name, "Unsafe member path")
            raw = archive.extractfile(member).read()
            require(identity(raw) == {k: expected[member.name][k] for k in ("bytes", "sha256")},
                    "Changed retained payload")
        output.mkdir(parents=True, exist_ok=False)
        for member in members:
            target = output / member.name
            target.parent.mkdir(parents=True, exist_ok=True)
            with target.open("xb") as stream:
                stream.write(archive.extractfile(member).read())
    return len(members)


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    exec(compile(path.read_bytes(), str(path), "exec"), module.__dict__)
    return module


def replay_values(directory, evidence, features, counts, projection):
    edge = load_module(directory / "edge_reader.py", "copied_edge_reader")
    rational = load_module(directory / "rational_reader.py", "copied_rational_reader")
    run = evidence / "swiss_model_divergence_20261007_v1"
    report = json.loads((run / "report.json").read_text())
    require(report["memberships"] == counts["memberships"] == projection["memberships"],
            "Member universes differ")
    computed, pairs = {}, {}
    for family, members in sorted(report["memberships"].items()):
        require(family.isalnum() or family.replace("_", "").isalnum(), "Unsafe family name")
        descriptor, family_pairs = edge.edge_distances(run / family / "inference.treefile", members)
        computed[family] = descriptor
        recorded = features["features"][family]
        require(set(recorded) == set(descriptor), "Changed descriptor fields")
        for key, value in descriptor.items():
            if key in ("members", "pairs", "unit"):
                require(recorded[key] == value, "Descriptor membership/unit differs")
            else:
                edge.close(recorded[key], value)
        pairs.update({(family, a, b): distance for (a, b), distance in family_pairs.items()})
    pair_count = edge.compare_pairs(run / "pairs.tsv", pairs)
    cutoff = statistics.median(v["median_pair_distance"] for v in computed.values())
    edge.close(report["median_family_distance"], cutoff)
    edge.close(features["median_family_distance"], cutoff)
    bins = dict(all=sorted(computed), lower_or_equal_median=sorted(f for f in computed
                if computed[f]["median_pair_distance"] <= cutoff), higher_than_median=sorted(f for f in computed
                if computed[f]["median_pair_distance"] > cutoff))
    require(bins == report["strata"] == features["strata"] == projection["bins"], "Changed cutoff/tie bins")
    rows, differences = rational.numerical_readback(projection, counts, report)
    rational.verify_tsv(directory / "scores.tsv", rows,
                        ("stratum", "cell", "families", "status", *rational.METRICS, "prediction_semantics"),
                        rational.METRICS)
    rational.verify_tsv(directory / "differences.tsv", differences,
                        ("stratum", "contrast", "families", "status", *(m + "_pp" for m in rational.METRICS)),
                        tuple(m + "_pp" for m in rational.METRICS))
    table = (directory / "TABLE.md").read_text()
    for row in rows:
        values = [f"{100 * float(row[m]):.3f}" for m in rational.METRICS]
        line = "| " + " | ".join([row["stratum"], str(row["families"]), row["cell"], *values]) + " |"
        require(table.count(line) == 1, "Incorrect human score row")
    for row in differences:
        values = [f"{float(row[m + '_pp']):+.3f}" for m in rational.METRICS]
        line = "| " + " | ".join([row["stratum"], str(row["families"]), row["contrast"], *values]) + " |"
        require(table.count(line) == 1, "Incorrect human contrast row")
    return dict(families=len(computed), proteins=sum(len(v["members"]) for v in computed.values()),
                pairs=pair_count, family_rows=len(projection["family_rows"]), score_rows=len(rows),
                differences=len(differences), median_family_distance=cutoff, bins=bins,
                family_descriptors=computed)


def replay(directory, index_sha, output):
    output = fresh(output)
    directory, index = checked_component(directory, index_sha)
    import Bio
    import numpy
    require(Bio.__version__ == "1.87", "Require the pinned Biopython 1.87 parser")
    require(numpy.__version__ == "2.2.6", "Require the pinned NumPy 2.2.6 dependency")
    features, counts, projection, counts_reader, projection_reader = [json.loads((directory / n).read_text())
        for n in ("features.json", "counts.json", "projection.json", "counts_reader.json", "projection_reader.json")]
    require(features["status"] == "features_verified" and counts_reader["family_rows_checked"] == 54
            and projection_reader["status"] == "projection_verified"
            and projection["cells"][1]["timing_eligible"] is False
            and projection["cells"][1]["timing_admitted"] is False, "Invalid inherited scope")
    output.mkdir(parents=True, exist_ok=False)
    try:
        restored = restore_payloads(directory / "evidence.tar.gz", features, output / "evidence")
        result = replay_values(directory, output / "evidence", features, counts, projection)
        require(tuple(result[k] for k in ("families", "proteins", "pairs", "family_rows", "score_rows", "differences"))
                == (18, 563, 10765, 54, 9, 6), "Incomplete selected population")
    except Exception as error:
        save(output / "failed_replay.json", dict(status="portable_replay_failed", retry=False,
             error_type=type(error).__name__, error=str(error), publication_ready=False))
        raise
    receipt = dict(status="retained_model_distance_arithmetic_replayed", index_sha256=index_sha,
                   component_files=len(index["files"]), restored_payloads=restored,
                   biopython_version=Bio.__version__, numpy_version=numpy.__version__, **result,
                   new_bootstrap_draws=0, native_inference_reproduced=False, raw_counts_readmitted=False,
                   scientific_timings_admitted=False, publication_ready=False, independent_confirmation=False,
                   limitations=["Reuses frozen edge/rational functions and Biopython parsing; not a second independent implementation.",
                                "Historical raw-count admission, inference, model fit, runtime and scheduler are not reexecuted.",
                                "No new uncertainty, default, isolated timing, generalization or public redistribution claim."])
    save(output / "replay.json", receipt)
    return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    builder = commands.add_parser("build")
    builder.add_argument("--repo", type=Path, required=True)
    builder.add_argument("--revision", required=True)
    builder.add_argument("--output", type=Path, required=True)
    runner = commands.add_parser("replay")
    runner.add_argument("directory", type=Path)
    runner.add_argument("--index-sha256", required=True)
    runner.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = build(args.repo, args.revision, args.output) if args.command == "build" else replay(
        args.directory, args.index_sha256, args.output)
    print(json.dumps(result, indent=2, sort_keys=True))
