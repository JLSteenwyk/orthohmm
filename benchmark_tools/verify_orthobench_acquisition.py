"""Verify a separately acquired upstream OrthoBench checkout, without executing it."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess

REVISION = "872d6f30592ab5ff837224db16a514b3f2bb916a"
PINS = {
    "orthobench_factorial_prepared_20260916.json": "5c325f4d77865e0c7571fe4bb4df0be0959977a49e1d22169f192518700f9382",
    "orthobench_paired_uncertainty_20260916.json": "660ead29c5b6ac0b8278cd1e62cdcdb0a513db81dda317e634b805e661d70ba9",
}


def expected_inputs(results):
    documents = {}
    for name, digest in PINS.items():
        data = (results / name).read_bytes()
        if hashlib.sha256(data).hexdigest() != digest:
            raise ValueError("Historical input manifest differs")
        documents[name] = json.loads(data)
    fasta = documents["orthobench_factorial_prepared_20260916.json"]["fasta_inputs"]
    refs = documents["orthobench_paired_uncertainty_20260916.json"]["inputs"]
    if (len(fasta), len(refs["references"]), len(refs["uncertain"])) != (12, 70, 11):
        raise ValueError("Unexpected historical input inventory")
    expected = {}
    for prefix, records in (("BENCHMARKS/Input", fasta), ("BENCHMARKS/RefOGs", refs["references"]),
            ("BENCHMARKS/RefOGs/low_certainty_assignments", refs["uncertain"])):
        for row in records:
            path = prefix + "/" + Path(row["path"]).name
            if path in expected:
                raise ValueError("Duplicate expected input")
            expected[path] = {k: row[k] for k in ("bytes", "sha256")}
    return expected


def verify(repo, results):
    def git(*args):
        return subprocess.check_output(["git", *args], cwd=repo)
    if git("rev-parse", "HEAD").decode().strip() != REVISION:
        raise ValueError("Require the frozen upstream commit")
    expected = expected_inputs(results)
    paths = [*expected, "BENCHMARKS/benchmark.py", "README.md"]
    rows = []
    for name in paths:
        path = repo / name
        if path.is_symlink() or not path.is_file():
            raise ValueError("Missing or indirect upstream file")
        content = path.read_bytes()
        committed = git("show", REVISION + ":" + name)
        if content != committed:
            raise ValueError("Working file differs from frozen upstream blob")
        identity = dict(bytes=len(content), sha256=hashlib.sha256(content).hexdigest())
        if name in expected and identity != expected[name]:
            raise ValueError("Upstream file differs from scientific input manifest")
        rows.append(dict(relative_path=name, **identity,
            git_blob_sha1=hashlib.sha1(b"blob " + str(len(content)).encode() + b"\0" + content).hexdigest()))
    expected_paths = set(expected)
    actual_paths = {str(p.relative_to(repo)) for p in (repo / "BENCHMARKS/Input").glob("*")}
    actual_paths.update(str(p.relative_to(repo)) for p in (repo / "BENCHMARKS/RefOGs").rglob("*.txt"))
    if actual_paths != expected_paths:
        raise ValueError("Extra or missing benchmark inputs")
    return dict(status="orthobench_acquisition_verified", upstream_revision=REVISION,
        upstream_url="https://github.com/davidemms/Open_Orthobench.git", files=rows,
        historical_manifest_pins=PINS, source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        redistribution_cleared=False, scientific_results_reproduced=False,
        limitations=["Verifies local bytes against a fixed upstream commit and retained scientific input hashes.",
            "Raw inputs and upstream scorer remain separately acquired, not bundled or redistribution-cleared.",
            "No upstream code, inference or scoring was executed; this is not a complete publication release."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkout", type=Path, required=True)
    parser.add_argument("--results", type=Path, default=Path(__file__).resolve().parent / "results")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    result = verify(args.checkout.resolve(), args.results.resolve())
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
