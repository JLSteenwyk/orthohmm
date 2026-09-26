"""Same-host traced relocation of external phylogeny tools; no redistribution grant."""

import argparse
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess

from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.verify_frozen_phylogeny_install import fixture, validate

PRIOR = "benchmark_tools/results/publication_frozen_phylogeny_smoke_20260926.json"
PRIOR_SHA = "8d258122281bed99d26e4c66bdc03aa3b8bfa3367db3de1bdb5daf86c815e4e0"


def copy_regular(source, destination):
    if source.is_symlink() or not source.is_file():
        raise ValueError("Require a regular, nonsymlink tool file")
    original = record(source)
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    copied = record(destination)
    mode = source.stat().st_mode & 0o777
    if (original["bytes"], original["sha256"], mode) != (
            copied["bytes"], copied["sha256"], destination.stat().st_mode & 0o777):
        raise ValueError("Relocated tool identity differs")
    return dict(original=original, relocated=copied, mode=mode)


def inspect_trace(path, forbidden):
    text = path.read_text()
    if not text or "execve(" not in text:
        raise ValueError("Missing executable trace")
    if any(str(prefix) in text for prefix in forbidden):
        raise ValueError("Original external-tool installation referenced in trace")
    return dict(trace=record(path), forbidden_prefixes=[str(p) for p in forbidden],
                original_prefix_matches=0,
                scope="strace file syscalls, including failed lookups; not OS-level filesystem isolation")


def run(repo, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    prior_record = record(repo / PRIOR)
    if prior_record["sha256"] != PRIOR_SHA:
        raise ValueError("Changed installed phylogeny evidence")
    prior = json.loads((repo / PRIOR).read_text())
    inputs = [prior_record, *prior["checked_records"], *prior["result"]["outputs"], prior["result"]["partition"],
              record(__file__), record(Path(__file__).with_name("verify_frozen_phylogeny_install.py"))]
    for item in inputs:
        check(item)
    mafft = Path(prior["result"]["manifest"]["tools"]["aligner"]["path"])
    fasttree = Path(prior["result"]["manifest"]["tools"]["tree_builder"]["path"])
    mafft_root = mafft.parent.parent
    helpers = mafft_root / "libexec/mafft"
    helper_paths = sorted(helpers.rglob("*"))
    if not helper_paths or any(p.is_symlink() or not (p.is_file() or p.is_dir()) for p in helper_paths):
        raise ValueError("Unexpected MAFFT helper tree")
    tracer = shutil.which("strace")
    if tracer is None:
        raise ValueError("strace unavailable")
    inputs.append(record(tracer))
    output.mkdir(parents=True)
    target_mafft, target_fasttree = output / "tools/bin/mafft", output / "tools/bin/FastTree"
    target_helpers = output / "tools/libexec/mafft"
    copies = [copy_regular(mafft, target_mafft), copy_regular(fasttree, target_fasttree)]
    copies.extend(copy_regular(p, target_helpers / p.relative_to(helpers)) for p in helper_paths if p.is_file())
    copies.extend(copy_regular(mafft_root / name, output / "tools/notices" / name)
                  for name in ("license", "license.extensions"))
    fasta = output / "input"
    for source in sorted(Path(prior["command"][4]).glob("*.fa")):
        copies.append(copy_regular(source, fasta / source.name))
    command = list(prior["command"])
    if command[1:4] != ["-I", "-m", "orthohmm"]:
        raise ValueError("Unexpected installed command shape")
    command[4] = str(fasta)
    inference = output / "inference"
    inference.mkdir()
    for flag, replacement in (("-o", inference), ("--aligner", target_mafft), ("--tree_builder", target_fasttree)):
        command[command.index(flag) + 1] = str(replacement)
    overrides = dict(prior["environment_overrides"], PATH="/usr/bin:/bin", MAFFT_BINARIES=str(target_helpers))
    env = dict(os.environ, **overrides)
    trace = output / "file-access.strace"
    wrapped = [tracer, "-f", "-qq", "-s", "4096", "-e", "trace=%file", "-o", str(trace), *command]
    report = dict(status="relocated_tools_running", checked_records=inputs, copied_files=copies,
        command=wrapped, cwd=str(output), environment_overrides=overrides, attempts=1, timeout_seconds=180,
        accuracy_evaluated=False, redistribution_cleared=False, publication_ready=False,
        limitations=["Local private copy, not a licensed public runtime bundle.",
            "Installed OrthoHMM/Python and OS libraries remain on the same host; not complete runtime relocation.",
            "Traced synthetic fixture only; no biological accuracy or controlled timing claim."])
    try:
        with (output / "inference.log").open("x") as log:
            child = subprocess.Popen(wrapped, cwd=output, env=env, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
            try:
                code = child.wait(timeout=180)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid, signal.SIGKILL)
                child.wait()
                raise
        report["returncode"] = code
        if code:
            raise RuntimeError(f"Relocated tool run failed: {code}")
        report["result"] = validate(inference, fixture())
        report["trace_audit"] = inspect_trace(trace, (mafft_root, fasttree.parent))
        original_phylo = Path(prior["result"]["summary"]["output_directory"])
        current_phylo = Path(report["result"]["summary"]["output_directory"])
        comparisons = []
        pairs = [(Path(prior["result"]["partition"]["path"]), Path(report["result"]["partition"]["path"]))]
        pairs.extend((original_phylo / name, current_phylo / name) for name in (
            "species_tree.rooted.nwk", "orthohmm_pairwise_orthologs.tsv", "orthohmm_root_hogs.tsv"))
        pairs.extend((p, current_phylo / "gene_trees" / p.name) for p in sorted((original_phylo / "gene_trees").glob("*.reconciled.nwk")))
        for original, relocated in pairs:
            if original.read_bytes() != relocated.read_bytes():
                raise ValueError("Relocated scientific output differs: " + original.name)
            comparisons.append(dict(original=record(original), relocated=record(relocated), byte_equal=True))
        report["comparisons"] = comparisons
        for item in [*inputs, *[r[k] for r in copies for k in ("original", "relocated")]]:
            check(item)
        report["status"] = "same_host_external_tools_relocated_and_traced"
    except Exception as error:
        report.update(status="relocated_tools_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        report["log"] = record(output / "inference.log")
        if trace.exists():
            report["trace"] = record(trace)
        with (output / "report.json").open("x") as stream:
            json.dump(report, stream, indent=2, sort_keys=True)
            stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.repo.resolve(), args.output.absolute())
