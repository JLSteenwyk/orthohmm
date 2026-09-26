"""Exercise inferred phylogeny in the separately installed frozen-source wheel."""

import argparse
import csv
import json
import os
from pathlib import Path
import random
import signal
import subprocess

from benchmark_tools.prepare_ob_candidate_neighborhood import record, check

AUDIT_SHA = "467036ee88d5e6f15d322aea4c94397591ef7ba521ac9f8ddde4b7cb783b97d2"
SEED = 20260926
ALPHABET = "ACDEFGHIKLMNPQRSTVWY"


def fixture():
    rng = random.Random(SEED)
    ancestors = ["".join(rng.choice(ALPHABET) for _ in range(200)) for _ in range(3)]
    result = {}
    for taxon in ("S1", "S2", "S3", "S4"):
        genes = {}
        for family, ancestor in enumerate(ancestors):
            sequence = list(ancestor)
            for position in rng.sample(range(200), 10):
                sequence[position] = rng.choice(ALPHABET.replace(sequence[position], ""))
            sequence = "".join(sequence)
            genes[f"{taxon}_family{family}"] = sequence
            if family == 2:
                genes[f"{taxon}_family{family}_duplicate"] = sequence
        result[taxon] = genes
    return result


def validate(directory, genes):
    manifests = list(directory.rglob("provenance_manifest.json"))
    if len(manifests) != 1:
        raise ValueError("Expected one phylogeny manifest")
    phylo = manifests[0].parent
    manifest = json.loads(manifests[0].read_text())
    summary = json.loads((phylo / "reconciliation_summary.json").read_text())
    if (manifest["species_tree_mode"] != "infer" or manifest["species_tree_rooting"] != "min_variance"
            or manifest["root_duplication_rule"] != "species_overlap"
            or manifest["pair_orthology_rule"] != "positive_paralogy"
            or set(manifest["species_tree_taxa"]) != set(genes)
            or summary["species_tree_families"] < 1 or summary["reconciled_families"] < 1
            or summary["checkpoint_hits"] != 0 or summary["species_tree_checkpoint_hit"] is not False):
        raise ValueError("Required fresh inferred phylogeny/reconciliation did not execute")
    universe = {gene for members in genes.values() for gene in members}
    assigned = []
    for line in (directory / "orthohmm_orthogroups.txt").read_text().splitlines():
        label, members = line.split(":", 1)
        if not label or not members.split():
            raise ValueError("Empty output group")
        assigned.extend(members.split())
    if len(assigned) != len(universe) or set(assigned) != universe:
        raise ValueError("Output does not cover input exactly once")
    pairs = set()
    with (phylo / "orthohmm_pairwise_orthologs.tsv").open() as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            a, b = row["gene_a"], row["gene_b"]
            sa, sb = row["species_a"], row["species_b"]
            pair = tuple(sorted((a, b)))
            if sa == sb or a not in genes.get(sa, {}) or b not in genes.get(sb, {}) or pair in pairs:
                raise ValueError("Invalid or duplicate cross-species pair")
            pairs.add(pair)
    if not pairs or len(pairs) != summary["ortholog_pairs"]:
        raise ValueError("Pair output/summary differs")
    trees = sorted(phylo.glob("gene_trees/*.reconciled.nwk"))
    if len(trees) != summary["reconciled_families"] or any(not p.read_text().strip() for p in trees):
        raise ValueError("Missing reconciled tree artifacts")
    return dict(manifest=manifest, summary=summary, genes=len(universe), pairs=len(pairs),
                outputs=[record(p) for p in sorted(phylo.rglob("*")) if p.is_file()],
                partition=record(directory / "orthohmm_orthogroups.txt"))


def run(repo, artifact, mafft, fasttree, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    audit_path = repo / "benchmark_tools/results/publication_frozen_overlay_install_20260926.json"
    identity = record(audit_path)
    if identity["sha256"] != AUDIT_SHA:
        raise ValueError("Changed installed-source evidence")
    audit = json.loads(audit_path.read_text())
    inputs = [identity, *audit["checked_records"], record(mafft), record(fasttree), record(__file__)]
    for item in inputs:
        check(item)
    python = artifact / "venv_clean/bin/python"
    probe = json.loads(subprocess.check_output([str(python), "-I", "-c",
        "import json,sys,orthohmm;print(json.dumps({'module':orthohmm.__file__,'prefix':sys.prefix}))"], text=True, cwd=artifact))
    module = Path(probe["module"]).resolve()
    if not module.is_relative_to(artifact / "venv_clean") or record(module) not in audit["checked_records"]:
        raise ValueError("Unexpected installed package")
    output.mkdir(parents=True)
    fasta = output / "input"
    fasta.mkdir()
    genes = fixture()
    for taxon, members in genes.items():
        (fasta / f"{taxon}.fa").write_text("".join(f">{name}\n{seq}\n" for name, seq in members.items()))
    inputs.extend(record(p) for p in sorted(fasta.iterdir()))
    destination = output / "inference"
    destination.mkdir()
    command = [str(python), "-I", "-m", "orthohmm", str(fasta), "-o", str(destination),
        "-c", "1", "--threads_per_worker", "1", "--search_mode", "builtin", "--clustering", "leiden",
        "--accuracy_profile", "high_sensitivity", "--phylogeny", "reconcile", "--species_tree_mode", "infer",
        "--species_tree_rooting", "min_variance", "--phylogeny_candidates", "satellite_v2",
        "--phylogeny_root_rule", "species_overlap", "--phylogeny_pair_rule", "positive_paralogy",
        "--aligner", str(mafft), "--tree_builder", str(fasttree)]
    overrides = dict(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", PYTHONHASHSEED="0")
    env = dict(os.environ, **overrides)
    report = dict(status="installed_phylogeny_running", command=command, cwd=str(output), checked_records=inputs,
        fixture_seed=SEED, environment_overrides=overrides, imported_package=probe, attempts=1,
        accuracy_evaluated=False, publication_ready=False, timeout_seconds=180,
        limitations=["Synthetic plumbing fixture, not biological accuracy, full benchmark or historical runtime equivalence.",
            "External MAFFT/FastTree remain local installations, not bundled or fully transitive runtime closure.",
            "Same host, not cross-platform portability or controlled timing."])
    log = output / "inference.log"
    try:
        with log.open("x") as stream:
            child = subprocess.Popen(command, cwd=output, env=env, stdout=stream, stderr=subprocess.STDOUT, start_new_session=True)
            try:
                code = child.wait(timeout=180)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid, signal.SIGKILL)
                child.wait()
                raise
        report["returncode"] = code
        if code:
            raise RuntimeError(f"Installed phylogeny exited {code}; no retry")
        report["result"] = validate(destination, genes)
        for item in inputs:
            check(item)
        report["status"] = "installed_inferred_phylogeny_fixture_verified"
    except Exception as error:
        report.update(status="installed_phylogeny_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        if log.exists():
            report["log"] = record(log)
        with (output / "report.json").open("x") as stream:
            json.dump(report, stream, indent=2, sort_keys=True)
            stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("repo", "artifact", "mafft", "fasttree", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    run(args.repo.resolve(), args.artifact.resolve(), args.mafft.resolve(), args.fasttree.resolve(), args.output.absolute())
