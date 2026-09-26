"""Verify retained baseline FastTree binaries without rebuilding or replacing them."""

import argparse
import json
from pathlib import Path
import subprocess

from benchmark_tools.acquire_publication_fasttree import identity, verify
from benchmark_tools.build_publication_fasttree import validate_help
from benchmark_tools.build_publication_mafft import execute
from benchmark_tools.audit_phylogeny_structure import audit as structure_audit
from benchmark_tools.audit_phylogeny_sequences import audit as sequence_audit
from benchmark_tools.audit_phylogeny_events import audit as event_audit
from benchmark_tools.audit_phylogeny_hierarchy import audit as hierarchy_audit


def run(repo, build_report, mafft, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    prior = json.loads(build_report.read_text())
    if (prior["status"] != "failed" or prior.get("error") != "Unexpected built FastTree identity/precision"
            or len(prior["builds"]) != 2 or not prior["repeat_build_byte_equal"]):
        raise ValueError("Expected retained help-validator failure after two successful builds")
    checked = [identity(build_report), *prior["sources"], *prior["builds"], prior["compiler"]["binary"],
               prior["linker"], identity(mafft)]
    for item in checked:
        verify(Path(item["path"]), item["sha256"])
    output.mkdir(parents=True)
    report = dict(status="verifying", checked_records=checked, source=identity(__file__),
                  helper=identity(Path(__file__).with_name("build_publication_fasttree.py")),
                  builds_repeated=False, historical_binary_modified=False, publication_ready=False)
    env = dict(prior["compile_environment"], MAFFT_BINARIES=str(mafft.parents[1] / "libexec/mafft"))
    binary = Path(prior["builds"][0]["path"])
    try:
        probe = subprocess.run([str(binary), "-help"], env=env, capture_output=True, text=True, timeout=30)
        validate_help(probe.returncode, probe.stdout, probe.stderr)
        report["help_probe"] = dict(returncode=probe.returncode, stdout=probe.stdout, stderr=probe.stderr)
        report["elf"] = subprocess.check_output(["/usr/bin/readelf", "-h", "-n", "-d", str(binary)], text=True, env=env)
        report["runtime_libraries"] = subprocess.check_output(["/usr/bin/ldd", str(binary)], text=True, env=env)
        fixture = output / "phylogeny"
        command = ["/usr/bin/python3", "-S", "-m", "benchmark_tools.verify_frozen_phylogeny_install",
            "--repo", str(repo), "--artifact", str(repo / "benchmarks/work/publication_frozen_overlay_20260926"),
            "--mafft", str(mafft), "--fasttree", str(binary), "--output", str(fixture)]
        report.update(command=command, environment=env)
        execute(command, repo, env, output / "fixture_driver.log", timeout=240)
        report["fixture_report"] = identity(fixture / "report.json")
        phylo = fixture / "inference/orthohmm_phylogeny"
        paths = {name: output / f"{name}.json" for name in ("structure", "sequence", "events", "hierarchy")}
        results = {}
        for name in paths:
            if name == "structure":
                result = structure_audit(phylo, fixture / "input")
            elif name == "sequence":
                result = sequence_audit(phylo, paths["structure"])
            elif name == "events":
                result = event_audit(phylo, paths["structure"])
            else:
                result = hierarchy_audit(phylo, paths["events"])
            with paths[name].open("x") as stream:
                json.dump(result, stream, indent=2, sort_keys=True)
                stream.write("\n")
            results[name] = dict(status=result["status"], report=identity(paths[name]))
        report["readbacks"] = results
        baseline = repo / "benchmarks/work/publication_frozen_phylogeny_smoke_20260926/inference"
        current = fixture / "inference"
        names = ["orthohmm_orthogroups.txt", "orthohmm_phylogeny/orthohmm_root_hogs.tsv",
                 "orthohmm_phylogeny/orthohmm_pairwise_orthologs.tsv",
                 "orthohmm_phylogeny/species_tree.rooted.nwk",
                 "orthohmm_phylogeny/gene_trees/Family0000002.reconciled.nwk"]
        report["comparisons"] = [dict(relative=name, baseline=identity(baseline / name), current=identity(current / name),
            byte_equal=(baseline / name).read_bytes() == (current / name).read_bytes()) for name in names]
        for item in checked:
            verify(Path(item["path"]), item["sha256"])
        report["status"] = "retained_baseline_fasttree_fixture_and_semantic_readbacks_verified"
    except Exception as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        report["logs"] = [identity(p) for p in sorted(output.glob("*.log"))]
        with (output / "report.json").open("x") as stream:
            json.dump(report, stream, indent=2, sort_keys=True)
            stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--build-report", type=Path, required=True)
    parser.add_argument("--mafft", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.repo.resolve(), args.build_report.resolve(), args.mafft.resolve(), args.output.absolute())
