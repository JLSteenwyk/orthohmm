"""Read back MAFFT build identities and compare the fixed phylogeny fixture."""

import argparse
import json
from pathlib import Path

from benchmark_tools.acquire_publication_fasttree import identity, verify
from benchmark_tools.probe_relocated_phylogeny_tools import PRIOR, PRIOR_SHA


def compare_files(left, right):
    a, b = identity(left), identity(right)
    return dict(original=a, rebuilt=b, byte_equal=(a["sha256"], a["bytes"]) == (b["sha256"], b["bytes"]))


def run(repo, build, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    current = json.loads(build.read_text())
    if current["status"] != "core_built_and_installed_phylogeny_fixture_passed":
        raise ValueError("Build/fixture did not pass")
    prior_record = verify(repo / PRIOR, PRIOR_SHA)
    prior = json.loads((repo / PRIOR).read_text())
    new = json.loads(Path(current["phylogeny_report"]["path"]).read_text())
    records = [identity(build), identity(__file__), prior_record, current["archive"], current["script"],
               current["launcher"], current["phylogeny_report"], current["compiler"]["file"],
               *current["source_files"], *current["built_helpers"], *current["retained_helpers"],
               *current["logs"], *prior["checked_records"], *new["checked_records"],
               *new["result"]["outputs"], new["result"]["partition"], new["log"],
               *[row["retained"] for row in current["retained_source_comparisons"]]]
    for record in records:
        observed = verify(Path(record["path"]), record["sha256"])
        if observed["bytes"] != record["bytes"]:
            raise ValueError("Recorded byte count changed")
    previous = {Path(r["path"]).name: Path(r["path"]) for r in current["retained_helpers"]}
    rebuilt = {Path(r["path"]).name: Path(r["path"]) for r in current["built_helpers"]}
    if set(previous) != set(rebuilt) or not previous:
        raise ValueError("Helper inventories differ")
    helpers = [compare_files(previous[name], rebuilt[name]) for name in sorted(previous)]
    old_dir = Path(prior["result"]["summary"]["output_directory"])
    new_dir = Path(new["result"]["summary"]["output_directory"])
    comparisons = [compare_files(prior["result"]["partition"]["path"], new["result"]["partition"]["path"])]
    comparisons.extend(compare_files(old_dir / name, new_dir / name) for name in (
        "species_tree.rooted.nwk", "orthohmm_pairwise_orthologs.tsv", "orthohmm_root_hogs.tsv",
        "gene_trees/Family0000002.reconciled.nwk"))
    report = dict(status="readback_complete", checked_records=records, helper_comparisons=helpers,
                  fixture_comparisons=comparisons, all_helpers_identical=all(r["byte_equal"] for r in helpers),
                  all_selected_fixture_outputs_identical=all(r["byte_equal"] for r in comparisons),
                  compiler_warning_occurrences=Path(build.parent / "build.log").read_text().count("warning:"),
                  scope="Same-host source build and synthetic fixture; not whole-dataset or cross-platform equivalence.")
    with output.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("repo", "build", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    run(args.repo.resolve(), args.build.resolve(), args.output.absolute())
