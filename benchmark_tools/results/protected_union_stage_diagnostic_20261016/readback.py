"""Read-only verification recipe for the rejected union pilot's stage diagnostic."""

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

from Bio import Phylo

REPO = Path(__file__).resolve().parents[3]
PROTOCOL = REPO / "benchmark_tools/hillclimb_protected_union_stage_diagnostic_20261016.json"
PLAN = json.loads(PROTOCOL.read_text())
FROZEN = REPO / PLAN["executor_source"]
sys.path.insert(0, str(FROZEN))

from benchmark_tools.bootstrap_recent_hillclimb_ob import selected
from benchmark_tools.compare_native_membership import partition_contrast
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_full_search_hillclimb import runtime
from benchmark_tools.run_membership_hillclimb import score_orthobench, write_json
from benchmark_tools.score_ygob_groups import membership, read_predictions
from benchmark_tools.summarize_matched_resources import verbose_time


def packed_pairs(path, index, owners):
    values, count = set(), 0
    with path.open() as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        assert reader.fieldnames == ["gene_a", "species_a", "gene_b", "species_b"]
        for row in reader:
            assert None not in row and all(v and v.strip() == v for v in row.values())
            g, h, s, t = (row[k] for k in ("gene_a", "gene_b", "species_a", "species_b"))
            assert g != h and s != t and g in index and h in index
            assert owners.setdefault(g, s) == s and owners.setdefault(h, t) == t
            a, b = sorted((index[g], index[h]))
            values.add((a << 32) | b)
            count += 1
    assert len(values) == count
    return values


def topology(path):
    tree = Phylo.read(path, "newick")
    leaves = [node.name for node in tree.get_terminals()]
    assert len(leaves) == len(set(leaves)) == 12
    taxa = set(leaves)
    clades = {frozenset(x.name for x in node.get_terminals()) for node in tree.get_nonterminals()}
    rooted = {group for group in clades if 1 < len(group) < 12}
    unrooted = {min(tuple(sorted(group)), tuple(sorted(taxa - group)))
                for group in clades if 1 < len(group) < 11}
    return taxa, rooted, unrooted


def main():
    out = REPO / PLAN["output"]
    assert not (out / "readback.json").exists() and not (out / "partition_contrasts.json").exists()
    assert subprocess.check_output(["git", "show", "eebc4364:benchmark_tools/" + PROTOCOL.name]) == PROTOCOL.read_bytes()
    accounting = subprocess.check_output(["sacct", "-j", "24379", "-X", "-P", "-o",
        "JobIDRaw,State,ExitCode,AllocCPUS,NodeList,Elapsed,ReqMem,Timelimit,AllocTRES"], text=True)
    assert "24379|COMPLETED|0:0|1|bizon|00:00:19|4G|00:15:00|billing=1,cpu=1,mem=4G,node=1" in accounting
    summary_path = REPO / "benchmark_tools/results/hillclimb_comparison_readback_v27_20261016/summary.json"
    summary = json.loads(summary_path.read_text())
    cells = {name: selected(summary, key) for name, key in (
        ("control", "canonical_collection:canonical"), ("union", "protected_union:protected_union"))}
    pins = [record(__file__), record(PROTOCOL), record(summary_path), summary["score_table"],
            record(REPO / PLAN["all_family_manifest"]["path"])]
    scores, official = {}, {}
    for label, pin in PLAN["partitions"].items():
        path = REPO / pin["path"]
        identity = record(path)
        assert all(identity[k] == pin[k] for k in ("bytes", "sha256"))
        raw_path = out / (label + ".json")
        raw = json.loads(raw_path.read_text())
        assert raw["partition"] == "all" and raw["inputs"]["predictions"] == identity
        check(raw["source"])
        assert Path(raw["source"]["path"]) == FROZEN / "benchmark_tools/score_orthobench_partition.py"
        reread = score_orthobench(REPO, path)
        assert reread["details"] == raw["score"] and len(raw["score"]["refog_records"]) == 70
        scores[label], official[label] = raw["score"], reread["official"]
        pins += [identity, record(raw_path), *reread["references"]]
        print(label, {k: scores[label][k] for k in ("f_score", "precision", "recall")}, flush=True)
    for label, cell in cells.items():
        assert scores[label + "_final"] == cell["details"]["family_details"]
        pins += cell["checked"]
    pairs = [("between_candidates", "control_candidate", "union_candidate"),
             ("between_finals", "control_final", "union_final"),
             ("control_candidate_to_final", "control_candidate", "control_final"),
             ("union_candidate_to_final", "union_candidate", "union_final")]
    changes = {label: {metric: scores[b][metric] - scores[a][metric]
                      for metric in ("f_score", "precision", "recall")} for label, a, b in pairs}
    for metric in changes["between_finals"]:
        assert abs(changes["between_finals"][metric] - (changes["between_candidates"][metric]
            + changes["union_candidate_to_final"][metric] - changes["control_candidate_to_final"][metric])) < 1e-12
    universe = set(membership(read_predictions(
        Path(cells["control"]["products"]["orthohmm_root_hogs.tsv"]["path"]), "root_hogs")))
    assert len(universe) == 251378
    contrasts = {label: partition_contrast(REPO / PLAN["partitions"][a]["path"],
        REPO / PLAN["partitions"][b]["path"], universe) for label, a, b in pairs}
    native_path = out / "native_contrast.json"
    native = json.loads(native_path.read_text())
    assert native["genes"] == 251378
    assert {k: v for k, v in contrasts["between_finals"].items()
            if k not in ("genes", "checked")} == native["root_membership"]
    pins += [record(native_path), *native["checked"]]
    index, owners = {gene: i for i, gene in enumerate(sorted(universe))}, {}
    sets = [packed_pairs(Path(cells[label]["products"]["orthohmm_pairwise_orthologs.tsv"]["path"]),
                         index, owners) for label in ("control", "union")]
    counts = dict(baseline_records=len(sets[0]), candidate_records=len(sets[1]),
                  common_records=len(sets[0] & sets[1]), removed_records=len(sets[0] - sets[1]),
                  added_records=len(sets[1] - sets[0]))
    assert all(native["pair_membership"][k] == v for k, v in counts.items())
    assert all(len(values) == cells[label]["coverage"]["pairs"]
               for values, label in zip(sets, ("control", "union")))
    affected = {i for value in sets[0] ^ sets[1] for i in (value >> 32, value & ((1 << 32) - 1))}
    assert len(affected) == native["pair_membership"]["affected_genes"]
    tree_paths = [Path(cells[label]["products"]["orthohmm_root_hogs.tsv"]["path"]).parent
                  / "species_tree.rooted.nwk" for label in ("control", "union")]
    trees = [topology(path) for path in tree_paths]
    assert trees[0][0] == trees[1][0]
    tree_result = dict(taxa=sorted(trees[0][0]), rooted_clade_counts=[len(t[1]) for t in trees],
        unrooted_bipartition_counts=[len(t[2]) for t in trees],
        rooted_symmetric_difference=len(trees[0][1] ^ trees[1][1]),
        unrooted_symmetric_difference=len(trees[0][2] ^ trees[1][2]),
        same_rooted_topology=trees[0][1] == trees[1][1], same_unrooted_topology=trees[0][2] == trees[1][2])
    producer = json.loads((REPO / "benchmarks/results/protected_union_orthobench_r1_v1_20261016/results.json").read_text())
    execution = json.loads(Path(producer["execution_plan"]["path"]).read_text())
    assert runtime() == execution["runtime"] and cells["union"]["gate"] == producer["gate"]
    assert producer["gate"]["passes_development_gate"] is False
    pins += [*execution["sources"], *execution["inputs"], *execution["runtime"]["files"],
             *summary["reporter_files"], *[record(path) for path in tree_paths],
             record(REPO / "benchmarks/results/protected_union_stage_preparation_tests_20261016.xml")]
    for name, pin in PLAN["helpers"].items():
        path = FROZEN / "benchmark_tools" / name
        assert path.stat().st_size == pin["bytes"] and hashlib.sha256(path.read_bytes()).hexdigest() == pin["sha256"]
        assert path.read_bytes() == subprocess.check_output([
            "git", "-C", str(FROZEN), "show", PLAN["executor_revision"] + ":benchmark_tools/" + name])
    pins = list({pin["path"]: pin for pin in pins}.values())
    for pin in pins:
        check(pin)
    assert runtime() == execution["runtime"]
    write_json(out / "partition_contrasts.json", contrasts)
    result = dict(status="protected_union_stage_readback_independently_passed", recipe=record(__file__),
        protocol=record(PROTOCOL), executor_revision=PLAN["executor_revision"], scheduler=accounting,
        unique_postflight_pins=len(pins), changes_percentage_points=changes, official_scores=official,
        score_points_percent={label: {k: score[k] for k in ("f_score", "precision", "recall", "exact_refogs")}
                              for label, score in scores.items()},
        native_pair_counts=counts, affected_pair_genes=len(affected), species_tree_topology=tree_result,
        artifacts=[record(path) for path in [*(out / (k + ".json") for k in PLAN["partitions"]),
            native_path, out / "partition_contrasts.json", REPO / "benchmarks/results/protected_union_stage_diagnostic_20261016.time.txt",
            REPO / "benchmarks/results/protected_union_stage_diagnostic_20261016-24379.out"]],
        resources=verbose_time((REPO / "benchmarks/results/protected_union_stage_diagnostic_20261016.time.txt").read_text()),
        preparation_tests=record(REPO / "benchmarks/results/protected_union_stage_preparation_tests_20261016.xml"),
        readback_setup_failures=["Initial ad hoc readback passed None instead of a complete universe to strict read_partition; it failed before any derived file was written.",
            "Second ad hoc readback assumed the tree was in the control products dictionary; it failed after score/pair checks but before any derived file was written. Corrected to the native root-products directory's actual tree file."],
        scope="Read-only post-hoc stage association. All 70 families and 11 low-certainty files retained; original table/gate unchanged. CLI Git metadata describes its working directory; frozen source identity is separately authenticated. Diagnostic resource record excludes this independent readback. No unique causal attribution, inference rerun, parameter retuning or independent confirmation.",
        defaults_promoted=False, goal_completed=False)
    write_json(out / "readback.json", result)
    print(json.dumps({k: result[k] for k in ("status", "unique_postflight_pins", "changes_percentage_points",
        "native_pair_counts", "species_tree_topology", "resources")}, indent=2))


if __name__ == "__main__":
    main()
