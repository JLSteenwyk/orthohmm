"""Score the frozen matched-recall graph contrast with paired seed-block uncertainty."""

import argparse
import json
from pathlib import Path

import numpy as np
from Bio import SeqIO

from benchmark_tools.audit_matched_graph import partition
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.simulation_conditions import group_pairs, score_pairs
from benchmark_tools.summarize_simulation_panel import validate_score


READBACK_SHA = "af7a9e0cdba97b7b84c9fc6925e577594802b7bf74c99ac19840baf2e42594f4"
CONDITIONS = ("baseline", "divergent", "divergent_turnover", "missing20", "taxon_count_control", "turnover", "uneven_taxa")
SEEDS = tuple(range(20261106, 20261111))
METRICS = ("f1", "precision", "recall")


def summarize(rows):
    indexed = {}
    for row in rows:
        key = row["condition"], row["seed"], row["arm"]
        if key in indexed:
            raise ValueError("Duplicate score row")
        validate_score(row["score"])
        indexed[key] = row
    expected = {(c, s, a) for c in CONDITIONS for s in SEEDS for a in ("hmm", "diamond")}
    if set(indexed) != expected:
        raise ValueError("Require both arms on all 35 reporting datasets")
    for c in CONDITIONS:
        for s in SEEDS:
            a, b = [indexed[c, s, arm] for arm in ("hmm", "diamond")]
            if a["truth"] != b["truth"] or any(a["score"][k] != b["score"][k] for k in ("input_genes", "eligible_true_pairs")):
                raise ValueError("Paired truth/input universes differ")
    arrays = {arm: np.array([[[indexed[c, s, arm]["score"][m] for m in METRICS]
                             for c in CONDITIONS] for s in SEEDS]) for arm in ("hmm", "diamond")}
    delta = 100 * (arrays["hmm"] - arrays["diamond"])
    rng = np.random.Generator(np.random.PCG64(20260927))
    weights = rng.multinomial(5, np.full(5, .2), size=20000)
    draws = np.einsum("bs,scm->bcm", weights, delta) / 5
    contrasts = {}
    for index, condition in enumerate((*CONDITIONS, "overall")):
        effect = delta.mean(axis=1) if condition == "overall" else delta[:, index, :]
        samples = draws.mean(axis=1) if condition == "overall" else draws[:, index, :]
        methods = {arm: (array.mean(axis=1) if condition == "overall" else array[:, index, :])
                   for arm, array in arrays.items()}
        contrasts[condition] = {}
        for metric_index, metric in enumerate(METRICS):
            values = effect[:, metric_index]
            contrasts[condition][metric] = dict(
                hmm_mean=float(methods["hmm"][:, metric_index].mean()),
                diamond_mean=float(methods["diamond"][:, metric_index].mean()),
                difference_percentage_points=float(values.mean()),
                seed_differences_percentage_points=values.tolist(),
                wins=int((values > 0).sum()), ties=int((values == 0).sum()), losses=int((values < 0).sum()),
                marginal_95_percent_ci=np.quantile(samples[:, metric_index], [.025, .975]).tolist(),
                inference_role="primary_family" if metric == "f1" else "exploratory")
            if metric == "f1":
                contrasts[condition][metric]["bonferroni_8_ci"] = np.quantile(samples[:, metric_index], [.025 / 8, 1 - .025 / 8]).tolist()
    return dict(contrasts=contrasts, seeds=list(SEEDS), conditions=list(CONDITIONS),
                bootstrap=dict(replicates=20000, seed=20260927, rng="PCG64 multinomial",
                               unit="paired seed block carrying all seven conditions", f1_contrasts=8,
                               numpy_version=np.__version__),
                statistic="Equal-weight mean seed pair metrics within condition; equal-condition overall mean",
                direction="HMM minus DIAMOND", records=rows)


def score_cell(cell):
    for item in [cell["truth"], cell["final_partition"], cell["execution"], cell["native_receipt"], *cell["inputs"]]:
        check(item)
    truth = json.loads(Path(cell["truth"]["path"]).read_text())
    owners = {}
    for item in cell["inputs"]:
        for seq in SeqIO.parse(item["path"], "fasta"):
            if seq.id in owners:
                raise ValueError("Duplicate FASTA gene")
            owners[seq.id] = Path(item["path"]).stem
    if len(owners) != cell["genes"] or len(owners) != truth["extant_genes"]:
        raise ValueError("Input/truth gene counts differ")
    groups = partition(Path(cell["final_partition"]["path"]), owners)
    predictions = set(group_pairs(groups, owners))
    score = score_pairs(predictions, truth["ortholog_pairs"], owners)
    validate_score(score)
    # Verify the set scorer against direct membership accounting over every cross-species pair.
    group_ids = {g: i for i, group in enumerate(groups) for g in group}
    true = {frozenset(pair) for pair in truth["ortholog_pairs"]}
    tp = sum(group_ids[a] == group_ids[b] for a, b in truth["ortholog_pairs"])
    fp = sum(frozenset(pair) not in true for pair in predictions)
    if (tp, fp, len(true) - tp) != (score["tp"], score["fp"], score["fn"]):
        raise ValueError("Independent pair count mismatch")
    return dict(condition=cell["condition"], seed=cell["seed"], arm=cell["arm"],
                score=score, truth=cell["truth"], prediction=cell["final_partition"],
                groups=len(groups), gene_coverage=1.0,
                nonsingleton_coverage=sum(len(g) for g in groups if len(g) > 1) / len(owners))


def render(result):
    lines = ["# Matched-Recall Graph Control", "", "All values below are percentages or percentage-point differences.",
             "Intervals are the prespecified Bonferroni-adjusted eight-contrast F1 intervals.", "",
             "| Condition | HMM F1 | DIAMOND F1 | Difference (pp) | Adjusted interval (pp) | Wins/ties/losses |",
             "|---|---:|---:|---:|---|---|"]
    for condition in (*CONDITIONS, "overall"):
        r = result["contrasts"][condition]["f1"]
        lo, hi = r["bonferroni_8_ci"]
        lines.append(f'| {condition} | {100*r["hmm_mean"]:.4f} | {100*r["diamond_mean"]:.4f} | {r["difference_percentage_points"]:+.4f} | [{lo:+.4f}, {hi:+.4f}] | {r["wins"]}/{r["ties"]}/{r["losses"]} |')
    lines += ["", "## Limits", "", *["- " + text for text in result["limitations"]], ""]
    return "\n".join(lines)


def run(readback_path, output):
    if output.exists():
        raise FileExistsError(output)
    readback_record = record(readback_path)
    if readback_record["sha256"] != READBACK_SHA:
        raise ValueError("Unrecognized graph admission")
    readback = json.loads(readback_path.read_text())
    if readback["status"] != "all_native_graphs_verified_pending_orthology_scoring" or readback["accuracy_evaluated"]:
        raise ValueError("Invalid graph admission status")
    check(readback["source"])
    check(readback["submission"])
    rows = [score_cell(cell) for cell in readback["cells"]]
    result = summarize(rows)
    result.update(status="complete_matched_graph_panel_scored", readback=readback_record,
                  source=record(__file__), scoring_dependency=record(Path(__file__).with_name("simulation_conditions.py")),
                  independent_validation=False, scientific_defaults_changed=False, publication_ready=False,
                  limitations=["Development-exposed simulations; not independent confirmation or real-data sensitivity matching.",
                               "Five paired seed blocks give limited tail resolution; percentile intervals have approximate coverage.",
                               "Initial-search graph control only: profile expansion, candidate expansion and phylogeny are off.",
                               "Cluster-derived cross-species pairs are not native reconciled ortholog predictions.",
                               "Matching search recall does not equalize score distributions, computational effort or hit identities.",
                               "Shared-host incremental resources are descriptive, not controlled efficiency evidence."])
    for cell in readback["cells"]:
        for item in [cell["truth"], cell["final_partition"], *cell["inputs"]]:
            check(item)
    output.mkdir(parents=True)
    with (output / "results.json").open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    with (output / "results.md").open("x") as stream:
        stream.write(render(result))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--readback", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.readback.resolve(), args.output.absolute())
