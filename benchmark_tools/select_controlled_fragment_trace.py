"""Freeze one deterministic representative per fixed fragment diagnostic bin."""

import argparse
import hashlib
import json
from pathlib import Path

from Bio import SeqIO
from benchmark_tools.export_native_qfo_three_cell_strata import record, require, write_tsv
from benchmark_tools import simulation_method_outputs as native


PINS = {
    "result": ("benchmark_tools/results/controlled_fragment_results_20261009_v1/report.json",
               "b1480e8eb2ace2aed743f5fc4e515e27ac4ba68b48e644b33989eda3cf8a7ab7"),
    "reader": ("benchmark_tools/results/controlled_fragment_execution_20261009_v1.json",
               "af665df1fff1ed454c9cce0286d9bcac831e935ad50b65aba5ee14e76d58edcf"),
    "panel": ("benchmarks/work/controlled_fragment_observations_20261009_v1/manifest.json",
              "8e58830911f60de5aec9ad616e5c02c12229702780e26216aa1371017637666d"),
    "protocol": ("benchmark_tools/results/CONTROLLED_FRAGMENT_TRACE_PROTOCOL_20261009.md",
                 "bb5fadb866b33ec6131e61cc3159695aaa3d9dac483d33851610904c20694a9b"),
}
SEEDS = tuple(range(20261101, 20261111))
METHODS = ("orthohmm_high_sensitivity", "orthohmm_satellite_v2",
           "orthofinder_full", "orthofinder_sequence_only")
ARMS = ("baseline", "fragment")
CATEGORIES = ("fragment_fn", "fragment_fp", "new_fn", "new_fp", "recovered_fn", "removed_fp", "retained_tp")
BIN_FIELDS = ("method", "fragment_endpoints", "category", "eligible_records", "status", "seed", "gene_a", "gene_b", "selection_sha256")


def checked(ref, evidence):
    path = Path(ref.get("absolute_path") or ref["path"]).resolve()
    actual = record(path)
    require(actual["bytes"] == ref["bytes"] and actual["sha256"] == ref["sha256"], "Changed trace input: " + str(path))
    require(str(path) not in evidence or evidence[str(path)] == actual, "Conflicting trace input identity")
    evidence[str(path)] = actual
    return path


def categories(truth, baseline, fragment):
    return dict(fragment_fn=truth - fragment, fragment_fp=fragment - truth,
        new_fn=(truth & baseline) - fragment, new_fp=(fragment - truth) - baseline,
        recovered_fn=(truth & fragment) - baseline, removed_fp=(baseline - truth) - fragment,
        retained_tp=truth & baseline & fragment)


def selection_key(method, category, seed, pair):
    text = f"{method}:{category}:{seed}:{pair[0]}:{pair[1]}"
    return hashlib.sha256(text.encode("utf-8")).hexdigest(), seed, pair[0], pair[1]


def select(datasets):
    require(len(datasets) == 10 and sorted(d["seed"] for d in datasets) == list(SEEDS), "Wrong ten-seed selection inventory")
    buckets = {(m, n, c): [] for m in METHODS for n in range(3) for c in CATEGORIES}
    for data in sorted(datasets, key=lambda d: d["seed"]):
        flags, owners, truth = (data[k] for k in ("flags", "owners", "truth"))
        require(set(flags) == set(owners) and all(type(v) is bool for v in flags.values()), "Invalid trace flags")
        require(set(data["predictions"]) == set(ARMS)
                and all(set(data["predictions"][a]) == set(METHODS) for a in ARMS), "Incomplete trace methods")
        for method in METHODS:
            baseline, fragment = (data["predictions"][a][method] for a in ARMS)
            for pair in truth | baseline | fragment:
                require(type(pair) is tuple and len(pair) == 2 and pair[0] < pair[1]
                        and set(pair) <= set(owners) and owners[pair[0]] != owners[pair[1]], "Invalid canonical trace pair")
            for category, pairs in categories(truth, baseline, fragment).items():
                for pair in pairs:
                    endpoints = sum(flags[g] for g in pair)
                    buckets[method, endpoints, category].append(selection_key(method, category, data["seed"], pair))
    bins, cases = [], []
    for (method, endpoints, category), records in buckets.items():
        require(len(records) == len(set(records)), "Repeated eligible trace identity")
        row = dict(method=method, fragment_endpoints=endpoints, category=category,
                   eligible_records=len(records), status="selected" if records else "empty_bin")
        if records:
            digest, seed, a, b = min(records)
            row.update(seed=seed, gene_a=a, gene_b=b, selection_sha256=digest)
            cases.append(dict(case_id=f"Case{len(cases):04d}", **row))
        else:
            row.update(seed=None, gene_a=None, gene_b=None, selection_sha256=None)
        bins.append(row)
    return bins, cases


def load(root):
    evidence, docs = {}, {}
    for key, (name, digest) in PINS.items():
        ref = record(root / name)
        require(ref["sha256"] == digest, "Changed frozen trace input: " + key)
        path = checked(ref, evidence)
        docs[key] = path.read_text() if key == "protocol" else json.loads(path.read_text())
    report, reader, panel = (docs[k] for k in ("result", "reader", "panel"))
    require(report["schema"] == "controlled_fragment_results_v1"
            and report["status"] == "complete_with_explicit_outcomes" and report["publication_ready"] is False,
            "Wrong retained fragment results")
    require(checked(report["manifest"], evidence) == root / PINS["panel"][0], "Mixed panel binding")
    terminal = reader["independent_result_readback"]["terminal"]
    observed = json.loads(terminal["output"])
    require(terminal["exit_code"] == 0 and observed["report_sha256"] == PINS["result"][1]
            and observed["status"] == "independent_counts_intervals_and_tables_verified"
            and observed["score_records"] == 80 and observed["stratum_records"] == 240,
            "Missing complete independent result readback")
    require(panel["schema"] == "controlled_fragment_observations_v1"
            and panel["condition"] == "fragment20_center60_v1"
            and [d["seed"] for d in panel["datasets"]] == list(SEEDS), "Changed prepared panel")
    parent_ref = [r for r in panel["baseline_pins"] if Path(r["absolute_path"]).name == "publication_variable_native_methods_20260916.json"]
    require(len(parent_ref) == 1, "Ambiguous retained baseline configuration")
    parent = json.loads(checked(parent_ref[0], evidence).read_text())
    baseline = {d["seed"]: d for d in parent["datasets"] if d["condition"] == "baseline"}
    require(set(baseline) == set(SEEDS), "Incomplete baseline configurations")
    records = {(r["arm"], r["seed"], r["method"]): r for r in report["records"]}
    require(len(records) == len(report["records"]) == 80
            and set(records) == {(a, s, m) for a in ARMS for s in SEEDS for m in METHODS}
            and all(r["status"] == "complete" for r in records.values()), "Incomplete retained score inventory")
    data = []
    for dataset in panel["datasets"]:
        seed = dataset["seed"]
        truth_path = checked(dataset["verified_inputs"]["truth"], evidence)
        require(str(truth_path) == dataset["truth"] and dataset["verified_inputs"]["truth"]["sha256"]
                == dataset["parent_inputs"]["truth"]["sha256"], "Changed fragment truth")
        truth_doc = json.loads(truth_path.read_text())
        owners = {}
        for ref in dataset["verified_inputs"]["inputs"]:
            path = checked(ref, evidence)
            for sequence in SeqIO.parse(path, "fasta"):
                require(sequence.id not in owners, "Duplicate input gene")
                owners[sequence.id] = path.stem
        coordinate_rows = json.loads(checked(dataset["coordinates"], evidence).read_text())
        flags = {r["gene"]: r["fragment"] for r in coordinate_rows}
        require(len(flags) == len(coordinate_rows) == len(owners) and set(flags) == set(owners), "Changed coordinate universe")
        ancestors = {}
        for family, genes in truth_doc["families"].items():
            for gene in genes:
                require(gene not in ancestors, "Repeated truth-family gene")
                ancestors[gene] = family
        require(set(ancestors) == set(owners), "Truth families differ from input universe")
        truth = {tuple(sorted(p)) for p in truth_doc["ortholog_pairs"]}
        require(len(truth) == len(truth_doc["ortholog_pairs"]), "Duplicate truth pair")
        row = dict(seed=seed, owners=owners, ancestors=ancestors, flags=flags, truth=truth,
                   predictions={a: {} for a in ARMS}, bindings={a: {} for a in ARMS})
        original = {r["method"]: r for r in dataset["baseline_admission"]}
        require(len(original) == len(dataset["baseline_admission"]) == 4, "Incomplete original admission binding")
        for arm in ARMS:
            config = baseline[seed] if arm == "baseline" else dataset
            for method in METHODS:
                outcome = records[arm, seed, method]
                retained = original[method] if arm == "baseline" else outcome
                execution_path = checked(retained["execution"], evidence)
                execution = json.loads(execution_path.read_text())
                parent_method = "orthofinder_full" if method == METHODS[3] else method
                executed = execution["methods"][parent_method]
                require(executed["status"] == "process_succeeded" and executed["exit_code"] == 0
                        and executed["argv"] == config["methods"][parent_method]["argv"], "Changed native method execution")
                artifacts = {str(Path(r["absolute_path"]).resolve()): r for r in executed["outputs"]}
                require(len(artifacts) == len(executed["outputs"]), "Duplicate executed artifacts")
                predictions_ref = retained["predictions" if arm == "baseline" else "prediction_artifacts"]
                for ref in predictions_ref:
                    path = checked(ref, evidence)
                    require(str(path) in artifacts and all(ref[k] == artifacts[str(path)][k] for k in ("bytes", "sha256")),
                            "Prediction absent from executed native inventory")
                output = Path(config["methods"][method]["output"]).resolve()
                pairs, paths = native.load_predictions(method, output, owners, sorted(set(owners.values())))
                require({str(p.resolve()) for p in paths} == {str(Path(r["absolute_path"]).resolve()) for r in predictions_ref},
                        "Native adapter found different prediction artifacts")
                canonical = {tuple(sorted(p)) for p in pairs}
                expected = outcome["score"]
                require((len(canonical & truth), len(canonical - truth), len(truth - canonical))
                        == (expected["tp"], expected["fp"], expected["fn"]), "Retained native counts do not reproduce")
                row["predictions"][arm][method] = canonical
                row["bindings"][arm][method] = dict(output=str(output), execution=record(execution_path),
                    configured=config["methods"][method], predicted_pairs=len(canonical))
        data.append(row)
    for path in (__file__, native.__file__):
        checked(record(path), evidence)
    return data, evidence


def run(root, output):
    root, output = Path(root).resolve(), Path(output).absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    require(output.parent.resolve() == root / "benchmark_tools/results", "Keep trace artifacts in established results")
    data, evidence = load(root)
    bins, cases = select(data)
    index = {d["seed"]: d for d in data}
    for row in cases:
        source = index[row["seed"]]
        pair = row["gene_a"], row["gene_b"]
        row.update(truth=pair in source["truth"], species=[source["owners"][g] for g in pair],
            ancestral_families=[source["ancestors"][g] for g in pair],
            fragment_flags=[source["flags"][g] for g in pair],
            comparator_predictions={a: {m: pair in source["predictions"][a][m] for m in METHODS} for a in ARMS})
    distinct = {(r["method"], r["seed"], r["gene_a"], r["gene_b"]) for r in cases}
    report = dict(schema="controlled_fragment_trace_selection_v1", protocol=record(root / PINS["protocol"][0]),
        source=record(__file__), status="selected_untraced", planned_bins=84, bins=bins, cases=cases,
        selected_records=len(cases), unique_method_pair_records=len(distinct), repeated_selected_records=len(cases)-len(distinct),
        bindings=[dict(seed=d["seed"], arms=d["bindings"]) for d in data], checked_inputs=list(evidence.values()),
        new_inference_or_scoring=False, new_bootstrap_draws=0, uncertainty_admitted=False,
        scientific_timings_admitted=False, independent_confirmation=False, publication_ready=False,
        limitations=["Retrospective deterministic representatives after aggregate exposure; not prevalence or independent effects.",
            "Empty categories retained; same pair may represent multiple categories.",
            "Initial HMM search on; missing significant hits cannot distinguish prefilter rejection from scoring.",
            "Synthetic fragments, changed full OrthoHMM inventory and shared-host limitations remain."])
    for ref in evidence.values():
        require(record(ref["path"]) == ref, "Input changed during selection")
    output.mkdir(parents=True)
    write_tsv(output / "bins.tsv", bins, BIN_FIELDS)
    with (output / "selection.json").open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "output"):
        parser.add_argument("--" + name, required=True, type=Path)
    args = parser.parse_args()
    result = run(args.root, args.output)
    print(json.dumps({k: result[k] for k in ("status", "planned_bins", "selected_records", "unique_method_pair_records")}))
