"""Gate all eight frozen cells before evaluating any OrthoBench outcome."""

import argparse
import csv
import io
import json
import math
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import SeqIO
from benchmark_tools.recover_factorial_postflight import PREPARED_HASH
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.run_orthobench_factorial_cell import select_cell, verify_prepared
from benchmark_tools.validate_factorial_native import validate
from benchmark_tools.score_ygob_groups import membership, read_predictions
from benchmark_tools.score_orthobench_partition import read_named_groups, score_partition
from benchmark_tools.bootstrap_orthobench_factorial import CELLS, factorial_bootstrap, render_report
from benchmark_tools.orthobench_stage_diagnostics import file_provenance, run_official_benchmark

REFERENCE_SNAPSHOT = "660ead29c5b6ac0b8278cd1e62cdcdb0a513db81dda317e634b805e661d70ba9"
OFFICIAL_HASH = "81eb1e660c17819549b07eea8a54b4fb42a89180cafeb4569d92195d282f5e6f"


def coverage_and_resources(groups, gene_species, metrics=None):
    seen = membership({str(i): list(group) for i, group in enumerate(groups)})
    if set(seen) != set(gene_species):
        raise ValueError("Coverage requires the complete input gene universe")
    sizes = [len(group) for group in groups]
    multispecies = [group for group in groups if len({gene_species[g] for g in group}) > 1]
    result = {"input_genes": len(gene_species), "assigned_genes": len(seen),
              "groups": len(groups), "singleton_groups": sizes.count(1),
              "genes_in_nonsingleton_groups": sum(n for n in sizes if n > 1),
              "multispecies_groups": len(multispecies),
              "genes_in_multispecies_groups": sum(map(len, multispecies)),
              "resource_scope": "Upstream cached partition; separate cell cost not measured",
              "wall_s": None, "user_cpu_s": None, "system_cpu_s": None,
              "mean_cpu_cores": None, "peak_process_tree_rss_bytes": None}
    if metrics is not None:
        if metrics.get("rss_measurement") != "sampled_sum_of_linux_proc_tree_rss":
            raise ValueError("Unexpected memory accounting convention")
        for key in ("wall_s", "user_cpu_s", "system_cpu_s", "peak_process_tree_rss_bytes"):
            value = metrics[key]
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
                raise ValueError("Invalid resource measurement")
            result[key] = value
        if result["wall_s"] <= 0:
            raise ValueError("Nonpositive wall time")
        result["mean_cpu_cores"] = (result["user_cpu_s"] + result["system_cpu_s"]) / result["wall_s"]
        result["resource_scope"] = "Incremental cached reconciliation; shared node; sampled sum of process-tree RSS"
    return result


def render_execution(result):
    lines = ["", "## Coverage And Incremental Resources", "",
             "All assigned genes includes singletons; multispecies coverage is not orthology accuracy.",
             "Times are shared-node cached reconciliation costs, not end-to-end efficiency comparisons.",
             "NA for candidate-only cells means no separate cell cost was measured, not zero cost.", "",
             "| Cell | Groups | Singleton Groups | Genes In Multispecies Groups | Wall (s) | Mean CPU Cores | Peak Tree RSS (GiB) |",
             "| --- | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for label in CELLS:
        row = result["coverage_resources"][label]
        timing = "NA | NA | NA" if row["wall_s"] is None else f"{row['wall_s']:.2f} | {row['mean_cpu_cores']:.2f} | {row['peak_process_tree_rss_bytes'] / 2**30:.3f}"
        lines.append(f"| {label} | {row['groups']} | {row['singleton_groups']} | {row['genes_in_multispecies_groups']} | {timing} |")
    lines += ["", "Original reconciliation batch failures are retained in the machine-readable scheduler records.",
              "Recovered postflight verification errors are not relabeled as native inference failures."]
    return "\n".join(lines) + "\n"


def require_terminal_array(accounting):
    rows = list(csv.DictReader(io.StringIO(accounting), delimiter="|"))
    expected = {f"21248_{i}" for i in range(4)}
    tasks = [r for r in rows if r["JobID"] in expected]
    if len(tasks) != 4 or {r["JobID"] for r in tasks} != expected:
        raise ValueError("Missing or duplicate factorial scheduler tasks")
    if any(r["State"] not in {"COMPLETED", "FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY", "NODE_FAIL"} for r in tasks):
        raise ValueError("Factorial still running; no partial scoring")
    return tasks


def compare_official(score, official):
    for key in ("f_score", "precision", "recall"):
        # The frozen official CLI prints one decimal place in percentage units.
        if key not in official or not math.isfinite(official[key]) or abs(score[key] - official[key]) > 0.05000001:
            raise ValueError(f"Official scorer disagreement: {key}")
    if official.get("exact_refogs") != score["exact_refogs"]:
        raise ValueError("Official exact-family count disagreement")


def load_reference_snapshot(path):
    snapshot = read_frozen(path, REFERENCE_SNAPSHOT)
    inputs = snapshot["inputs"]
    records = inputs["references"] + inputs["uncertain"]
    for record in records:
        if file_provenance(Path(record["path"])) != record:
            raise ValueError("Reference resource changed")
    directory = Path(inputs["references"][0]["path"]).parent
    if {p.resolve() for p in directory.glob("RefOG*.txt")} != {Path(r["path"]).resolve() for r in inputs["references"]}:
        raise ValueError("Reference file set changed")
    if {p.resolve() for p in (directory / "low_certainty_assignments").glob("RefOG*.txt")} != {Path(r["path"]).resolve() for r in inputs["uncertain"]}:
        raise ValueError("Low-certainty file set changed")
    names = sorted(Path(r["path"]).name for r in inputs["references"])
    if names != [f"RefOG{i:03d}.txt" for i in range(1, 71)]:
        raise ValueError("Expected full 70-RefOG panel")
    references = read_named_groups(directory, names)
    uncertain = {Path(r["path"]).name: set(Path(r["path"]).read_text().split()) for r in inputs["uncertain"]}
    official = directory.parent / "benchmark.py"
    if file_provenance(official)["sha256"] != OFFICIAL_HASH:
        raise ValueError("Official scorer source changed")
    return references, uncertain, official, records


def assemble(root, output):
    accounting = subprocess.check_output(["sacct", "-j", "21248", "--parsable2",
                                          "--format=JobID,State,ExitCode"], text=True)
    tasks = require_terminal_array(accounting)
    if output.exists():
        raise FileExistsError(output)
    results = root / "benchmark_tools/results"
    prepared = read_frozen(results / "orthobench_factorial_prepared_20260916.json", PREPARED_HASH)
    # A genuine failure or integrity error halts assembly for explicit diagnosis;
    # it cannot silently remove a cell or turn failed scheduling into zero accuracy.
    native = {f"p{i // 2}_c{i % 2}_r1": validate(root, i) for i in range(4)}
    cell, _, launcher = select_cell(prepared, 0)
    verify_prepared(prepared, cell, launcher, 21161)
    genes = set()
    gene_species = {}
    for item in prepared["fasta_inputs"]:
        for record in SeqIO.parse(item["path"], "fasta"):
            if record.id in genes:
                raise ValueError("Duplicate FASTA ID")
            genes.add(record.id)
            gene_species[record.id] = Path(item["path"]).name
    predictions, prediction_sources, coverage = {}, {}, {}
    for cell in prepared["cells"]:
        path = Path(cell["prediction"])
        if cell["reconciliation"]:
            groups = read_predictions(path, "root_hogs")
        else:
            with path.open() as handle:
                groups = {str(i): line.split() for i, line in enumerate(handle) if line.strip()}
            arm = prepared["candidate_arms"][f"p{int(cell['profile_expansion'])}_c{int(cell['candidate_expansion'])}"]
            if file_provenance(path) != arm["candidate_partition"]:
                raise ValueError("Candidate checkpoint changed")
        if set(membership(groups)) != genes:
            raise ValueError("Cell does not partition the complete input universe")
        predictions[cell["label"]] = [set(g) for g in groups.values()]
        prediction_sources[cell["label"]] = file_provenance(path)
        metrics = None
        if cell["reconciliation"]:
            metrics = json.loads(Path(native[cell["label"]]["native_metrics"]["path"]).read_text())
        coverage[cell["label"]] = coverage_and_resources(predictions[cell["label"]], gene_species, metrics)
    if set(predictions) != set(CELLS):
        raise ValueError("Incomplete factorial design")
    references, uncertain, official, reference_records = load_reference_snapshot(results / "orthobench_paired_uncertainty_20260916.json")
    if not set().union(*references.values()).issubset(genes):
        raise ValueError("Reference genes absent from inference universe")
    output.mkdir(parents=True)
    scores, official_scores = {}, {}
    for label in CELLS:
        converted = output / f"{label}.txt"
        with converted.open("x") as handle:
            for group in predictions[label]:
                handle.write(" ".join(sorted(group)) + "\n")
        scores[label] = score_partition(predictions[label], references, uncertain)
        official_scores[label] = run_official_benchmark(official, converted)
        compare_official(scores[label], official_scores[label])
    for record in [*reference_records, *prediction_sources.values()]:
        if file_provenance(Path(record["path"])) != record:
            raise ValueError("Scoring input changed during evaluation")
    result = factorial_bootstrap({label: {"status": "complete", "refog_records": scores[label]["refog_records"]} for label in CELLS})
    result.update(schema_version=1, scores=scores, official_scores=official_scores,
                  native_validation=native, scheduler=tasks, accounting_raw=accounting,
                  predictions=prediction_sources, references=reference_records, coverage_resources=coverage,
                  official_scorer=file_provenance(official), assembler=file_provenance(Path(__file__)),
                  publication_ready=False, timing_scope="Incremental cached reconciliation on a shared node")
    with (output / "results.json").open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    with (output / "results.md").open("x") as handle:
        handle.write(render_report(result) + render_execution(result))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    assemble(args.root.resolve(), args.output.resolve())
    print("Eight cells scored with official crosschecks and prespecified paired uncertainty")


if __name__ == "__main__":
    main()
