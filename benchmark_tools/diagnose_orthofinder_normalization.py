"""Reproduce native score normalization on all divergent variable-panel seeds.

Run with the frozen OrthoFinder interpreter. No tool files, inputs, predictions
or scientific scores are changed. Includes failed and successful seeds.
"""

import argparse
from contextlib import redirect_stdout
import importlib.metadata
import io
import json
from pathlib import Path
import sys
from types import SimpleNamespace
import warnings

from Bio import SeqIO
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.validate_simulation_outputs import verify_process, validate_graph_weights, NativeOutputFailure


def species_lengths(path, species):
    lengths = {}
    for record in SeqIO.parse(path, "fasta"):
        parts = record.id.split("_")
        if len(parts) != 2 or not all(p.isdigit() for p in parts) or int(parts[0]) != species:
            raise ValueError("Unexpected native species/sequence identifier")
        index = int(parts[1])
        if index in lengths or not record.seq:
            raise ValueError("Duplicate or empty native sequence")
        lengths[index] = len(record.seq)
    if not lengths or sorted(lengths) != list(range(len(lengths))):
        raise ValueError("Native sequence indices must be contiguous")
    return np.array([lengths[i] for i in range(len(lengths))], dtype=float)


def design_summary(products):
    values = np.asarray(products, dtype=float)
    if not np.isfinite(values).all() or (values <= 0).any():
        raise ValueError("Invalid length products")
    design = np.column_stack((np.log10(values), np.ones(len(values))))
    return {"observations": len(values), "unique_length_products": np.unique(values).tolist(),
            "design_rank": int(np.linalg.matrix_rank(design)) if len(values) else 0}


def diagnose_pair(matrix, lengths, i, j, waterfall):
    li, lj, scores = waterfall.scnorm.GetLengthArraysForMatrix(matrix, lengths[i], lengths[j])
    top_l, top_s = waterfall.scnorm.GetTopPercentileOfScores(li * lj, scores, 95)
    result = {"query_species": i, "target_species": j, "raw_hits": matrix.nnz,
              "all_hit_design": design_summary(li * lj), "selected_fit_design": design_summary(top_l)}
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        if len(top_s) > 1:
            try:
                params = waterfall.scnorm.CalculateFittingParameters(top_l, top_s)
                result["fit_parameters"] = [float(p) if np.isfinite(p) else None for p in params]
            except Exception as error:
                result["fit_error"] = f"{type(error).__name__}: {error}"
    result["fit_warnings"] = [str(w.message) for w in caught]
    with warnings.catch_warnings(record=True) as caught, redirect_stdout(io.StringIO()) as stdout:
        warnings.simplefilter("always")
        try:
            normalized = waterfall.WaterfallMethod.NormaliseScores(matrix, lengths, i, j).tocsr()
            result.update(normalized_stored=normalized.nnz,
                          nonfinite_normalized=int(np.count_nonzero(~np.isfinite(normalized.data))))
        except Exception as error:
            result["normalization_error"] = f"{type(error).__name__}: {error}"
    result["normalization_warnings"] = [str(w.message) for w in caught]
    result["normalization_stdout"] = stdout.getvalue()
    return result


def diagnose_dataset(dataset, waterfall, blast_reader):
    config = dataset["methods"]["orthofinder_full"]
    status_path = Path(dataset["methods"]["orthohmm_high_sensitivity"]["output"]).parent / "execution/status.json"
    status = json.loads(status_path.read_text())
    inventory = verify_process(config, status["methods"]["orthofinder_full"])
    roots = list(Path(config["output"]).glob("**/WorkingDirectory"))
    if len(roots) != 1:
        raise ValueError("Expected one native working directory")
    root = roots[0]
    fastas = {i: root / f"Species{i}.fa" for i in range(8)}
    lengths = {i: species_lengths(p, i) for i, p in fastas.items()}
    if not {p.resolve() for p in fastas.values()}.issubset(inventory):
        raise ValueError("Unverified native FASTA input")
    info = SimpleNamespace(nSeqsPerSpecies={i: len(values) for i, values in lengths.items()})
    report = {"dataset": dataset["label"], "seed": dataset["seed"], "pairs": [],
              "execution_evidence": file_provenance(status_path),
              "fastas": [file_provenance(p) for p in fastas.values()],
              "unique_proteome_lengths": {str(i): len(np.unique(values)) for i, values in lengths.items()}}
    for i in range(8):
        for j in range(8):
            path = root / f"Blast{i}_{j}.txt.gz"
            if path.resolve() not in inventory or path.with_suffix("").exists():
                raise ValueError("Missing or ambiguous inventoried native BLAST file")
            matrix = blast_reader.GetBLAST6Scores(info, [str(root)], i, j,
                                                  qExcludeSelfHits=True, qDoubleBlast=True, q_allow_empty=False)
            pair = diagnose_pair(matrix, lengths, i, j, waterfall)
            pair["input"] = file_provenance(path)
            report["pairs"].append(pair)
    graphs = list(Path(config["output"]).glob("**/OrthoFinder_graph.txt"))
    if len(graphs) != 1 or graphs[0].resolve() not in inventory:
        raise ValueError("Expected one verified native graph")
    report["graph"] = file_provenance(graphs[0])
    try:
        validate_graph_weights(graphs[0])
        report["native_graph_finite"] = True
    except NativeOutputFailure as error:
        report.update(native_graph_finite=False, native_graph_error=str(error))
    report["nonfinite_pairs"] = [[p["query_species"], p["target_species"]] for p in report["pairs"] if p.get("nonfinite_normalized", 0)]
    verify_process(config, status["methods"]["orthofinder_full"])
    return report


def main():
    from orthofinder.tools import waterfall
    from orthofinder.utils import blast_file_processor
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--manifest-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("Refusing to overwrite diagnostic evidence")
    manifest = read_frozen(args.manifest, args.manifest_sha256)
    if importlib.metadata.version("orthofinder") != "3.1.5":
        raise ValueError("Expected native OrthoFinder 3.1.5")
    names = sorted({d.metadata["Name"] for d in importlib.metadata.distributions() if d.metadata["Name"]})
    environment = {"python": sys.version, "packages": {n: importlib.metadata.version(n) for n in names}}
    if environment != manifest["environments"]["orthofinder"]:
        raise ValueError("OrthoFinder interpreter/package inventory changed")
    expected = {r["absolute_path"]: r for r in manifest["orthofinder_distribution"]}
    native_sources = []
    for module in (waterfall, blast_file_processor):
        path = Path(module.__file__).resolve()
        actual = file_provenance(path)
        if actual["sha256"] != expected[str(path)]["sha256"]:
            raise ValueError("Native normalization source differs from frozen installation")
        native_sources.append(actual)
    datasets = [d for d in manifest["datasets"] if d["condition"] == "divergent"]
    if len(datasets) != 10 or {d["seed"] for d in datasets} != set(range(20261101, 20261111)):
        raise ValueError("Require all ten divergent variable-length seeds")
    result = {"schema_version": 1, "diagnostic_only": True, "accuracy_computed": False,
              "source": file_provenance(Path(__file__)), "native_sources": native_sources,
              "manifest": file_provenance(args.manifest), "command": [sys.executable, *sys.argv],
              "python": sys.version, "packages": {n: importlib.metadata.version(n) for n in ("orthofinder", "numpy", "scipy", "biopython")},
              "datasets": [diagnose_dataset(d, waterfall, blast_file_processor) for d in datasets]}
    with args.output.open("x") as handle:
        handle.write(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("Diagnosed all ten seeds and 640 ordered native score matrices; no scientific scores changed")


if __name__ == "__main__":
    main()
