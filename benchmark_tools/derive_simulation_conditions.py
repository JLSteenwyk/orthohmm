"""Export frozen missingness/sampling transforms from validated native truth."""

import argparse
import json
from pathlib import Path
import sys

from Bio import Phylo

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.benchmark_production import file_record
from benchmark_tools.simulation_conditions import project_truth, select_removals
from benchmark_tools.zombi_truth import validate_run


def derive(run, output, seed):
    if output.exists():
        raise FileExistsError("Refusing to overwrite transformed inputs")
    truth, sequences = validate_run(run)
    tree_path = run / "T/ExtantTree.nwk"
    tree = Phylo.read(tree_path, "newick")
    owners = {gene: species for gene, (species, _) in sequences.items()}
    selections = {condition: select_removals(owners, tree, seed, condition)
                  for condition in ("missing20", "uneven_taxa", "taxon_count_control")}
    output.mkdir(parents=True)
    report = {"schema_version": 1, "seed": seed, "accuracy_computed": False,
              "baseline_run": str(run.resolve()), "baseline_truth": truth,
              "baseline_species_tree": file_record(tree_path, run),
              "sources": [file_record(p, p.parent) for p in
                  (Path(__file__).resolve(), Path(__file__).with_name("simulation_conditions.py").resolve(),
                   Path(__file__).with_name("zombi_truth.py").resolve())],
              "conditions": {}}
    for condition, selection in selections.items():
        if selection["status"] != "selected":
            report["conditions"][condition] = selection
            continue
        retained, pairs = project_truth(truth["ortholog_pairs"], owners, selection)
        target = output / condition
        inputs = target / "input"
        inputs.mkdir(parents=True)
        for species in selection["retained_species"]:
            with (inputs / f"{species}.fasta").open("w") as handle:
                for gene in sorted(retained):
                    if retained[gene] == species:
                        handle.write(f">{gene}\n{sequences[gene][1]}\n")
        projected = {"schema_version": 1, "condition": condition, "selection": selection,
                     "scope": truth["scope"], "extant_genes": len(retained),
                     "species": selection["retained_species"], "ortholog_pairs": pairs,
                     "ortholog_pair_count": len(pairs),
                     "families": {family: [g for g in genes if g in retained]
                                  for family, genes in truth["families"].items()},
                     "inputs": [file_record(p, target) for p in sorted(inputs.glob("*.fasta"))],
                     "species_tree_note": "No supplied tree exported; primary methods infer trees from retained inputs"}
        path = target / "truth.json"
        path.write_text(json.dumps(projected, indent=2, sort_keys=True) + "\n")
        report["conditions"][condition] = {"status": "complete", "selection": selection,
            "genes": len(retained), "eligible_true_pairs": len(pairs), "truth": file_record(path, output)}
    (output / "manifest.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-run", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, required=True)
    args = parser.parse_args()
    report = derive(args.baseline_run, args.output, args.seed)
    print(json.dumps({c: {k: r[k] for k in ("status", "genes", "eligible_true_pairs") if k in r}
                      for c, r in report["conditions"].items()}, sort_keys=True))


if __name__ == "__main__":
    main()
