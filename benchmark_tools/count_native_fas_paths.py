"""Count feature paths in retained annotations using the native FAS container."""

import argparse
from copy import deepcopy
import json
from pathlib import Path
import sys
import tempfile
from unittest.mock import patch

from probe_native_fas_omission import fingerprint


def native_options(annotation):
    from greedyFAS import calcFAS, calcFASmulti
    from greedyFAS.mainFAS import greedyFAS

    with tempfile.TemporaryDirectory(prefix="orthohmm-fas-options-") as tmp:
        argv = ["fas.runMultiTaxa", "--input", "unused.txt", "-a", str(annotation.parent),
                "-o", tmp, "--bidirectional", "--tsv", "--domain", "--no_config", "--json",
                "--mergeJson", "--outName", "computed_results", "--max_cardinality", "40",
                "--paths_limit", "15", "--pairLimit", "30000", "--cpus", "4"]
        with patch.object(sys, "argv", argv):
            args = calcFASmulti.get_options()
        args.seed = args.query = annotation.name
        args.out_name = "option_probe"
        args.seed_id = args.query_id = args.ref_proteome = args.ref_2 = args.pairwise = None
        args.silent = True
        toolpath = Path(calcFASmulti.__file__).with_name("pathconfig.txt").read_text().strip()
        captured = []
        with patch.object(greedyFAS, "fc_start", side_effect=lambda o: captured.append(deepcopy(o))):
            calcFAS.fas([args, toolpath])
        if len(captured) != 1:
            raise AssertionError("Native option capture failed")
    keys = ("MS_uni", "input_linearized", "input_normal", "eFeature", "eInstance",
            "max_overlap", "max_overlap_percentage", "paths_limit", "max_cardinality",
            "priority_threshold", "priority_mode")
    return {k: captured[0][k] for k in keys}, fingerprint(Path(toolpath) / "annoTools.txt")


def run(annotation):
    from importlib.metadata import version
    from greedyFAS import calcFAS, calcFASmulti
    from greedyFAS.mainFAS import greedyFAS as engine, fasPathing

    annotation = annotation.resolve()
    source = fingerprint(annotation)
    options, tools = native_options(annotation)
    refs = [fingerprint(p) for p in (__file__, calcFAS.__file__, calcFASmulti.__file__,
                                    engine.__file__, fasPathing.__file__,
                                    Path(__file__).with_name("probe_native_fas_omission.py"))] + [tools, source]
    data = json.loads(annotation.read_text())
    rows = []
    for protein in sorted(data["feature"]):
        linear, features, *_ = engine.su_lin_query_protein(protein, data["feature"], data["clan"], options)
        _, count = engine.pb_region_paths(engine.pb_region_mapper(linear, features,
            options["max_overlap"], options["max_overlap_percentage"]))
        if type(count) is not int or count < 0:
            raise ValueError("Unexpected native path count")
        rejected = options["paths_limit"] > 1 and count > options["paths_limit"]
        rows.append(dict(protein=protein, linearized_instances=len(features), paths=count,
                         rejected_by_path_limit=rejected))
    if refs != [fingerprint(r["path"]) for r in refs]:
        raise AssertionError("Source or annotation changed")
    return dict(status="native_annotation_feature_paths_counted", annotation=source, sources=refs,
        greedyfas_version=version("greedyFAS"), options=options, proteins=rows,
        summary=dict(proteins=len(rows), rejected=sum(r["rejected_by_path_limit"] for r in rows),
                     maximum_paths=max((r["paths"] for r in rows), default=0)),
        historical_omissions_attributed=False,
        limitations=["Native graph counts for one annotation file, not a reconstructed historical FAS sample.",
            "No scoring, endpoint replacement, method comparison or missing-pair attribution is performed.",
            "The native option builder is intercepted before inference solely to capture its effective options."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotation", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    result = run(args.annotation)
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(result["summary"], sort_keys=True))
