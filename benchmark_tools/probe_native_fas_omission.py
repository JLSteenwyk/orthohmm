"""Controlled branch probe for the retained FAS container, not historical attribution."""

import argparse
from contextlib import ExitStack
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import tempfile
from unittest.mock import patch


def fingerprint(path):
    path = Path(path).resolve()
    data = path.read_bytes()
    return dict(path=str(path), bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def probe():
    from importlib.metadata import version
    from greedyFAS.mainFAS import greedyFAS as engine, fasOutput
    from greedyFAS import calcFAS, calcFASmulti

    path = shutil.which("fas_benchmark.py")
    if path is None:
        raise RuntimeError("Run inside the retained QfO FAS container")
    # The native scorer imports helpers beside its executable.
    import sys
    sys.path.insert(0, str(Path(path).parent))
    spec = importlib.util.spec_from_file_location("native_qfo_fas", path)
    scorer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(scorer)
    refs = [fingerprint(p) for p in (__file__, path, engine.__file__, fasOutput.__file__,
                                    calcFAS.__file__, calcFASmulti.__file__)]
    limit = 10 ** 15
    options = dict(MS_uni=1, max_overlap=0, max_overlap_percentage=0.4,
        paths_limit=limit, priority_threshold=30, max_cardinality=40,
        priority_mode=True, domain=False, pairwise=[("P1", "P2")], progress=False)
    rows = []
    # Inject only the graph preparation result. Native cutoff, NA construction,
    # serialization and loader execute unchanged; no real annotation is evaluated.
    with ExitStack() as stack:
        stack.enter_context(patch.object(engine, "su_lin_query_protein", return_value=([], {}, {}, {}, {})))
        stack.enter_context(patch.object(engine, "pb_region_mapper", return_value={}))
        for paths, cutoff, rejected in ((limit - 1, limit, False), (limit, limit, False),
                                        (limit + 1, limit, True), (limit + 1, 1, False)):
            with patch.object(engine, "pb_region_paths", return_value=({}, paths)):
                value = engine.fc_prep_query("P2", {}, {}, dict(options, paths_limit=cutoff), {})
                if (value is None) != rejected:
                    raise AssertionError("Native path-limit boundary differs")
                rows.append(dict(injected_paths=paths, configured_limit=cutoff, rejected=rejected))
        with patch.object(engine, "pb_region_paths", return_value=({}, limit + 1)):
            forward = engine.fc_main({}, {"P1": {}}, {"P2": {}}, {}, options, {}, {})
            reverse = engine.fc_main({}, {"P2": {}}, {"P1": {}}, {},
                                     dict(options, pairwise=[("P2", "P1")]), {}, {})
        if forward != [("P1", "P2", ("NA",) * 5, "NA")]:
            raise AssertionError("Native rejected pair did not produce NA")
    synthetic = []
    native_options = dict(options, input_linearized=["pfam"], input_normal=[],
                          eFeature=0.001, eInstance=0.01)
    for regions, rejected in ((49, False), (50, True)):
        annotation = dict(length=5000, pfam={f"pfam_{i}": dict(evalue=0,
            instance=[[100 * j + 1, 100 * j + 50, 0] for j in range(regions)]) for i in range(2)})
        proteome = {"P1": annotation, "P2": annotation}
        linear, features, *_ = engine.su_lin_query_protein("P1", proteome, {}, native_options)
        _, count = engine.pb_region_paths(engine.pb_region_mapper(linear, features, 0, 0.4))
        prepared = engine.fc_prep_query("P1", {}, proteome, native_options, {})
        if count != 2 ** regions or (prepared is None) != rejected:
            raise AssertionError("Unmocked synthetic feature graph disagrees")
        synthetic.append(dict(regions=regions, alternatives_per_region=2,
                              native_path_count=count, rejected=rejected))
        if rejected:
            actual = engine.fc_main({}, proteome, proteome, {}, native_options, {}, {})
            if actual != forward:
                raise AssertionError("Unmocked synthetic graph did not produce NA")
    with tempfile.TemporaryDirectory(prefix="orthohmm-fas-probe-") as tmp:
        output = str(Path(tmp) / "native")
        fasOutput.write_json_out(output, True, (forward, reverse))
        raw = json.loads(Path(output + ".json").read_text())
        loaded = scorer.load_precomputed_fas_scores(Path(output + ".json"))
        if raw != {"P1_P2": ["NA", "NA"]} or loaded != {}:
            raise AssertionError("Native NA serialization or omission differs")
        controls = {"P1_P2": ["0.2", "0.8"], "P3_P4": ["NA", "0.8"],
                    "P5_P6": ["0.2", "NA"]}
        control_file = Path(tmp) / "controls.json"
        control_file.write_text(json.dumps(controls))
        control_loaded = scorer.load_precomputed_fas_scores(control_file)
        if control_loaded != {("P1", "P2"): 0.5}:
            raise AssertionError("Numeric or mixed-NA control differs")
    if refs != [fingerprint(r["path"]) for r in refs]:
        raise AssertionError("Probe source changed")
    return dict(status="native_fas_path_limit_omission_mechanism_verified",
        greedyfas_version=version("greedyFAS"), sources=refs, boundary_cases=rows,
        unmocked_synthetic_architectures=synthetic,
        rejected_pair_json=raw, rejected_pair_loader_entries=len(loaded),
        control_json=controls, control_retained=[dict(pair=list(k), value=v) for k, v in control_loaded.items()],
        historical_omissions_attributed=False, benchmark_scores_changed=False,
        limitations=["Boundary cases mock graph preparation; separate synthetic architecture cases use native graph construction.",
            "Synthetic architecture checks do not calculate any real protein's feature-path count.",
            "Confirms a native mechanism, not the identities or causes of the historical missing scores.",
            "No historical RNG state, sampled missing.txt or temporary FAS JSON is reconstructed."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    result = probe()
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
