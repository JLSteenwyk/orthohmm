"""Prospective fixture-schema correction; reuse the frozen native context kernels."""

import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import sys
import tempfile

import probe_native_fas_context as base


BASE_SHA = "d722944a7cadeee3e25433fc472cdf724ba58309aaac29952338e5c7575bae06"
TOOL_MAPPING = {name: name.lower() for name in base.TOOLS}


def fixtures():
    annotations, owners = base.fixtures()
    for document in annotations.values():
        for gene, features in document["feature"].items():
            document["feature"][gene] = {TOOL_MAPPING.get(name, name): value for name, value in features.items()}
    return annotations, owners


def probe(protocol):
    from importlib.metadata import version
    from greedyFAS import calcFAS, calcFASmulti
    from greedyFAS.mainFAS import greedyFAS as engine, fasInput, fasOutput, fasPathing, fasScoring, fasWeighting

    require, fingerprint = base.require, base.fingerprint
    protocol_ref = fingerprint(protocol)
    require(protocol_ref["sha256"] == base.PROTOCOL_SHA, "Changed prospective protocol")
    executable = shutil.which("fas_benchmark.py")
    require(executable is not None and version("greedyFAS") == "1.18.7", "Use retained native container")
    sys.path.insert(0, str(Path(executable).parent))
    spec = importlib.util.spec_from_file_location("native_fas_context_scorer", executable)
    scorer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(scorer)
    paths = dict(scorer=executable, calcFAS=calcFAS.__file__, calcFASmulti=calcFASmulti.__file__, engine=engine.__file__)
    refs = {name: fingerprint(path) for name, path in paths.items()}
    require(all(refs[name]["sha256"] == digest for name, digest in base.PINS.items()), "Changed historical native source")
    for module in (fasInput, fasOutput, fasPathing, fasScoring, fasWeighting):
        refs[module.__name__] = fingerprint(module.__file__)
    refs["driver_original"] = fingerprint(base.__file__)
    require(refs["driver_original"]["sha256"] == BASE_SHA, "Changed frozen context kernels")
    refs["driver"] = fingerprint(__file__)
    tool_config = Path(calcFAS.__file__).parent / "annoTools.txt"
    refs["native_tool_config"] = fingerprint(tool_config)
    linear, normal = fasInput.featuretypes(str(tool_config))
    require(linear == [TOOL_MAPPING[name] for name in base.TOOLS[:2]] and
            normal == [TOOL_MAPPING[name] for name in base.TOOLS[2:]], "Changed native tool configuration")
    annotations, owners = fixtures()
    fixture_sha = hashlib.sha256(json.dumps(annotations, sort_keys=True).encode()).hexdigest()
    options = dict(MS_uni=1, input_linearized=linear, input_normal=normal,
                   eFeature=.001, eInstance=.01, max_overlap=0, max_overlap_percentage=.4)
    counts = {}
    for gene in ("C", "F", "G"):
        feature = annotations[owners[gene]]["feature"]
        linear_features, features, *_ = engine.su_lin_query_protein(gene, feature, {}, options)
        _, counts[gene] = engine.pb_region_paths(engine.pb_region_mapper(linear_features, features, 0, .4))
    require(counts == {"C": 2 ** 16, "F": 2 ** 16, "G": 2 ** 50}, "Changed native architecture contexts")
    rows = []
    with tempfile.TemporaryDirectory(prefix="native-fas-context-lowercase-") as directory:
        annotation_dir = Path(directory) / "annotations"
        annotation_dir.mkdir()
        annotation_refs = []
        for tax, document in annotations.items():
            path = annotation_dir / (tax + ".json")
            path.write_text(json.dumps(document, sort_keys=True))
            annotation_refs.append(fingerprint(path))
        for context in base.contexts():
            row = base.execute_context(scorer, context, annotation_dir, owners)
            rows.append(row)
            print(context["name"], row["terminal"]["returncode"], flush=True)
            if row["error"] is not None or row["terminal"]["returncode"] != 0:
                break
        require(all(fingerprint(ref["path"]) == ref for ref in annotation_refs), "Annotation fixtures mutated")
    assessment_error = None
    try:
        assessment = base.compare(rows)
    except ValueError as exc:
        assessment_error = str(exc)
        assessment = dict(planned_contexts=len(base.contexts()), executed_contexts=len(rows), complete=False,
                          all_tested_contexts_invariant=False)
    require(all(fingerprint(ref["path"]) == ref for ref in refs.values()), "Native source changed")
    require(fingerprint(protocol) == protocol_ref, "Protocol changed")
    return dict(schema="native_fas_context_probe_lowercase_v2", status="controlled_native_context_invariance_passed"
                if assessment["all_tested_contexts_invariant"] else "controlled_native_context_invariance_not_established",
                protocol=protocol_ref, sources=refs, greedyfas_version=version("greedyFAS"),
                fixture_annotations=annotations, fixture_sha256=fixture_sha, owners=owners,
                tool_key_normalization=TOOL_MAPPING, native_path_counts=counts,
                contexts=rows, assessment=assessment, assessment_error=assessment_error,
                historical_scores_rerun=False, historical_omissions_attributed=False, native_sampling_law_admitted=False,
                native_intervals_admitted=False, publication_ready=False,
                limitations=["Controlled annotations and successful contexts only; finite checks are not proof for all native pairs.",
                             "Fixed numeric return/status is a design assumption; batch crashes are outside this check.",
                             "No historical random state, sample, missing identities, benchmark score or annotation is reconstructed.",
                             "No confidence construction or independent biological generalization is admitted."])


def main():
    original_probe = base.probe
    try:
        base.probe = probe
        return base.main()
    finally:
        base.probe = original_probe


if __name__ == "__main__":
    raise SystemExit(main())
