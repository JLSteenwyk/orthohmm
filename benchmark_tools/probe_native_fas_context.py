"""Check native FAS pair outcomes across controlled batch contexts, not benchmarks."""

import argparse
from copy import deepcopy
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
from unittest.mock import patch


PROTOCOL_SHA = "e47b22d456dea7a41fcdf9248981d0c2d4eafec5c6dcb6ace5ad2153297b70dd"
PINS = {
    "scorer": "1045c57f4d0f4787bec3d1f0690799df63c68dcc67ccfd5a925c338e3d33661d",
    "calcFASmulti": "7f9621c0dc6737c8b928c0cf35575101decc4af1cf039e4f3c044e80ce460b75",
    "calcFAS": "40a920fcce24d76e7441f5f9368aa6232d779815956742c725ca16a169980c2b",
    "engine": "f82a78d2787d30f106a58c760cae49d28f28d55956ad188b0b259ec8b9f0f7e8",
}
FOCAL = (("A", "D"), ("B", "E"), ("C", "F"), ("G", "H"))
COMPANIONS = (("I", "J"), ("K", "L"))
TOOLS = ("Pfam", "SMART", "fLPS", "COILS2", "SEG", "SignalP", "TMHMM")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def fingerprint(path):
    path = Path(path).resolve(strict=True)
    data = path.read_bytes()
    return dict(path=str(path), bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def architecture(features, length=6000):
    result = dict(length=length, **{tool: {} for tool in TOOLS})
    for name, instances in features.items():
        result["Pfam"][name] = dict(evalue=0, instance=[[start, stop, 0] for start, stop in instances])
    return result


def fixtures():
    simple = architecture({"pfam_F1": [(10, 60)]})
    partial = architecture({"pfam_F1": [(10, 60)], "pfam_F2": [(110, 160)]})
    priority = {name: [(110 * i + 10, 110 * i + 60) for i in range(16)] for name in ("pfam_F1", "pfam_F2")}
    rejected = {name: [(110 * i + 10, 110 * i + 60) for i in range(50)] for name in ("pfam_F1", "pfam_F2")}
    features = dict(A=simple, D=deepcopy(simple), B=partial, E=deepcopy(simple),
                    C=architecture(priority), F=architecture(priority, length=6400),
                    G=architecture(rejected), H=deepcopy(simple),
                    I=deepcopy(simple), J=deepcopy(simple), K=deepcopy(partial), L=deepcopy(simple))
    owners = dict.fromkeys(("A", "B", "C", "G", "I"), "T1")
    owners.update(dict.fromkeys(("D", "E", "F", "H"), "T2"))
    owners.update(J="T3", K="T4", L="T3")
    annotations = {tax: dict(feature={gene: deepcopy(value) for gene, value in features.items() if owners[gene] == tax},
                            count={}, clan={}, version={"context_fixture": 1}) for tax in sorted(set(owners.values()))}
    return annotations, owners


def contexts():
    full = FOCAL + COMPANIONS
    rows = [dict(name="singleton_" + str(i), pairs=[pair], workers=1, hash_seed="0") for i, pair in enumerate(FOCAL)]
    rows.extend((dict(name="full_forward_one", pairs=list(full), workers=1, hash_seed="0"),
                 dict(name="full_reverse_one", pairs=list(reversed(full)), workers=1, hash_seed="0"),
                 dict(name="full_rotated_two", pairs=list(full[2:] + full[:2]), workers=2, hash_seed="0"),
                 dict(name="full_reverse_two_hash1", pairs=list(reversed(full)), workers=2, hash_seed="1")))
    rows.extend(dict(name="subset_" + str(i), pairs=[COMPANIONS[0], pair, COMPANIONS[1]], workers=2, hash_seed="0")
                for i, pair in enumerate(FOCAL))
    return rows


def canonical_raw(raw, pairs):
    expected = {"_".join(pair): tuple(pair) for pair in pairs}
    require(set(raw) == set(expected), "Incomplete or unexpected native JSON pairs")
    normalized = {}
    for name, directions in raw.items():
        require(isinstance(directions, list) and len(directions) == 2, "Invalid directional scores")
        if any(value == "NA" for value in directions):
            value = None
        else:
            numbers = [float(value) for value in directions]
            require(all(math.isfinite(value) and 0 <= value <= 1 for value in numbers), "Invalid native score")
            value = sum(numbers) / 2
        normalized[name] = dict(pair=list(expected[name]), directions=directions, returned=value is not None, value=value)
    return normalized


def compare(rows):
    plan = contexts()
    require(len(rows) <= len(plan), "Unexpected contexts")
    baseline, variations, failures = {}, [], []
    evaluations = 0
    for row, context in zip(rows, plan):
        require(row["context"] == {**context, "pairs": [list(p) for p in context["pairs"]]}, "Changed context plan")
        if row["error"] is not None or row["terminal"]["returncode"] != 0:
            failures.append(context["name"])
            continue
        normalized = canonical_raw(row["raw"], context["pairs"])
        require(row["normalized"] == normalized, "Changed raw/normalized binding")
        evaluations += len(normalized)
        for name, value in normalized.items():
            pair = tuple(value["pair"])
            if pair == FOCAL[0]:
                require(value["value"] == 1, "Identical numeric control failed")
            elif pair == FOCAL[1]:
                require(value["returned"] and 0 < value["value"] < 1, "Nonidentical numeric control failed")
            elif pair == FOCAL[2]:
                require(value["returned"], "Priority numeric control failed")
            elif pair == FOCAL[3]:
                require(not value["returned"] and value["directions"] == ["NA", "NA"], "NA control failed")
            if name not in baseline:
                baseline[name] = value
            elif value != baseline[name]:
                variations.append(dict(context=context["name"], pair=value["pair"], baseline=baseline[name], observed=value))
    complete = len(rows) == len(plan) and not failures
    return dict(planned_contexts=len(plan), executed_contexts=len(rows), complete=complete,
                pair_evaluations=evaluations, failures=failures, variations=variations,
                all_tested_contexts_invariant=complete and not variations)


def run_child(command, hash_seed):
    environment = dict(os.environ, PYTHONHASHSEED=hash_seed, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    child = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                             env=environment, start_new_session=True)
    timed_out = False
    try:
        stdout, stderr = child.communicate(timeout=120)
    except subprocess.TimeoutExpired:
        timed_out = True
        os.killpg(child.pid, signal.SIGKILL)
        stdout, stderr = child.communicate()
    return subprocess.CompletedProcess(command, child.returncode, stdout, stderr), timed_out


def execute_context(scorer, context, annotations, owners):
    captured = []

    def invalid_constant(value):
        raise ValueError("Nonfinite JSON constant: " + value)

    def native_subprocess(command, **kwargs):
        require(not kwargs and not captured, "Changed native subprocess contract")
        argv = [str(value) for value in command]
        record = dict(command=argv, returncode=None, timed_out=False,
                      stdout=None, stderr=None, raw=None, raw_text=None)
        captured.append(record)
        terminal, timed_out = run_child(argv, context["hash_seed"])
        record.update(returncode=terminal.returncode, timed_out=timed_out,
                      stdout=terminal.stdout, stderr=terminal.stderr)
        output = Path(argv[argv.index("-o") + 1]) / "computed_results.json"
        if terminal.returncode == 0 and output.exists():
            record["raw_text"] = output.read_text()
            record["raw"] = json.loads(record["raw_text"], parse_constant=invalid_constant)
        return terminal

    error = None
    normalized = None
    stage = "native_helper_and_io"
    try:
        # Only capture the I/O boundary; native helper, CLI, score functions and loader stay unchanged.
        with patch.object(scorer.subprocess, "run", native_subprocess):
            lookup = scorer.compute_fas_scores_for_pairs(list(context["pairs"]), owners, annotations, context["workers"])
        stage = "native_output_validation"
        require(len(captured) == 1, "Native child was not captured")
        if captured[0]["returncode"] == 0:
            require(captured[0]["raw"] is not None, "Successful child lacks native JSON")
            normalized = canonical_raw(captured[0]["raw"], context["pairs"])
            expected = {tuple(value["pair"]): value["value"] for value in normalized.values() if value["returned"]}
            require(lookup == expected, "Native loader differs from directional JSON")
    except Exception as exc:
        error = dict(stage=stage, type=type(exc).__name__, message=str(exc))
    terminal = captured[0] if captured else dict(command=None, returncode=None, timed_out=False,
                                               stdout=None, stderr=None, raw=None, raw_text=None)
    raw = terminal.pop("raw")
    raw_text = terminal.pop("raw_text")
    return dict(context={**context, "pairs": [list(pair) for pair in context["pairs"]]},
                terminal=terminal, raw=raw, raw_text=raw_text, normalized=normalized, error=error)


def probe(protocol):
    from importlib.metadata import version
    from greedyFAS import calcFAS, calcFASmulti
    from greedyFAS.mainFAS import greedyFAS as engine, fasInput, fasOutput, fasPathing, fasScoring, fasWeighting
    import shutil

    protocol_ref = fingerprint(protocol)
    require(protocol_ref["sha256"] == PROTOCOL_SHA, "Changed prospective protocol")
    executable = shutil.which("fas_benchmark.py")
    require(executable is not None and version("greedyFAS") == "1.18.7", "Use retained native container")
    sys.path.insert(0, str(Path(executable).parent))
    spec = importlib.util.spec_from_file_location("native_fas_context_scorer", executable)
    scorer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(scorer)
    paths = dict(scorer=executable, calcFAS=calcFAS.__file__, calcFASmulti=calcFASmulti.__file__, engine=engine.__file__)
    refs = {name: fingerprint(path) for name, path in paths.items()}
    require(all(refs[name]["sha256"] == digest for name, digest in PINS.items()), "Changed historical native source")
    for module in (fasInput, fasOutput, fasPathing, fasScoring, fasWeighting):
        refs[module.__name__] = fingerprint(module.__file__)
    refs["driver"] = fingerprint(__file__)
    annotations, owners = fixtures()
    fixture_sha = hashlib.sha256(json.dumps(annotations, sort_keys=True).encode()).hexdigest()
    options = dict(MS_uni=1, input_linearized=["Pfam", "SMART"], input_normal=list(TOOLS[2:]),
                   eFeature=.001, eInstance=.01, max_overlap=0, max_overlap_percentage=.4)
    counts = {}
    for gene in ("C", "F", "G"):
        feature = annotations[owners[gene]]["feature"]
        linear, features, *_ = engine.su_lin_query_protein(gene, feature, {}, options)
        _, count = engine.pb_region_paths(engine.pb_region_mapper(linear, features, 0, .4))
        counts[gene] = count
    require(counts == {"C": 2 ** 16, "F": 2 ** 16, "G": 2 ** 50}, "Changed native architecture contexts")
    rows = []
    with tempfile.TemporaryDirectory(prefix="native-fas-context-") as directory:
        annotation_dir = Path(directory) / "annotations"
        annotation_dir.mkdir()
        annotation_refs = []
        for tax, document in annotations.items():
            path = annotation_dir / (tax + ".json")
            path.write_text(json.dumps(document, sort_keys=True))
            annotation_refs.append(fingerprint(path))
        for context in contexts():
            row = execute_context(scorer, context, annotation_dir, owners)
            rows.append(row)
            print(context["name"], row["terminal"]["returncode"], flush=True)
            if row["error"] is not None or row["terminal"]["returncode"] != 0:
                break
        require(all(fingerprint(ref["path"]) == ref for ref in annotation_refs), "Annotation fixtures mutated")
    assessment_error = None
    try:
        assessment = compare(rows)
    except ValueError as exc:
        assessment_error = str(exc)
        assessment = dict(planned_contexts=len(contexts()), executed_contexts=len(rows), complete=False,
                          all_tested_contexts_invariant=False)
    require(all(fingerprint(ref["path"]) == ref for ref in refs.values()), "Native source changed")
    require(fingerprint(protocol) == protocol_ref, "Protocol changed")
    return dict(schema="native_fas_context_probe_v1", status="controlled_native_context_invariance_passed"
                if assessment["all_tested_contexts_invariant"] else "controlled_native_context_invariance_not_established",
                protocol=protocol_ref, sources=refs, greedyfas_version=version("greedyFAS"),
                fixture_annotations=annotations, fixture_sha256=fixture_sha, owners=owners, native_path_counts=counts,
                contexts=rows, assessment=assessment, assessment_error=assessment_error,
                historical_scores_rerun=False, historical_omissions_attributed=False, native_sampling_law_admitted=False,
                native_intervals_admitted=False, publication_ready=False,
                limitations=["Controlled annotations and successful contexts only; finite checks are not proof for all native pairs.",
                             "Fixed numeric return/status is a design assumption; batch crashes are outside this check.",
                             "No historical random state, sample, missing identities, benchmark score or annotation is reconstructed.",
                             "No confidence construction or independent biological generalization is admitted."])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    allowed = sorted(os.sched_getaffinity(0))
    require(len(allowed) >= 2, "Two CPUs required for prospective worker check")
    os.sched_setaffinity(0, allowed[:2])
    import resource
    resource.setrlimit(resource.RLIMIT_AS, (4 * 1024 ** 3, 4 * 1024 ** 3))
    result = probe(args.protocol)
    result["resources"] = dict(cpu_affinity=allowed[:2], per_process_address_space_bytes=4 * 1024 ** 3,
                               maximum_child_workers=2, child_timeout_seconds=120, numerical_threads=1)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"status": result["status"], "assessment": result["assessment"]}))
    return 0 if result["assessment"]["all_tested_contexts_invariant"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
