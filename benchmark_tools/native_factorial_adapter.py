"""Explicit P/C/R ablations of the frozen pipeline, without editing its package.

P-off still uses high-sensitivity initial HMM search and multipass grouping.
Only C-on/R-off needs an AST change, removing the reconciliation gate from
candidate expansion. Namespace copies keep production globals unchanged.
This library does not authorize or measure full-dataset execution.
"""

import ast
from dataclasses import asdict, replace
import hashlib
import inspect
from pathlib import Path
import re
from types import FunctionType


CORE_SHA = "2afb89b9dc683e64e58208188f720e07701d4760ff09d1ac3a53c7a8075b84bb"
GATE = ast.parse('phylogeny == "reconcile" and phylogeny_candidates in {"satellite_v1", "satellite_v2"}',
                 mode="eval").body


def factors(cell):
    if not isinstance(cell, str) or not re.fullmatch(r"p[01]_c[01]_r[01]", cell):
        raise ValueError("Require an explicit factorial cell")
    return dict(cell=cell, profile_expansion=cell[1] == "1",
                candidate_expansion=cell[4] == "1", reconciliation=cell[7] == "1")


def execution_tree(source, cell):
    flags = factors(cell)
    if hashlib.sha256(source).hexdigest() != CORE_SHA:
        raise ValueError("Require exact frozen pipeline source bytes")
    tree = ast.parse(source)
    functions = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "_execute"]
    if len(functions) != 1:
        raise ValueError("Frozen execution function missing or ambiguous")
    function = functions[0]
    matching = [n for n in ast.walk(function) if isinstance(n, ast.If)
                and ast.dump(n.test) == ast.dump(GATE)]
    if len(matching) != 1:
        raise ValueError("Candidate gate missing or ambiguous")
    if flags["candidate_expansion"] and not flags["reconciliation"]:
        matching[0].test = matching[0].test.values[1]
    return ast.fix_missing_locations(ast.Module(body=[function], type_ignores=[]))


def entrypoint(module, cell):
    flags = factors(cell)
    path = Path(module.__file__).resolve()
    source = path.read_bytes()
    tree = execution_tree(source, cell)
    original = module.resolve_accuracy_profile("high_sensitivity")
    active = replace(original, profile_expansion=flags["profile_expansion"])
    if (original.kmer_k, original.max_candidates_per_query, original.multipass_graph,
            original.profile_expansion, original.leiden_seed) != (4, 100, True, True, 4):
        raise ValueError("Frozen high-sensitivity profile differs")
    changes = {k for k in asdict(original) if asdict(original)[k] != asdict(active)[k]}
    if changes != (set() if flags["profile_expansion"] else {"profile_expansion"}):
        raise ValueError("Unexpected accuracy-profile change")
    namespace = dict(module.__dict__)
    def resolve(name):
        if name != "high_sensitivity":
            raise ValueError("Factorial ablation cannot substitute standard search")
        return active
    namespace["resolve_accuracy_profile"] = resolve
    exec(compile(tree, str(path), "exec"), namespace)
    compiled = namespace["_execute"]
    if inspect.signature(compiled) != inspect.signature(module._execute):
        raise ValueError("Frozen execution signature differs")
    report = dict(flags, frozen_pipeline_sha256=CORE_SHA,
        active_accuracy_profile=asdict(active), original_accuracy_profile=asdict(original),
        candidate_gate_decoupled=flags["candidate_expansion"] and not flags["reconciliation"],
        production_source_modified=False, production_globals_modified=False)
    def execute(**kwargs):
        expected = dict(accuracy_profile="high_sensitivity", search_mode="builtin", clustering="leiden",
            cpm_resolution=.1, refinement_profile="default", evalue_threshold=.0001,
            start=None, phylogeny="reconcile" if flags["reconciliation"] else "off",
            phylogeny_candidates="satellite_v2" if flags["candidate_expansion"] else "seed",
            species_tree_mode="infer", species_tree_rooting="min_variance",
            phylogeny_root_rule="species_overlap", phylogeny_pair_rule="positive_paralogy")
        if any(kwargs.get(k) != v for k, v in expected.items()):
            raise ValueError("Arguments differ from the explicit frozen factorial configuration")
        kwargs["metrics"].add_metadata(native_factorial=report)
        return compiled(**kwargs)
    namespace["_execute"] = execute
    original_entry = module.execute
    adapted = FunctionType(original_entry.__code__, namespace, original_entry.__name__,
                           original_entry.__defaults__, original_entry.__closure__)
    adapted.__kwdefaults__ = original_entry.__kwdefaults__
    return adapted, report
