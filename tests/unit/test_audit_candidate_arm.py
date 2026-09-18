import json

import pytest

from benchmark_tools.audit_candidate_arm import audit
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def fixture(tmp_path, expanded):
    seed = tmp_path / "seed.txt"
    seed.write_text("a b\nc\nd\n")
    directory = tmp_path / "arm"
    working = directory / "orthohmm_working_res"
    working.mkdir(parents=True)
    partition = working / "orthohmm_edges_clustered.txt"
    partition.write_text("a b c\nd\n" if expanded else seed.read_text())
    arm = {"seed_partition": record(seed), "candidate_expansion": expanded}
    if expanded:
        (working / "phylogeny_candidate_superfamilies.txt").write_bytes(partition.read_bytes())
        (working / "phylogeny_candidate_seeds.tsv").write_text(
            "candidate_family\tseed_families\nFamily0000000\tSeed0000000,Seed0000001\nFamily0000001\tSeed0000002\n")
        trace = [{"source_genes": ["c"], "target_genes": ["a", "b"], "source_size": 1,
                  "target_size": 2, "source_seed_families": 1, "target_seed_families": 1}]
        (working / "phylogeny_candidate_merges.json").write_text(json.dumps(trace))
        arm["expansion"] = {"profile": "satellite_v2", "membership_policy": "high_confidence_pair",
            "candidate_checkpoint": str(working / "phylogeny_candidate_superfamilies.txt"),
            "seed_sidecar": str(working / "phylogeny_candidate_seeds.tsv"),
            "merge_trace_sidecar": str(working / "phylogeny_candidate_merges.json"),
            "seed_families": 3, "candidate_families": 2, "merges": 1}
    refresh(arm, directory)
    return arm, record(seed), directory


def refresh(arm, directory):
    working = directory / "orthohmm_working_res"
    arm["candidate_partition"] = record(working / "orthohmm_edges_clustered.txt")
    if arm["candidate_expansion"]:
        arm["membership_constraints"] = record(working / "phylogeny_candidate_merges.json")
    arm["output_files"] = [record(p) for p in sorted(directory.rglob("*")) if p.is_file()]


@pytest.mark.parametrize("expanded", [False, True])
def test_valid_content(tmp_path, expanded):
    args = fixture(tmp_path, expanded)
    result = audit(*args, set("abcd"), expanded)
    assert result["genes"] == 4 and result["seed_families"] == 3
    assert result["candidate_families"] == (2 if expanded else 3)
    assert result["merges"] == int(expanded)
    assert result["accuracy_evaluated"] is False


@pytest.mark.parametrize("problem", ["seed", "inventory", "hash", "missing", "duplicate", "split", "superfamilies",
    "sidecar", "sidecar_path", "profile", "trace_partial", "trace_duplicate", "trace_missing", "trace_foreign",
    "trace_size", "trace_redundant", "summary", "unexpected_file"])
def test_bad_expanded_content(tmp_path, problem):
    arm, seed, directory = fixture(tmp_path, True)
    work = directory / "orthohmm_working_res"
    partition = work / "orthohmm_edges_clustered.txt"
    trace_path = work / "phylogeny_candidate_merges.json"
    trace = json.loads(trace_path.read_text())
    if problem == "seed":
        arm["seed_partition"] = {}
    elif problem == "inventory":
        arm["output_files"].pop()
    elif problem == "hash":
        partition.write_text("a b\nc d\n")
    elif problem in ("missing", "duplicate", "split"):
        partition.write_text({"missing": "a b c\n", "duplicate": "a b c\nc d\n", "split": "a c\nb d\n"}[problem])
    elif problem == "superfamilies":
        (work / "phylogeny_candidate_superfamilies.txt").write_text("changed\n")
    elif problem == "sidecar":
        (work / "phylogeny_candidate_seeds.tsv").write_text("candidate_family\tseed_families\n")
    elif problem == "sidecar_path":
        arm["expansion"]["seed_sidecar"] = "foreign"
    elif problem == "profile":
        arm["expansion"]["profile"] = "satellite_v1"
    elif problem == "trace_partial":
        trace[0]["target_genes"] = ["a"]
    elif problem == "trace_duplicate":
        trace[0]["source_genes"] = ["c", "c"]
    elif problem == "trace_missing":
        trace = []
    elif problem == "trace_foreign":
        trace[0]["source_genes"] = ["z"]
    elif problem == "trace_size":
        trace[0]["target_size"] = 3
    elif problem == "trace_redundant":
        trace.append(dict(trace[0]))
    elif problem == "summary":
        arm["expansion"]["merges"] = 2
    elif problem == "unexpected_file":
        (work / "extra.txt").write_text("unexpected")
    if problem.startswith("trace_"):
        trace_path.write_text(json.dumps(trace))
    if problem not in ("inventory", "hash"):
        refresh(arm, directory)
    with pytest.raises(ValueError):
        audit(arm, seed, directory, set("abcd"), True)


def test_off_arm_cannot_change_even_to_an_equivalent_partition(tmp_path):
    arm, seed, directory = fixture(tmp_path, False)
    (directory / "orthohmm_working_res/orthohmm_edges_clustered.txt").write_text("b a\nc\nd\n")
    refresh(arm, directory)
    with pytest.raises(ValueError, match="exact seed"):
        audit(arm, seed, directory, set("abcd"), False)


@pytest.mark.parametrize("next_iteration", [False, True])
def test_multiple_attachments_and_prior_iteration_components(tmp_path, next_iteration):
    arm, seed, directory = fixture(tmp_path, True)
    work = directory / "orthohmm_working_res"
    for name in ("orthohmm_edges_clustered.txt", "phylogeny_candidate_superfamilies.txt"):
        (work / name).write_text("a b c d\n")
    (work / "phylogeny_candidate_seeds.tsv").write_text(
        "candidate_family\tseed_families\nFamily0000000\tSeed0000000,Seed0000001,Seed0000002\n")
    path = work / "phylogeny_candidate_merges.json"
    trace = json.loads(path.read_text())
    trace.append({"source_genes": ["d"], "target_genes": ["a", "b", "c"] if next_iteration else ["a", "b"],
                  "source_size": 1, "target_size": 3 if next_iteration else 2,
                  "source_seed_families": 1, "target_seed_families": 2 if next_iteration else 1})
    path.write_text(json.dumps(trace))
    arm["expansion"].update(candidate_families=1, merges=2)
    refresh(arm, directory)
    assert audit(arm, seed, directory, set("abcd"), True)["merges"] == 2
    if next_iteration:
        trace.reverse()
        path.write_text(json.dumps(trace))
        refresh(arm, directory)
        with pytest.raises(ValueError, match="preceding trace"):
            audit(arm, seed, directory, set("abcd"), True)
