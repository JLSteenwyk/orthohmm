import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmark_tools.audit_candidate_arm import audit
from benchmark_tools.prepare_ob_candidate_neighborhood import ARMS, record
from benchmark_tools.prepare_qfo_candidate_neighborhood import (
    ENVIRONMENT, PLAN_SHA, PROTOCOL_SHA, prepare_arms, require_environment,
)


def fixture(tmp_path, fail_variant=False):
    params = {"min_norm": .03, "min_margin": 1.5}
    seed = tmp_path / "seed.txt"
    seed.write_text("a b\nc\n")
    seen = []
    def merge(*args, **kwargs):
        seen.append(kwargs)
        if fail_variant and kwargs["min_norm"] == .024:
            raise RuntimeError("test native failure")
    engine = SimpleNamespace(merge_supported_satellite_candidate_clusters=merge)
    def expand(directory, names, species, hits, profile):
        engine.merge_supported_satellite_candidate_clusters(**params, merge_trace=[])
        w = Path(directory) / "orthohmm_working_res"
        (w / "phylogeny_candidate_superfamilies.txt").write_bytes((w / "orthohmm_edges_clustered.txt").read_bytes())
        (w / "phylogeny_candidate_merges.json").write_text("[]\n")
        (w / "phylogeny_candidate_seeds.tsv").write_text(
            "candidate_family\tseed_families\nFamily0000000\tSeed0000000\nFamily0000001\tSeed0000001\n")
        return {"parameters": dict(params), "profile": profile, "membership_policy": "high_confidence_pair",
                "candidate_checkpoint": str(w / "phylogeny_candidate_superfamilies.txt"),
                "seed_sidecar": str(w / "phylogeny_candidate_seeds.tsv"),
                "merge_trace_sidecar": str(w / "phylogeny_candidate_merges.json"),
                "seed_families": 2, "candidate_families": 2, "merges": 0}
    engine._expand_phylogeny_candidates = expand
    trace = tmp_path / "baseline_trace.json"
    trace.write_text("[]\n")
    baseline = {"seed_partition": record(seed), "candidate_partition": record(seed),
                "membership_constraints": record(trace), "expansion": {"parameters": params}}
    out = tmp_path / "output"
    out.mkdir()
    return engine, merge, seen, baseline, out, {"arms": []}


def test_all_candidate_arms_and_real_content_audit(tmp_path):
    engine, original, seen, baseline, output, report = fixture(tmp_path)
    prepare_arms(engine, baseline, list("abc"), [0, 1, 2], (), output, report, audit)
    assert engine.merge_supported_satellite_candidate_clusters is original
    assert [r["label"] for r in report["arms"]] == [k for k, _ in ARMS]
    assert len(seen) == 5
    assert report["arms"][0]["baseline_byte_equivalent"] is True
    for row, (label, delta) in zip(report["arms"], ARMS):
        assert row["status"] == "candidate_prepared_unscored"
        assert row["applied_parameters"] == {**baseline["expansion"]["parameters"], **delta}
        assert row["candidate_arm"]["content_audit"]["genes"] == 3
        assert row["candidate_arm"]["content_audit"]["accuracy_evaluated"] is False
    assert json.loads((output / "progress.json").read_text()) == report


@pytest.mark.parametrize("key", ["candidate_partition", "membership_constraints"])
def test_control_difference_prevents_variants(tmp_path, key):
    engine, original, seen, baseline, output, report = fixture(tmp_path)
    baseline[key] = {**baseline[key], "sha256": "changed"}
    with pytest.raises(ValueError, match="Unchanged corrected control"):
        prepare_arms(engine, baseline, list("abc"), [0, 1, 2], (), output, report, audit)
    assert len(seen) == 1 and not (output / "norm_low").exists()
    assert engine.merge_supported_satellite_candidate_clusters is original


def test_variant_failure_preserves_control_and_restores_engine(tmp_path):
    engine, original, seen, baseline, output, report = fixture(tmp_path, fail_variant=True)
    with pytest.raises(RuntimeError, match="test native failure"):
        prepare_arms(engine, baseline, list("abc"), [0, 1, 2], (), output, report, audit)
    assert len(seen) == 2
    assert engine.merge_supported_satellite_candidate_clusters is original
    assert report["arms"][0]["status"] == "candidate_prepared_unscored"
    assert report["arms"][1]["status"] == "preparing"
    assert not (output / "norm_high").exists()
    assert (output / "control/orthohmm_working_res/orthohmm_edges_clustered.txt").is_file()


@pytest.mark.parametrize("key", [None, "SLURM_JOB_ID", "SLURM_CPUS_PER_TASK", "SLURM_JOB_NODELIST", *ENVIRONMENT])
def test_environment(key):
    env = {**ENVIRONMENT, "SLURM_JOB_ID": "123", "SLURM_CPUS_PER_TASK": "2", "SLURM_JOB_NODELIST": "bizon"}
    if key is None:
        require_environment(env)
    else:
        env[key] = ""
        with pytest.raises(ValueError):
            require_environment(env)


def test_actual_protocol_and_plan_pins():
    root = Path(__file__).resolve().parents[2] / "benchmark_tools/results"
    for name, digest in [("qfo_parameter_neighborhood_plan_20260919.json", PLAN_SHA),
                         ("QFO_PARAMETER_NEIGHBORHOOD_PROTOCOL_20260919.md", PROTOCOL_SHA)]:
        assert hashlib.sha256((root / name).read_bytes()).hexdigest() == digest
    plan = json.loads((root / "qfo_parameter_neighborhood_plan_20260919.json").read_text())
    selected = [r for r in plan["arms"] if not r["label"].startswith("cpm_")]
    assert selected == [{"label": label, "delta": delta} for label, delta in ARMS]
