"""Read-only Git/producer/source/selection audit, not a new inference experiment."""

import collections
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys


PRIMARY = {"OrthoBench: OrthoBench F1", "QfO: VGNC F", "QfO: SwissTrees F", "QfO: TreeFam-A F"}
ENDPOINTS = PRIMARY | {"QfO: EC", "QfO: GO", "QfO: FAS", "ThreeKingdoms: Three Kingdoms F1"}
NATIVE = {"hmm_viterbi.so": "f21d902852d727d6270b65ba900db52e7216cd874ee2d2b89fa72b6f16dcd5c0",
    "kmer_prefilter.so": "348b1c715e2e305dae7e889bae1d2498e609e48ce5d50cd668db71b548c60b88",
    "pair_align.so": "81c26c6ce92b1da37285824ea1b8dcea007a7fc18a11ee5fa3d150628fa800c7"}


def git(*args):
    return subprocess.check_output(["git", *args])


def pin(path):
    path = Path(path).resolve()
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1048576), b""):
            digest.update(block)
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": digest.hexdigest()}


def require_pin(item):
    assert pin(item["path"]) == item, item["path"]


def tree(revision):
    entries = {}
    for entry in git("ls-tree", "-r", "-z", revision).split(b"\0"):
        if entry:
            metadata, path = entry.split(b"\t", 1)
            entries[os.fsdecode(path)] = metadata.split()[2].decode()
    return entries


def blob_id(path):
    content = Path(path).read_bytes()
    return hashlib.sha1(b"blob " + str(len(content)).encode() + b"\0" + content).hexdigest()


def main(root):
    destination = root / "protocol_readback.json"
    if destination.exists():
        raise FileExistsError("Preserve prior protocol audits")
    specification = root / "protocol_specification.json"
    spec = json.loads(specification.read_text())
    summary_path = Path(spec["comparison"])
    summary = json.loads(summary_path.read_text())
    methods = {m["key"]: m for m in summary["methods"]}
    head = git("rev-parse", "HEAD").decode().strip()
    assert git("rev-parse", "--show-object-format").strip() == b"sha1"
    head_tree = tree(head)
    checked = {str(p.resolve()): pin(p) for p in (specification, summary_path, Path(__file__))}
    inventory = []
    changed = git("diff", "--name-only", spec["baseline_revision"] + ".." + head, "--", "benchmark_tools").decode().splitlines()
    for name in sorted(changed):
        path = Path(name)
        if (path.parent != Path("benchmark_tools") or path.suffix != ".json"
                or not path.name.startswith(("hillclimb_", "confirmation_"))):
            continue
        item = pin(path)
        assert blob_id(path) == head_tree[name], name
        checked[item["path"]] = item
        first = git("log", "--diff-filter=A", "--format=%H", "--", name).decode().splitlines()
        assert len(first) == 1
        inventory.append(dict(source=item, first_committed_revision=first[0], git_blob=head_tree[name],
                              protocol=json.loads(path.read_text())))
    source_trees = {}
    rounds = []
    for round_spec in spec["rounds"]:
        name = round_spec["round"]
        screening = summary["round_screening"][name]
        control = methods[round_spec["control"]]
        plans = []
        for planned in round_spec["plans"]:
            revision = git("rev-parse", planned["revision"]).decode().strip()
            data = Path(planned["path"]).read_bytes()
            assert git("show", revision + ":" + planned["path"]) == data, planned
            first = next(p["first_committed_revision"] for p in inventory if p["source"]["path"] == str(Path(planned["path"]).resolve()))
            assert subprocess.run(["git", "merge-base", "--is-ancestor", first, revision]).returncode == 0
            plans.append(dict(path=planned["path"], frozen_execution_revision=revision,
                first_committed_revision=first, exact_plan_present_in_execution_revision=True))
        observations = []
        vectors = {}
        for label, declared in screening["comparisons"].items():
            key = round_spec["control"] if label == screening["control"] else label if ":" in label else name + ":" + label
            method = methods[key]
            deltas = {endpoint: value - control["scores"][endpoint] for endpoint, value in method["scores"].items()
                      if endpoint in control["scores"]}
            assert deltas == declared["deltas_from_matched_control"], (name, label)
            complete = set(deltas) == ENDPOINTS
            passes = complete and any(deltas.get(k, 0) > 0 for k in PRIMARY) and all(v >= -.005 for v in deltas.values())
            assert complete == declared["complete_eight_endpoint_contrast"]
            assert passes == declared["passes_point_estimate_screen"], (name, label)
            if set(method["scores"]) == ENDPOINTS:
                vectors[label] = method["scores"]
            observations.append(dict(label=label, method=key, deltas=deltas, complete=complete, passes=passes,
                endpoints_declining_more_than_screen=[k for k, v in deltas.items() if v < -.005],
                endpoint_count=len(method["scores"])))
        domination = {label: sorted(other for other, values in vectors.items() if other != label
            and all(values[k] >= vector[k] for k in ENDPOINTS) and any(values[k] > vector[k] for k in ENDPOINTS))
            for label, vector in vectors.items()}
        assert domination == screening["pareto"]["dominated_by"]
        assert sorted(k for k, v in domination.items() if not v) == screening["pareto"]["nondominated_complete_labels"]
        rounds.append(dict(round=name, plans=plans, comparisons=observations,
            pareto_dominators=domination, interpretation=round_spec["gate_interpretation"]))
    producers, execution_plans, source_checks = {}, {}, {}
    frozen_code_checks = set()
    for method in summary["methods"][8:]:
        for cell in method["cells"]:
            if cell["status"] != "independent_content_readback_verified":
                continue
            assert cell["coverage"]["genes"] == {"OrthoBench": 251378, "QfO": 984137, "ThreeKingdoms": 443217}[cell["dataset"]]
            assert cell["coverage"]["species"] == {"OrthoBench": 12, "QfO": 78, "ThreeKingdoms": 12}[cell["dataset"]]
            item = cell.get("source_snapshot") or cell.get("independent_readback") or cell.get("independent_admission")
            if not isinstance(item, dict):
                # Canonical OrthoBench/Three Kingdoms cells have direct checked evidence instead of a producer snapshot.
                assert method["key"] == "canonical_collection:canonical"
                continue
            if item["path"] in producers:
                continue
            require_pin(item)
            checked[item["path"]] = item
            producer = json.loads(Path(item["path"]).read_text())
            plan_pin = producer.get("execution_plan")
            producers[item["path"]] = dict(source=item, status=producer.get("status"), execution_plan=plan_pin)
            if not isinstance(plan_pin, dict) or plan_pin["path"] in execution_plans:
                continue
            require_pin(plan_pin)
            checked[plan_pin["path"]] = plan_pin
            plan = json.loads(Path(plan_pin["path"]).read_text())
            revision, source_root = plan["source_revision"], Path(plan["source_root"])
            source_trees.setdefault(revision, tree(revision))
            code_count = native_count = 0
            for source in plan["sources"]:
                if source["path"] not in source_checks:
                    require_pin(source)
                    source_checks[source["path"]] = source
                else:
                    assert source_checks[source["path"]] == source
                path = Path(source["path"])
                try:
                    relative = path.relative_to(source_root).as_posix()
                except ValueError:
                    continue
                if path.suffix in (".py", ".c"):
                    assert relative in source_trees[revision], (revision, relative)
                    assert blob_id(path) == source_trees[revision][relative], (revision, relative)
                    frozen_code_checks.add((revision, relative, source["sha256"]))
                    code_count += 1
                elif path.name in NATIVE:
                    assert source["sha256"] == NATIVE[path.name]
                    native_count += 1
            execution_plans[plan_pin["path"]] = dict(source=plan_pin, revision=revision, source_root=str(source_root),
                checked_source_records=len(plan["sources"]), git_verified_code_records=code_count,
                original_native_libraries_verified=native_count)
    pilot = methods["protected_union:protected_union"]["cells"][0]
    baseline = next(c for c in methods["canonical_collection:canonical"]["cells"] if c["dataset"] == "OrthoBench")
    strict_gate = (pilot["scores"]["OrthoBench: OrthoBench F1"] > baseline["scores"]["OrthoBench: OrthoBench F1"]
        and all(pilot["details"][k] - baseline["details"][k] >= -.005 for k in ("precision", "recall")))
    assert strict_gate is False and len(methods["protected_union:protected_union"]["scores"]) == 1
    for item in [*checked.values(), *source_checks.values()]:
        require_pin(item)
    result = dict(status="core_round_protocol_git_source_and_screening_readback_verified",
        audited_head=head, specification=pin(specification), recipe=pin(__file__),
        protocol_inventory=inventory, rounds=rounds, producers=list(producers.values()),
        execution_plans=list(execution_plans.values()), checked=list(checked.values()),
        unique_frozen_source_records_checked=len(source_checks), unique_git_code_identities_checked=len(frozen_code_checks),
        protected_union_strict_gate_passed=False,
        limitations=["All core-round plan bytes are present at their declared executed driver/scientific revisions. This is Git ancestry, not historical push/execution attestation.",
            "All ten complete/partial round screening and point-Pareto results were independently recomputed with the original 0.005 tolerance. The stricter union pilot gate remains failed.",
            "Producer status strings retain historical snapshots, including running; they are not current scheduler liveness claims. Actual cell admission and retained terminal accounting take precedence.",
            "This is a source/selection audit, not new raw scoring. Full per-gate execution evidence and every supporting diagnostic plan still require the final requirement audit.",
            "Initial audit-specification lookup guessed the QfO native recovery driver as 402d3232; it lacked that later plan. The pinned actual execution manifest identifies 6b64c6fa, used here without changing any historical source."],
        remaining=["Complete supporting-diagnostic and per-gate execution-evidence review against the original eight requirements.",
                   "Issue final separate-endpoint default recommendation and next justified experiments, then perform completion audit."],
        scientific_inference_rerun=False, defaults_promoted=False, goal_completed=False)
    with destination.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"protocols": len(inventory), "rounds": len(rounds), "producers": len(producers),
        "execution_plans": len(execution_plans), "source_records": len(source_checks), "git_code_identities": len(frozen_code_checks)}))


if __name__ == "__main__":
    main(Path(sys.argv[1]))
