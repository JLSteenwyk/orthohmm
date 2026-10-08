"""Prospective native12 history and runtime compatibility, without launching."""

from pathlib import Path
import time

from benchmark_tools.check_threadripper_runtime import RuntimeChecker
from benchmark_tools import native11_composed_review_binding as binding
from benchmark_tools.native_factorial_allocated_execution import reviewed_history
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.review_native11_postterminal_runtime import compare_inventory
from benchmark_tools.run_native_factorial_cost import read
from benchmark_tools.snapshot_runtime_trees import inventory
from benchmark_tools.validate_native_factorial_outputs import require


HISTORY_SCHEMA = "native12_composed_history_v1"
TREE_SCHEMA = "native12_current_runtime_inventory_v1"
LATE_ADDITIONS = [
    dict(path="/usr/bin/lftp", kind="file", mode=493, bytes=1700488,
         sha256="a9ec9b5ebe74f12062260d513ee717d03c3da2c7ea0bceb79eb5a67674f3a249"),
    dict(path="/usr/bin/lftpget", kind="file", mode=493, bytes=1301,
         sha256="51682b3030b25a3dc5b4b9b1c19cb14a63d2350cb18cb4b661df06269392ff05"),
]


def history_gate(request, plan, composed_review, prior_scheduler):
    require(type(request.get("index")) is int and request["index"] == 11
        and request.get("job_id") == 23985
        and isinstance(request.get("history"), list) and len(request["history"]) == 11
        and len(prior_scheduler) == 11
        and [row.get("review") for row in prior_scheduler] == request["history"],
        "Require the actual ordered eleven-identity historical prefix")
    require(len(plan.get("runs", [])) == 13
        and plan["runs"][11].get("index") == 11
        and plan["runs"][11].get("cell") == "p1_c1_r0"
        and type(plan["runs"][12].get("index")) is int
        and plan["runs"][12]["index"] == 12
        and plan["runs"][12].get("cell") == "p1_c1_r1"
        and plan["runs"][12].get("dataset") == "qfo_corrected"
        and composed_review.get("schema") == binding.composed.SCHEMA
        and composed_review.get("status") == "native_success"
        and composed_review.get("index") == 11
        and composed_review.get("job_id") == 23985
        and composed_review.get("review_job_id") == 24034
        and composed_review.get("plan") == request.get("plan")
        and composed_review.get("amendment") == request.get("amendment")
        and all(composed_review.get(k) is True for k in (
            "terminal_reviewed", "composed_full_review_complete", "native_outputs_validated",
            "primary_resources_replayed", "shared_host_resources_reviewed"))
        and all(composed_review.get(k) is False for k in (
            "original_ordinary_full_review_success", "current_original_inventory_equality",
            "next_identity_authorized", "automatic_retry", "publication_ready"))
        and composed_review.get("historical_failures_retained") == [23986, 24033],
        "Require complete new-type native11 review and unchanged next frozen identity")


def history_binding(request_ref, review_ref):
    context = binding.native_binding(request_ref, review_ref)
    request, execution, plan, _, composed_review = context[:5]
    prior_scheduler = reviewed_history(request, request["amendment"], execution, plan)
    history_gate(request, plan, composed_review, prior_scheduler)
    refs = [*request["history"], *context[8], record(__file__), record(binding.__file__)]
    for ref in refs:
        check(ref)
    return context, dict(schema=HISTORY_SCHEMA, status="prefix_revalidated_for_native12_preparation",
        original_request=request_ref, composed_review=review_ref,
        historical_prefix=request["history"], historical_scheduler=prior_scheduler,
        native11_scheduler=context[7], composed_binding=context[9], next_unrun_index=12,
        original_review_translated=False, original_review_authorization_changed=False,
        historical_failures_retained=[23986, 24033], automatic_retry=False,
        execution_request_prepared=False, production_execution_launched=False,
        next_identity_authorized=False, scientific_timings_admitted=False,
        publication_ready=False, evidence=refs)


def require_unrun(run, plan, request_path):
    require(type(run.get("index")) is int and run["index"] == 12
        and run.get("cell") == "p1_c1_r1", "Only the final frozen native identity is supported")
    paths = [Path(run["output_root"]), Path(plan["panel_root"]) / "sessions/run_12",
             Path(request_path)]
    require(all(path.is_absolute() and path.resolve() == path for path in paths),
        "Require direct absolute unrun paths")
    require(all(not path.exists() and not path.is_symlink() for path in paths),
        "Native12 already has output/session/request; no retry or resume")
    return paths


def runtime_gate(report, specifications):
    require(report.get("schema") == "native11_composed_runtime_gate_v1"
        and report.get("status") == "historical_brackets_and_current_exact_additions_revalidated"
        and report.get("current_original_inventory_equality") is False
        and report.get("current_original_entries_unchanged") is True
        and report.get("current_private_inventory_equality") is True
        and report.get("continuous_runtime_integrity_established") is False
        and report.get("current_additions") == LATE_ADDITIONS,
        "Require exact known additions, not a general runtime-mutation exemption")
    rows = report.get("current_inventories")
    require(isinstance(rows, list) and len(rows) == len(specifications) == 2,
        "Require the two original runtime trees in order")
    for number, (row, (path, digest)) in enumerate(zip(rows, specifications)):
        require(row.get("manifest", {}).get("path") == str(path)
            and row["manifest"].get("sha256") == digest,
            "Prospective runtime tree binding differs")
        delta = row.get("comparison", {})
        require(delta.get("missing") == [] and delta.get("changed") == []
            and delta.get("metadata_differences") == {}
            and delta.get("added") == (LATE_ADDITIONS if number == 0 else [])
            and delta.get("equal") is (number == 1)
            and type(delta.get("expected_records")) is int and delta["expected_records"] > 0
            and type(delta.get("observed_records")) is int
            and delta["observed_records"] == delta["expected_records"] + (2 if number == 0 else 0),
            "Retained runtime comparison does not establish the exact current identity")


def current_trees(specifications, runtime_ref):
    report = read(runtime_ref)
    runtime_gate(report, specifications)
    source = record(binding.composed.__file__)
    require(source["sha256"] == binding.COMPOSED_SOURCE_SHA and report.get("source") == source
        and report.get("runtime_component") == record(binding.composed.RUNTIME),
        "Composed runtime source or component differs")
    results = []
    for number, (path, digest) in enumerate(specifications):
        manifest_ref = record(path)
        require(manifest_ref["sha256"] == digest
            and manifest_ref == report["current_inventories"][number]["manifest"],
            "Original runtime inventory bytes changed")
        manifest = read(manifest_ref)
        comparison = compare_inventory(manifest, inventory(manifest["roots"]))
        require(comparison == report["current_inventories"][number]["comparison"],
            "Current runtime differs from the explicitly bound prospective identity")
        check(manifest_ref)
        results.append(dict(schema=TREE_SCHEMA, manifest=manifest_ref,
            status="prospective_current_runtime_identity_matches", comparison=comparison,
            original_inventory_equality=number == 1,
            prospective_inventory_equality=True, scientific_execution_authorized=False,
            continuous_runtime_integrity_established=False))
    check(runtime_ref)
    check(source)
    check(report["runtime_component"])
    return results


class ProspectiveRuntimeChecker:
    """Fresh exact tree checks plus the unchanged private-interpreter lookup probe."""

    def __init__(self, plan, runtime_ref, output):
        self.runtime_ref = runtime_ref
        self.output = Path(output)
        self.source = record(__file__)
        lookup = read(plan["runtime_lookup"])
        self.specifications = read(lookup["binding"])["runtime_specs"]
        require(lookup["baseline"] == plan["baseline"], "Lookup baseline differs")
        runtime_gate(read(runtime_ref), self.specifications)
        self.checker = RuntimeChecker(Path(plan["runtime_lookup"]["path"]),
            plan["runtime_lookup"]["sha256"], self.output, tree_checker=self.check_trees)

    @property
    def count(self):
        return self.checker.count

    def check_trees(self, specifications):
        require([[str(path), digest] for path, digest in specifications] == self.specifications,
            "Fresh tree specifications differ")
        started = time.monotonic()
        result = current_trees(specifications, self.runtime_ref)
        check(self.source)
        save(self.output / f"tree_check_{self.count:02d}.json", dict(
            schema="native12_current_runtime_check_v1", source=self.source,
            runtime_basis=self.runtime_ref, inventories=result,
            check_wall_s=time.monotonic() - started, observed_unix_ns=time.time_ns(),
            original_os_inventory_equality=False, current_identity_revalidated=True,
            next_identity_authorized=False, continuous_runtime_integrity_established=False))
        return result

    def __call__(self, specifications):
        check(self.source)
        result = self.checker(specifications)
        check(self.source)
        return result
