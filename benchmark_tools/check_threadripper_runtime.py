"""Compose pinned runtime-tree and native Python lookup checks, without admission."""

from pathlib import Path
import subprocess
import time

from benchmark_tools.measure_native_scaling_run import check_manifests
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.inspect_native_python_lookup import compare_lookup
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.probe_dgx_step_separation import save


class RuntimeChecker:
    def __init__(self, receipt_path, receipt_sha256, output, *, tree_checker=check_manifests,
                 runner=subprocess.run):
        self.receipt_path = Path(receipt_path)
        self.receipt_sha256 = receipt_sha256
        self.output = Path(output)
        if not self.output.is_absolute() or self.output.resolve() != self.output:
            raise ValueError("Require direct absolute lookup evidence directory")
        self.output.mkdir(parents=True, exist_ok=False)
        self.tree_checker, self.runner = tree_checker, runner
        self.count = 0

    def __call__(self, specifications):
        self.count += 1
        receipt = read_frozen(self.receipt_path, self.receipt_sha256)
        if receipt["status"] != "native_lookup_repeated_identity_match":
            raise ValueError("Require established native lookup baseline")
        binding_row, baseline_row = receipt["binding"], receipt["baseline"]
        binding = read_frozen(Path(binding_row["path"]), binding_row["sha256"])
        baseline = read_frozen(Path(baseline_row["path"]), baseline_row["sha256"])
        expected = [(str(p), sha) for p, sha in binding["runtime_specs"]]
        if [(str(p), sha) for p, sha in specifications] != expected:
            raise ValueError("Runtime specifications differ from lookup binding")
        source = Path(__file__).with_name("inspect_native_python_lookup.py").resolve()
        if record(source) != receipt["source"]:
            raise ValueError("Lookup inspector source differs from validated version")
        runtime = self.tree_checker(specifications)
        directory = self.output / f"check_{self.count:02d}"
        command = [baseline["tool_entrypoints"]["orthohmm_python"]["absolute_path"], "-B", str(source),
            "--baseline", baseline_row["path"], "--baseline-sha256", baseline_row["sha256"],
            "--binding", binding_row["path"], "--binding-sha256", binding_row["sha256"],
            "--output", str(directory)]
        save(self.output / f"command_{self.count:02d}.json", dict(command=command))
        started = time.monotonic()
        with (self.output / f"process_{self.count:02d}.log").open("x") as log:
            process = self.runner(command, stdout=log, stderr=subprocess.STDOUT, timeout=240)
        save(self.output / f"process_{self.count:02d}.json", dict(exit_code=process.returncode,
                                                                 wall_s=time.monotonic()-started))
        if process.returncode:
            raise RuntimeError("Native lookup inspector failed")
        comparisons = {}
        for name in ("orthohmm", "orthofinder"):
            expected_row = receipt["interpreters"][name]["reports"][-1]
            prior = read_frozen(Path(expected_row["path"]), expected_row["sha256"])
            observed_row = record(directory / (name + ".json"))
            observed = read_frozen(Path(observed_row["path"]), observed_row["sha256"])
            comparisons[name] = dict(compare_lookup(prior, observed), report=observed_row)
        check(receipt["source"])
        read_frozen(self.receipt_path, self.receipt_sha256)
        read_frozen(Path(binding_row["path"]), binding_row["sha256"])
        read_frozen(Path(baseline_row["path"]), baseline_row["sha256"])
        result = dict(status="runtime_and_lookup_checked", runtime=runtime, lookup=comparisons,
            scientific_execution_authorized=False,
            limitations=["Caller must bind run/baseline/collector and enforce resource and quiet-host eligibility.",
                         "Identity checks are outside native timing; neither continuous import enforcement nor admission."])
        save(self.output / f"checked_{self.count:02d}.json", result)
        return result
