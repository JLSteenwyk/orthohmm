"""Build and audit a fresh BPO checkpoint; does not admit the source search."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.convert_orthomcl_blast import convert_blast
from benchmark_tools.audit_orthomcl_bpo_content import audit
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.prepare_qfo_corrected_orthomcl_native import TOOL, RUNTIME_SHA, HELPERS_SHA
from benchmark_tools.run_qfo_corrected_blast import environment
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.snapshot_runtime_trees import verify


def native_step(argv, directory, env, name):
    with (directory / (name + ".stdout")).open("xb") as out, (directory / (name + ".stderr")).open("xb") as err:
        done = subprocess.run(argv, cwd=directory, env=env, stdout=out, stderr=err)
    if done.returncode:
        raise ValueError(f"Native {name} failed with exit code {done.returncode}")
    if (directory / (name + ".stderr")).stat().st_size:
        raise ValueError(f"Native {name} diagnostics require review")


def prepare(root, blast, fasta, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    results = root / "benchmark_tools/results"
    paths = [results / "qfo_corrected_orthomcl_perl_runtime_20260918.json",
             results / "qfo_corrected_orthomcl_system_helpers_20260918.json"]
    manifests = [read_frozen(p, sha) for p, sha in zip(paths, (RUNTIME_SHA, HELPERS_SHA))]
    for manifest in manifests:
        verify(manifest)
    helpers = Path(__file__).resolve().parent
    checked = [record(p) for p in [blast, fasta, *paths, Path(__file__),
        *[helpers / name for name in ("build_orthomcl_bpo_indexes.pl", "validate_orthomcl_bpo_indexes.pl",
            "run_orthomcl_perl_script.pl", "convert_orthomcl_blast.py", "audit_orthomcl_bpo_content.py")]]]
    output.mkdir(parents=True, exist_ok=False)
    env = environment()
    report = {"status": "running", "checked_records": checked, "environment": env,
              "commands": [], "search_admitted": False, "accuracy_admitted": False,
              "publication_ready": False}
    try:
        bpo = output / "all.bpo"
        count = convert_blast(blast, fasta, bpo, evalue_cutoff=1e-5)
        content = audit(blast, fasta, bpo, output / "content.json")
        if count != content["content"]["bpo_pair_records"]:
            raise ValueError("Converter and content auditor record counts differ")
        report["content"] = content["content"]
        prefix = [str(TOOL / "venv_orthomcl/bin/perl"), "-I" + str(TOOL),
                  str(helpers / "run_orthomcl_perl_script.pl")]
        indexes = output / "indexes"
        commands = [
            ("build", prefix + [str(helpers / "build_orthomcl_bpo_indexes.pl"), str(bpo), str(indexes)]),
            ("validate", prefix + [str(helpers / "validate_orthomcl_bpo_indexes.pl"), str(bpo),
                                   str(indexes / "all_bpo.idx"), str(indexes / "all_bpo.se")])]
        for name, argv in commands:
            report["commands"].append({"stage": name, "argv": argv, "cwd": str(output)})
            native_step(argv, output, env, name)
        index_result = json.loads((output / "validate.stdout").read_text())
        if (index_result["status"] != "native_bpo_indexes_verified"
                or type(index_result["records"]) is not int or index_result["records"] != count
                or index_result["bpo_bytes"] != bpo.stat().st_size
                or index_result["offset_entries_including_eof"] != count + 1):
            raise ValueError("Native index validation summary differs from audited BPO")
        report["index_validation"] = index_result
        for manifest in manifests:
            verify(manifest)
        for item in checked:
            check(item)
        report["outputs"] = [record(p) for p in sorted(output.rglob("*")) if p.is_file()]
        report["status"] = "bpo_checkpoint_content_and_indexes_verified"
        report["limitations"] = [
            "Source search provenance, terminal execution and query failures must be admitted separately.",
            "Fixed legacy E-value cutoff 1e-5; native index validation is not orthology accuracy validation.",
            "Native indexing loads offsets and query ranges into memory; no constant-memory claim.",
            "Partial files are retained on failure; existing output directories are never resumed or overwritten.",
            "Production scheduler and Python runtime binding remain the responsibility of the execution wrapper."]
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        with (output / "report.json").open("x") as stream:
            json.dump(report, stream, indent=2, sort_keys=True)
            stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "blast", "fasta", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    report = prepare(*(getattr(args, name).resolve() for name in ("root", "blast", "fasta", "output")))
    print(json.dumps({"status": report["status"], "content": report["content"]}))
