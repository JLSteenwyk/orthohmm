"""One no-site-import interpreter control for retained Memcheck startup reports."""

import argparse
import json
from pathlib import Path
import subprocess
import time

from benchmark_tools.diagnose_cpm_refinement_allocator import environment, memcheck_command
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.summarize_cpm_memcheck import digest, summarize

STATUS_SHA = "33d543f8d47b1df8175030aa209a1ea19406545f37401ba1b7c3f467ebea3d3a"


def run(status_path, output):
    if digest(status_path) != STATUS_SHA:
        raise ValueError("Changed retained diagnostic status")
    original = json.loads(status_path.read_text())
    child = original["scientific_child_command"]
    launcher = Path(child[2]).parents[1]
    python = Path(child[0])
    wanted = {str(python.resolve()), "/usr/lib/x86_64-linux-gnu/libc.so.6"}
    selected = {r["path"]: r for r in original["checked_records"] if r["path"] in wanted}
    if set(selected) != wanted:
        raise ValueError("Missing retained interpreter/libc identities")
    checked = [record(status_path), record(__file__), *selected.values(),
               original["memchecker"]["binary"], *original["memchecker"]["components"]]
    for item in checked:
        if record(item["path"]) != item:
            raise ValueError("Changed control executable/input")
    output.mkdir(parents=True, exist_ok=False)
    env = environment(launcher, native_memcheck=True)
    command = memcheck_command(original["memchecker"]["binary"]["path"],
                               [str(python), "-B", "-S", "-c", "pass"], output / "memcheck.xml")
    report = dict(status="startup_control_running", command=command, cwd=str(launcher),
                  attempts=1, timeout_seconds=120, checked_records=checked,
                  environment_overrides={k: env[k] for k in (
                      "PYTHONPATH", "PYTHONHASHSEED", "PYTHONNOUSERSITE", "OMP_NUM_THREADS",
                      "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "PYTHONMALLOC", "PYTHONFAULTHANDLER")},
                  removed_environment_keys=["PYTHONHOME", "LD_PRELOAD", "LD_LIBRARY_PATH",
                                            "VALGRIND_LIB", "VALGRIND_OPTS"],
                  accuracy_admitted=False, publication_ready=False)
    start = time.monotonic()
    try:
        with (output / "child.log").open("xb") as log:
            result = subprocess.run(command, cwd=launcher, env=env, timeout=120,
                                    stdout=log, stderr=subprocess.STDOUT)
        report["returncode"] = result.returncode
        xml = output / "memcheck.xml"
        report["summary"] = summarize(xml, digest(xml))
        report["xml"] = record(xml)
        for item in checked:
            if record(item["path"]) != item:
                raise ValueError("Control input changed during execution")
        report["status"] = "startup_control_observed"
    except Exception as error:
        report.update(status="startup_control_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        report["wall_seconds_descriptive_only"] = time.monotonic() - start
        report["limitations"] = [
            "No-site-import startup is not a refinement or scientific-runtime validation.",
            "Reproduced reports do not prove false positives or explain the earlier SIGSEGV.",
            "One local diagnostic, not comparative timing or repair of failed admission."]
        with (output / "report.json").open("x") as stream:
            json.dump(report, stream, indent=2, sort_keys=True)
            stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--status", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    run(args.status.resolve(), args.output.absolute())
