"""Isolated indexing-only fixtures for the pinned Proteinortho symbol policy."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


RUNTIME_SHA = "873e089dd4ffaf78540b2d46458f750ff5c5bcc04156e7c6f68c1f473570a628"
IMAGE_SHA = "990be9d066a02302fd49520abfa28d63bfdb5ad0c4792a99ac75fc6152cbe561"
SCRIPT_SHA = "b948dcb6d94059580cd16da0426bc3c1b533ffb88d86560bdfe591ef9abb8ec9"


def sequences():
    clean = "ACDEFGHIKLMNPQRSTVWY" * 4
    return dict(clean=clean, internal_asterisk=clean[:20] + "*" + clean[20:], terminal_asterisk=clean + "*")


def validate(case, returncode, log, unchanged):
    if not unchanged:
        return False
    if case == "clean":
        return returncode == 0 and "Invalid symbol" not in log
    return returncode != 0 and "Invalid symbol" in log and "sanitze with" in log


def run(image, output):
    output.mkdir(parents=True, exist_ok=False)
    runtime = Path("/usr/local/bin/singularity")
    checked = [record(runtime), record(image), record(__file__)]
    if checked[0]["sha256"] != RUNTIME_SHA or checked[1]["sha256"] != IMAGE_SHA:
        raise ValueError("Unrecognized container/runtime")
    env = dict(PATH="/usr/local/bin:/usr/bin:/bin", HOME=str(output), USER="bizon", LOGNAME="bizon", LC_ALL="C", LANG="C")
    source_command = [str(runtime), "exec", "--cleanenv", "--containall", str(image), "cat", "/usr/local/bin/proteinortho"]
    source = subprocess.run(source_command, env=env, cwd=output, capture_output=True, check=True, timeout=60)
    if hashlib.sha256(source.stdout).hexdigest() != SCRIPT_SHA:
        raise ValueError("Unexpected native source")
    (output / "proteinortho.source.pl").write_bytes(source.stdout)
    checked.append(record(output / "proteinortho.source.pl"))
    excerpts = [{"line": i + 1, "text": line} for i, line in enumerate(source.stdout.decode().splitlines())
                if i + 1 == 767 or 5080 <= i + 1 <= 5097]
    rows = []
    for case, sequence in sequences().items():
        directory = output / case
        directory.mkdir()
        (directory / "a.fa").write_text(">a\n" + sequence + "\n")
        (directory / "b.fa").write_text(">b\n" + sequences()["clean"] + "\n")
        inputs = [record(directory / name) for name in ("a.fa", "b.fa")]
        command = [str(runtime), "exec", "--cleanenv", "--containall", "--bind", str(directory) + ":/work",
                   "--pwd", "/work", str(image), "proteinortho", "-step=1", "-cpus=1", "-project=fixture", "/work/a.fa", "/work/b.fa"]
        process = subprocess.Popen(command, env=env, cwd=directory, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, text=True, start_new_session=True)
        try:
            stdout, stderr = process.communicate(timeout=90)
            log = stdout + stderr
            (directory / "native.log").write_text(log)
            unchanged = all(record(r["path"]) == r for r in inputs)
            rows.append(dict(case=case, command=command, inputs=inputs, returncode=process.returncode,
                             log=record(directory / "native.log"), inputs_unchanged=unchanged,
                             expectation_passed=validate(case, process.returncode, log, unchanged)))
        except subprocess.TimeoutExpired as error:
            os.killpg(process.pid, signal.SIGKILL)
            stdout, stderr = process.communicate()
            (directory / "native.log").write_text(stdout + stderr)
            rows.append(dict(case=case, command=command, inputs=inputs, expectation_passed=False,
                             error="Timed out; no automatic retry", timeout_seconds=error.timeout))
            break
    for item in checked:
        check(item)
    result = dict(status="proteinortho_symbol_policy_fixture", source=record(__file__), checked_records=checked,
                  source_command=source_command, source_excerpts=excerpts, environment=env, cells=rows,
                  all_expectations_passed=len(rows) == 3 and all(r["expectation_passed"] for r in rows),
                  historical_cause_proven=False, benchmark_inference_rerun=False, publication_ready=False,
                  limitations=["Indexing-only synthetic fixtures; no pair search, clustering or benchmark accuracy measured.",
                               "Pinned current runtime behavior does not establish historical preprocessing actions or binary identity.",
                               "No input normalization or matched-input promotion applied to retained benchmark results."])
    with (output / "receipt.json").open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--image", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.image.resolve(), args.output.absolute())
    raise SystemExit(0 if result["all_expectations_passed"] else 1)
