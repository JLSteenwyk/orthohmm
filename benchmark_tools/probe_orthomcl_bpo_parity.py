"""Compare the streaming BPO converter with native OrthoMCL 1.4 on small fixtures."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.convert_orthomcl_blast import convert_blast
from benchmark_tools.prepare_qfo_corrected_orthomcl import SOFTWARE
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_qfo_corrected_blast import environment

TOOL = SOFTWARE / "ORTHOMCLV1.4"


def expected_indexes(path):
    offsets, ranges, closed = [0], {}, set()
    current, start, count = None, None, 0
    with path.open("rb") as stream:
        for count, line in enumerate(stream, 1):
            fields = line.rstrip(b"\r\n").split(b";")
            if len(fields) != 8 or fields[0] != str(count).encode():
                raise ValueError("Malformed or nonsequential BPO record")
            query = fields[1].decode("ascii")
            if not query:
                raise ValueError("Empty BPO query")
            if query != current:
                if query in closed:
                    raise ValueError("Noncontiguous BPO query")
                if current is not None:
                    ranges[current] = f"{start};{count - 1}"
                closed.add(query)
                current, start = query, count
            offsets.append(stream.tell())
    if current is None:
        raise ValueError("Empty BPO")
    ranges[current] = f"{start};{count}"
    return offsets, ranges


def compare(native_bpo, converted_bpo, native):
    offsets, ranges = expected_indexes(native_bpo)
    if native["offsets"] != offsets or native["query_ranges"] != ranges:
        raise ValueError("Native BPO index differs from independent byte/range calculation")
    return {"bpo_byte_identical": native_bpo.read_bytes() == converted_bpo.read_bytes(),
            "native_pair_records": len(offsets) - 1, "native_queries": len(ranges),
            "offset_entries_including_eof": len(offsets), "indexes_verified": True}


def fixture(output):
    fasta, blast = output / "all.fa", output / "all.blast"
    fasta.write_text("".join(f">{key}\n" + "A" * length + "\n"
                             for key, length in (("A", 120), ("B", 100), ("C", 80), ("D", 30))))
    blast.write_text(
        "A\tA\t100.00\t120\t0\t0\t1\t120\t1\t120\t0.0\t200\n"
        "A\tB\t80.00\t10\t2\t0\t1\t10\t20\t29\te-20\t100\n"
        "A\tB\t25.00\t8\t4\t1\t11\t18\t3\t8\t2e-10\t50\n"
        "A\tC\t80.00\t10\t1\t1\t1\t9\t1\t10\t1e-5\t30\n"
        "A\tD\t100.00\t10\t0\t0\t1\t10\t1\t10\t2e-5\t20\n"
        "B\tA\t66.67\t3\t1\t0\t1\t3\t1\t3\t1e-6\t20\n"
        "B\tB\t100.00\t100\t0\t0\t1\t100\t1\t100\t0.0\t200\n"
        "C\tC\t100.00\t80\t0\t0\t1\t80\t1\t80\te-200\t500\n")
    return fasta, blast


def run(output, guarded=False):
    if output.exists():
        raise FileExistsError(output)
    perl = TOOL / "venv_orthomcl/bin/perl"
    driver = Path(__file__).with_name("probe_orthomcl_native_bpo.pl")
    index_checker = Path(__file__).with_name("validate_orthomcl_bpo_indexes.pl")
    checked = [record(p) for p in (perl, TOOL / "orthomcl_module.pm", driver, Path(__file__),
                                   index_checker, Path(__file__).with_name("convert_orthomcl_blast.py"))]
    wrapper = Path(__file__).with_name("run_orthomcl_perl_script.pl")
    if guarded:
        checked.append(record(wrapper))
    output.mkdir(parents=True, exist_ok=False)
    fasta, blast = fixture(output)
    checked.extend([record(fasta), record(blast)])
    argv = [str(perl), "-I" + str(TOOL), *([str(wrapper)] if guarded else []),
            str(driver), str(fasta), str(blast), str(output / "native")]
    report = {"status": "running", "command": argv, "cwd": str(output), "environment": environment(),
              "checked_records": checked, "accuracy_admitted": False, "publication_ready": False,
              "launch_policy": "absolute_module_paths_only" if guarded else "legacy_default"}
    try:
        with (output / "process.log").open("xb") as log:
            done = subprocess.run(argv, cwd=output, env=report["environment"], stdout=log,
                                  stderr=subprocess.STDOUT, timeout=60)
        report["exit_code"] = done.returncode
        if done.returncode:
            raise ValueError("Native BPO fixture failed")
        native = json.loads((output / "native/native.json").read_text())
        if native["orthomcl_version"] != "1.4":
            raise ValueError("Wrong native OrthoMCL version")
        report["loaded_modules"] = [record(p) for p in sorted(set(native["loaded_modules"].values()))]
        count = convert_blast(blast, fasta, output / "streaming.bpo", progress_every=0)
        content = compare(output / "native/native.bpo", output / "streaming.bpo", native)
        content["streaming_pair_records"] = count
        report["index_checker_command"] = [str(perl), *([str(wrapper)] if guarded else []), str(index_checker),
            *[str(output / "native" / name) for name in ("native.bpo", "native.idx", "native.se")]]
        index_result = subprocess.run(report["index_checker_command"], cwd=output,
            env=report["environment"], capture_output=True, timeout=60)
        (output / "index_validation.json").write_bytes(index_result.stdout)
        (output / "index_validation.log").write_bytes(index_result.stderr)
        if index_result.returncode or index_result.stderr:
            raise ValueError("Streaming native index checker failed")
        index_content = json.loads(index_result.stdout)
        if (type(index_content["records"]) is not int or index_content["records"] != content["native_pair_records"]
                or index_content["queries"] != content["native_queries"]
                or index_content["offset_entries_including_eof"] != content["offset_entries_including_eof"]):
            raise ValueError("Independent index checker counts differ")
        report["streaming_index_validation"] = index_content
        report.update(content=content, native_versions={k: v for k, v in native.items() if k.endswith("version")})
        for item in [*checked, *report["loaded_modules"]]:
            check(item)
        report["status"] = "native_bpo_fixture_parity_verified" if content["bpo_byte_identical"] else "native_bpo_fixture_difference_requires_review"
        report["limitations"] = [
            "Small synthetic conversion fixtures, not complete production BPO or biological validation.",
            "Loaded Perl modules are recorded after this probe, not a hermetic production runtime freeze.",
            "Byte parity on these HSP/gap/cutoff cases does not prove parity for all possible inputs.",
            "Independent offsets include the native EOF sentinel; query ranges use inclusive one-based similarity IDs."]
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        report["outputs"] = [record(p) for p in sorted(output.rglob("*")) if p.is_file()]
        with (output / "report.json").open("x") as stream:
            json.dump(report, stream, indent=2, sort_keys=True)
            stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--guarded", action="store_true")
    args = parser.parse_args()
    report = run(args.output.resolve(), args.guarded)
    print(json.dumps({"status": report["status"], "content": report.get("content")}))
