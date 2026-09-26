"""One-attempt orchestration of a frozen matched-search graph cell."""

import argparse
import json
import os
from pathlib import Path
import shutil
import signal
import time

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_search_sensitivity_cell import execute, interrupted, INSTALL_SHA


MANIFEST_SHA = "d0012cc2ad65dbd1b5272c340d76edfd6c459c2c98213955495c063acd7ecafb"


def run(manifest_path, install_path, index, output):
    manifest_record, install_record = record(manifest_path), record(install_path)
    if manifest_record["sha256"] != MANIFEST_SHA or install_record["sha256"] != INSTALL_SHA:
        raise ValueError("Unrecognized frozen manifest/runtime")
    manifest, installed = json.loads(manifest_path.read_text()), json.loads(install_path.read_text())
    if len(manifest["cells"]) != 70 or not 0 <= index < 70:
        raise ValueError("Invalid cell index/panel size")
    cell = manifest["cells"][index]
    output.mkdir(parents=True, exist_ok=False)
    result = dict(status="failed", index=index, cell=cell, attempt=1,
                  manifest=manifest_record, installation=install_record, stages=[],
                  slurm_job_id=os.environ.get("SLURM_JOB_ID"),
                  slurm_array_job_id=os.environ.get("SLURM_ARRAY_JOB_ID"),
                  slurm_array_task_id=os.environ.get("SLURM_ARRAY_TASK_ID"),
                  timing_comparability="incremental_shared_host_descriptive_only")
    try:
        staging = next(Path(r["path"]) for r in installed["checked_records"] if Path(r["path"]).name == "staging.json")
        python = staging.parent / "venv_clean/bin/python"
        worker = Path(__file__).with_name("matched_graph_worker.py").resolve()
        checked = [manifest_record, install_record, manifest["protocol"], manifest["search_result"],
                   *installed["checked_records"], cell["numeric"], *cell["sources"],
                   record(__file__), record(worker), record(python), record("/usr/bin/time"),
                   record(Path(__file__).with_name("run_search_sensitivity_cell.py")),
                   record(Path(__file__).with_name("prepare_ob_candidate_neighborhood.py"))]
        result["checked_records"] = checked
        for item in checked:
            check(item)
        private = output / "numeric.json"
        shutil.copyfile(cell["numeric"]["path"], private)
        private_record = record(private)
        if (private_record["bytes"], private_record["sha256"]) != (cell["numeric"]["bytes"], cell["numeric"]["sha256"]):
            raise ValueError("Private numeric copy differs")
        result["private_numeric"] = private_record
        env = dict(PATH="/usr/bin:/bin", HOME=str(output), LC_ALL="C", OMP_NUM_THREADS="1",
                   OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", PYTHONHASHSEED="0")
        result["environment"] = env
        command = [str(python), "-I", str(worker), "--numeric", str(private), "--output", str(output / "graph")]
        execute(command, output, "graph", env, time.monotonic() + 3300, result["stages"])
        native = json.loads((output / "graph/receipt.json").read_text())
        if native["settings"] != manifest["graph_settings"] or native["numeric"] != private_record:
            raise ValueError("Native graph settings/input mismatch")
        for item in [*checked, private_record, *native["outputs"], native["checkpoint_manifest"]]:
            check(item)
        result["native_receipt"] = record(output / "graph/receipt.json")
        result["status"] = "native_completed_pending_independent_readback"
    except BaseException as error:
        result["error"] = f"{type(error).__name__}: {error}"
        raise
    finally:
        result["logs"] = [record(p) for p in sorted(output.iterdir())
                          if p.is_file() and (p.name.endswith(".log") or p.name.endswith(".time.txt"))]
        with (output / "execution.json").open("x") as stream:
            json.dump(result, stream, indent=2, sort_keys=True)
            stream.write("\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--installation", type=Path, required=True)
    parser.add_argument("--index", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    signal.signal(signal.SIGTERM, interrupted)
    signal.signal(signal.SIGINT, interrupted)
    run(args.manifest.resolve(), args.installation.resolve(), args.index, args.output.absolute())
