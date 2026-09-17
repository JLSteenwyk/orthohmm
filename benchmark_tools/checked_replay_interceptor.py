"""Preserve and check each frozen isolated-clustering payload without changing graph generation."""

import json
from pathlib import Path
import shutil
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent))
from prepare_ob_candidate_neighborhood import check, record
from checked_replay_payload_worker import FILES, STAGES


class CheckedReplaySubprocess:
    def __init__(self, original, root, output, worker_path, validate_result):
        self.original, self.root, self.output = original, root, output
        self.worker_path, self.validate_result = worker_path, validate_result
        self.calls = []

    def __getattr__(self, name):
        return getattr(self.original, name)

    def run(self, command, *args, **kwargs):
        prefix = [sys.executable, "-m", "orthohmm.leiden_worker"]
        if not isinstance(command, list) or command[:3] != prefix:
            return self.original.run(command, *args, **kwargs)
        if (len(command) != 4 or args or kwargs != {"check": True} or len(self.calls) >= len(STAGES)
                or any(row["status"] != "checked" for row in self.calls)):
            raise ValueError("Unexpected frozen native worker invocation")
        index = len(self.calls)
        directory = self.output / f"cluster_{index}_{STAGES[index]}"
        directory.mkdir(parents=True, exist_ok=False)
        original_payload = Path(command[3]).resolve()
        if {p.name for p in original_payload.iterdir()} != set(FILES):
            raise ValueError("Unexpected original payload file inventory")
        payload = directory / "payload"
        shutil.copytree(original_payload, payload)
        original_records = [record(original_payload / name) for name in FILES]
        inputs = [record(payload / name) for name in FILES]
        if any(any(a[key] != b[key] for key in ("bytes", "sha256")) for a, b in zip(original_records, inputs)):
            raise ValueError("Payload copy changed content")
        metadata = json.loads((payload / "metadata.json").read_text())
        manifest = {"stage": STAGES[index], "index": index, "accuracy_evaluated": False,
                    "inputs": inputs, "output_directory": metadata["output_directory"],
                    "original_command": command, "original_payload": original_records}
        manifest_path = directory / "payload_manifest.json"
        manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        adapted = [sys.executable, str(self.worker_path), "--root", str(self.root), "--payload", str(payload),
                   "--manifest", str(manifest_path), "--manifest-sha256", record(manifest_path)["sha256"]]
        row = {"index": index, "stage": STAGES[index], "status": "running", "command": adapted,
               "manifest": record(manifest_path), "accuracy_evaluated": False}
        self.calls.append(row)
        started = time.monotonic()
        try:
            with (directory / "worker.log").open("x") as log:
                result = self.original.run(adapted, stdout=log, stderr=self.original.STDOUT, check=True)
            for item in [*inputs, *original_records]:
                check(item)
            # Validate the worker before allowing the frozen replay to use its output.
            row["validation"] = self.validate_result(payload, manifest)
            partition = Path(metadata["output_directory"]) / "orthohmm_working_res/orthohmm_edges_clustered.txt"
            retained = directory / "partition.txt"
            shutil.copyfile(partition, retained)
            if record(partition)["sha256"] != record(retained)["sha256"]:
                raise ValueError("Partition copy changed")
            row.update(status="checked", partition=record(retained), exit_code=result.returncode)
            return result
        except BaseException as error:
            row.update(status="failed", error_type=type(error).__name__, error=str(error))
            raise
        finally:
            row["wall_s"] = time.monotonic() - started
            (directory / "execution.json").write_text(json.dumps(row, indent=2, sort_keys=True) + "\n")
