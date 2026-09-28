"""Audit all retained FAS annotation graphs without running pairwise scoring."""

import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path

from count_native_fas_paths import run
from probe_native_fas_omission import fingerprint


def one_file(task):
    source, directory = task
    result = run(Path(source["path"]))
    if result["annotation"] != source:
        raise ValueError("Annotation differs from frozen environment")
    output = directory / (Path(source["path"]).stem + ".json")
    with output.open("x") as handle:
        json.dump(result, handle, sort_keys=True)
        handle.write("\n")
    return dict(annotation=source, result=fingerprint(output), summary=result["summary"],
                options=result["options"], greedyfas_version=result["greedyfas_version"])


def panel(environment, directory):
    env_ref = fingerprint(environment)
    env = json.loads(environment.read_text())
    inputs = [r for r in env["reference_files"] if "/fas_annotations/" in r["path"]
              and r["path"].endswith(".json")]
    if not inputs or len({r["path"] for r in inputs}) != len(inputs):
        raise ValueError("Missing or duplicate frozen annotation inventory")
    refs = [env_ref, fingerprint(__file__), fingerprint(Path(__file__).with_name("count_native_fas_paths.py")),
            fingerprint(Path(__file__).with_name("probe_native_fas_omission.py")), *inputs]
    if refs != [fingerprint(r["path"]) for r in refs]:
        raise ValueError("Frozen annotation inputs changed")
    directory.mkdir(parents=True, exist_ok=False)
    rows = []
    with ProcessPoolExecutor(max_workers=2) as pool:
        for row in pool.map(one_file, [(r, directory) for r in sorted(inputs, key=lambda r: r["path"])]):
            rows.append(row)
            print(json.dumps(dict(annotation=Path(row["annotation"]["path"]).name, **row["summary"])), flush=True)
    if refs != [fingerprint(r["path"]) for r in refs]:
        raise ValueError("Inputs changed during panel")
    if any(r["options"] != rows[0]["options"] or r["greedyfas_version"] != rows[0]["greedyfas_version"] for r in rows):
        raise ValueError("Inconsistent effective options or version")
    return dict(status="retained_annotation_path_panel_complete", checked_inputs=refs, workers=2,
        files=rows, summary=dict(annotation_files=len(rows),
            protein_records=sum(r["summary"]["proteins"] for r in rows),
            rejected_records=sum(r["summary"]["rejected"] for r in rows)),
        historical_omissions_attributed=False, benchmark_scores_changed=False,
        limitations=["Protein records may recur across annotation files; totals are not asserted unique proteins.",
            "Counts do not reconstruct historical sampled pairs or attribute their omissions.",
            "This is a shared-host graph diagnostic, not controlled runtime evidence."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--environment", required=True, type=Path)
    parser.add_argument("--directory", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    result = panel(args.environment.resolve(), args.directory.resolve())
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
