"""Pin the current presentation inventory; no scientific results are recomputed."""

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "benchmark_tools/results"
HISTORICAL = RESULTS / "publication_main_with_figures_20261002_v3/assembly.json"
HISTORICAL_SHA = "861b869a4e359106843692c50a3e43e7805c177b32cce151fcf9cc2ab271176e"
MAIN = RESULTS / "publication_main_print_20261004/document.pdf"
MAIN_SHA = "32edd18547572fb61f84af7d6627cc2ef96eb0330db8e52f1350bf1f7f56ae98"


def record(path):
    content = path.read_bytes()
    return {"path": str(path.relative_to(ROOT)), "bytes": len(content),
            "sha256": hashlib.sha256(content).hexdigest()}


def accepted(reference):
    relative = Path(reference["path"].split("/benchmark_tools/results/", 1)[1])
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError("Unsafe retained results path")
    observed = record(RESULTS / relative)
    if (observed["bytes"], observed["sha256"]) != (reference["bytes"], reference["sha256"]):
        raise ValueError("Retained PDF checksum mismatch")
    return observed


def build(output):
    history, main = record(HISTORICAL), record(MAIN)
    if history["sha256"] != HISTORICAL_SHA or main["sha256"] != MAIN_SHA:
        raise ValueError("Unexpected historical receipt or current manuscript")
    prior = json.loads(HISTORICAL.read_text())
    figures = [{"number": item["number"], "title": item["title"], "caption": item["caption"],
                "pdf": accepted(item["pdf"]), "link_aliases": [item["pdf"]["path"]]}
               for item in prior["figures"]]
    if [item["number"] for item in figures] != list(range(1, 15)):
        raise ValueError("Unexpected retained figure inventory")
    provenance = [history, record(RESULTS / "publication_main_final_visual_review_20261004.json")]
    additional = [
        ("ob_complete_strata_figure_20260928", "ob_complete_strata.pdf",
         "All-method OrthoBench error strata",
         "Fourteen frozen descriptive strata retain all eight methods, precision/recall/F1 and empty bins. "
         "Development exposure and small bins limit interpretation. Length, copy number and composition "
         "are descriptors, not validated fragments, duplication histories or causal error mechanisms."),
        ("threadripper_shared_resource_figure_20261004_v26", "shared_threadripper_resources.pdf",
         "Completed shared-host resource panel",
         "All 27 attempts are reviewed: 25 measured, 24 eligible; six cells have three eligible repeats. "
         "Three cells lack a third eligible repeat and remain unavailable. Matched CPU/RAM limits do not "
         "remove unknown, potentially tool-dependent contention; no isolated efficiency ranking follows."),
    ]
    for directory, filename, title, caption in additional:
        manifest_path = RESULTS / directory / "manifest.json"
        manifest = json.loads(manifest_path.read_text())
        path = RESULTS / directory / filename
        ref = record(path)
        outputs = manifest["outputs"]
        matching = [item for item in outputs if Path(item["path"]).name == filename]
        if len(matching) != 1 or (ref["bytes"], ref["sha256"]) != (
                matching[0]["bytes"], matching[0]["sha256"]):
            raise ValueError("Figure does not match producer manifest: " + filename)
        figures.append({"number": len(figures) + 1, "title": title, "caption": caption,
                        "pdf": ref, "link_aliases": [str(path)]})
        provenance.append(record(manifest_path))
    selection = {"schema": "publication_figure_selection_v1", "main": main, "main_pages": 12,
                 "figures": figures, "provenance": provenance,
                 "scientific_settings_or_results_changed": False, "publication_ready": False}
    with output.open("x") as stream:
        json.dump(selection, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps(record(output), indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    build(parser.parse_args().output.resolve())
