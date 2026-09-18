"""Stage verified BPO/index/GG copies in the exact native OrthoMCL mode-4 layout."""

import json
from pathlib import Path
import re
import shutil

from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.prepare_orthomcl_bpo_checkpoint import native_step
from benchmark_tools.prepare_qfo_corrected_orthomcl_native import TOOL
from benchmark_tools.run_qfo_corrected_blast import environment
from benchmark_tools.orthomcl_matrix_to_pairwise import load_species

NAMES = {"bpo": "all.bpo", "offsets": "all_bpo.idx", "ranges": "all_bpo.se", "species": "all.gg"}


def stage(records, directory, expected_proteins, expected_species, expected_indexes):
    if set(records) != set(NAMES):
        raise ValueError("Require exactly BPO, offsets, query ranges and species records")
    if (any(type(v) is not int or v < 1 for v in (expected_proteins, expected_species))
            or not directory.is_absolute() or not re.fullmatch(r"[A-Za-z0-9_./-]+", str(directory))):
        raise ValueError("Require positive scope and absolute shell-safe staging path")
    if directory.exists() or directory.is_symlink():
        raise FileExistsError(directory)
    if len({item["path"] for item in records.values()}) != 4:
        raise ValueError("Require distinct native input sources")
    for item in records.values():
        if type(item["bytes"]) is not int or item["bytes"] <= 0:
            raise ValueError("Empty or invalid native source record")
        check(item)
    directory.mkdir(parents=True, exist_ok=False)
    result = {"status": "staging", "sources": records, "staged": {}, "accuracy_admitted": False,
              "publication_ready": False}
    try:
        for key, name in NAMES.items():
            destination = directory / name
            partial = directory / (name + ".partial")
            shutil.copyfile(records[key]["path"], partial)
            copied = record(partial)
            if any(copied[field] != records[key][field] for field in ("bytes", "sha256")):
                raise ValueError("Copied native input differs: " + key)
            partial.rename(destination)
            result["staged"][key] = record(destination)
        owners = load_species(directory / "all.gg")
        if len(owners) != expected_proteins or len(set(owners.values())) != expected_species:
            raise ValueError("Staged species mapping differs from admitted scope")
        validation = directory / "validation"
        validation.mkdir()
        helpers = Path(__file__).resolve().parent
        argv = [str(TOOL / "venv_orthomcl/bin/perl"), str(helpers / "run_orthomcl_perl_script.pl"),
                str(helpers / "validate_orthomcl_bpo_indexes.pl"), str(directory / "all.bpo"),
                str(directory / "all_bpo.idx"), str(directory / "all_bpo.se")]
        result["index_validation_command"] = argv
        native_step(argv, validation, environment(), "indexes")
        indexes = json.loads((validation / "indexes.stdout").read_text())
        if indexes != expected_indexes:
            raise ValueError("Staged index validation differs from admission")
        for item in [*records.values(), *result["staged"].values()]:
            check(item)
        result.update(status="native_inputs_staged_and_indexes_verified", index_validation=indexes,
                      input_proteins=len(owners), species=len(set(owners.values())))
        result["limitations"] = ["Requires independently admitted input records and runtime checks from the caller.",
                                  "Copies are separate files, not shared writable links; no cache resume is permitted.",
                                  "Native execution and final-group admission remain separate requirements."]
    except BaseException as error:
        result.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        with (directory / "staging.json").open("x") as stream:
            json.dump(result, stream, indent=2, sort_keys=True)
            stream.write("\n")
    return result
