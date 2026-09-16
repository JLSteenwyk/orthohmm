"""Exercise profile construction in the exact inference checkout and interpreter."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys


PROBE = r'''
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

root = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(root))
result = {"status": "failed", "root": str(root), "python": sys.executable}
try:
    from orthohmm.search import msa_profile, msa_center_star
    from orthohmm.search.matrices import get_matrix, get_background_freqs
    actual = Path(msa_profile.__file__).resolve()
    if actual != root / "orthohmm/search/msa_profile.py":
        raise ValueError("Profile implementation imported from wrong checkout")
    result["profile_source"] = str(actual)
    result["profile_source_sha256"] = hashlib.sha256(actual.read_bytes()).hexdigest()
    library = root / "orthohmm/search/csrc/pair_align.so"
    result["pair_align"] = {"path": str(library), "exists": library.is_file()}
    if library.is_file():
        result["pair_align"]["sha256"] = hashlib.sha256(library.read_bytes()).hexdigest()
    msa_center_star._load_pair_align()
    sequences = ["ACDEFGHIKLMNPQRSTVWY", "ACDEFGHIKLMNPQRSTVWY", "ACDEFGHIKLMNPQKSTVWY"]
    profile = msa_profile.build_msa_profile(
        sequences, ["probe_a", "probe_b", "probe_c"],
        get_matrix("BLOSUM62"), get_background_freqs("BLOSUM62"))
    if profile is None or profile.length < 1:
        raise ValueError("Nonempty synthetic cluster produced no profile")
    if profile.match_emissions.shape != (profile.length, 20):
        raise ValueError("Invalid profile emission dimensions")
    if not np.isfinite(profile.match_emissions).all():
        raise ValueError("Nonfinite profile emissions")
    result.update(status="passed", profile_length=profile.length)
except Exception as error:
    result.update(error_type=type(error).__name__, error=str(error))
print(json.dumps(result, sort_keys=True))
raise SystemExit(0 if result["status"] == "passed" else 1)
'''


def probe_profile_runtime(root, python=sys.executable):
    root = Path(root).resolve()
    env = os.environ.copy()
    env.update(PYTHONPATH=str(root), OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1")
    completed = subprocess.run([str(python), "-c", PROBE, str(root)], cwd=root,
                               env=env, capture_output=True, text=True, timeout=60)
    try:
        report = json.loads(completed.stdout)
    except (ValueError, TypeError) as error:
        raise ValueError("Profile runtime probe did not return valid JSON") from error
    if not isinstance(report, dict) or report.get("status") not in {"passed", "failed"}:
        raise ValueError("Invalid profile runtime probe status")
    report.update(exit_code=completed.returncode, stderr=completed.stderr)
    if completed.returncode != 0 and report["status"] == "passed":
        raise ValueError("Profile runtime probe reported success with nonzero exit")
    return report


def require_profile_runtime(root, python=sys.executable):
    report = probe_profile_runtime(root, python)
    if report["status"] != "passed" or report["exit_code"] != 0:
        raise ValueError("Profile runtime unavailable: " + json.dumps(report, sort_keys=True))
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError("Refusing to overwrite runtime evidence")
    report = probe_profile_runtime(args.root, args.python)
    with args.output.open("x") as handle:
        handle.write(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return 0 if report["status"] == "passed" and report["exit_code"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
