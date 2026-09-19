"""Same-host isolated CPU wheel smoke test; not benchmark equivalence."""

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import SeqIO
from benchmark_tools.audit_wheel_contents import audit
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check

PROBE = """
import ctypes, importlib.metadata, json, sys
from pathlib import Path
import orthohmm
prefix = Path(sys.prefix).resolve()
module = Path(orthohmm.__file__).resolve()
assert module.is_relative_to(prefix)
libraries = sorted(module.parent.joinpath('search/csrc').glob('*.so'))
assert [p.name for p in libraries] == ['hmm_viterbi.so', 'kmer_prefilter.so', 'pair_align.so']
for path in libraries:
    ctypes.CDLL(str(path))
dependencies = []
for dist in importlib.metadata.distributions():
    assert Path(dist.locate_file('')).resolve().is_relative_to(prefix)
    dependencies.append({'name':dist.metadata['Name'], 'version':dist.version})
print(json.dumps({'module':str(module), 'prefix':str(prefix), 'python':sys.version,
                  'dependencies':sorted(dependencies,key=lambda d:d['name'].lower()),
                  'native_libraries_loaded':[str(p) for p in libraries]}))
"""


def partition(path, genes):
    groups = []
    for line in path.read_text().splitlines():
        group, members = line.split(":", 1)
        if not group or not members.split():
            raise ValueError("Empty native group")
        groups.append(members.split())
    assigned = [g for members in groups for g in members]
    if len(genes) != len(set(genes)) or sorted(assigned) != sorted(genes):
        raise ValueError("Native partition does not cover each input exactly once")
    return dict(groups=len(groups), genes=len(assigned), partition=record(path))


def verify(root, python, wheel, output):
    output.mkdir(parents=True, exist_ok=False)
    wheel_record = record(wheel)
    inventory = audit(wheel, root / "LICENSE.md")
    native = [r for r in inventory["entries"] if r["path"].endswith(".so")]
    if (not inventory["embedded_project_license_matches"] or len(native) != 3
            or any("cuda" in r["path"].lower() or r["cuda_runtime_symbol_witnesses"] for r in native)):
        raise ValueError("Require CPU-only project wheel")
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", PYTHONHASHSEED="0")
    probe = json.loads(subprocess.check_output([str(python), "-I", "-c", PROBE], cwd=output, env=env, text=True))
    site = Path(probe["module"]).parent.parent
    installed = []
    with zipfile.ZipFile(wheel) as archive:
        for name in archive.namelist():
            if name.startswith("orthohmm/") and not name.endswith("/"):
                path = site / name
                if path.read_bytes() != archive.read(name):
                    raise ValueError("Installed package differs from wheel: " + name)
                installed.append(record(path))
    inputs = output / "inputs"
    inputs.mkdir()
    genes, originals = [], []
    for path in sorted((root / "tests/samples").glob("*.fa")):
        originals.append(record(path))
        shutil.copyfile(path, inputs / path.name)
        genes.extend(r.id for r in SeqIO.parse(path, "fasta"))
    if not genes:
        raise ValueError("Missing smoke inputs")
    copies = [record(p) for p in sorted(inputs.iterdir())]
    runs = []
    for profile in ("standard", "high_sensitivity"):
        destination = output / profile
        destination.mkdir()
        argv = [str(python), "-I", "-m", "orthohmm", str(inputs), "-o", str(destination),
                "-c", "1", "--search_mode", "builtin", "--clustering", "leiden", "--accuracy_profile", profile]
        log = output / (profile + ".log")
        with log.open("x") as stream:
            subprocess.run(argv, cwd=output, env=env, stdout=stream, stderr=subprocess.STDOUT, timeout=180, check=True)
        runs.append(dict(profile=profile, argv=argv, log=record(log), **partition(destination / "orthohmm_orthogroups.txt", genes)))
    for item in [wheel_record, *installed, *originals, *copies]:
        check(item)
    result = dict(status="same_host_cpu_wheel_smoke_verified", wheel=wheel_record,
                  inventory=inventory, runtime=probe, installed=installed, original_inputs=originals,
                  inputs=copies, runs=runs, source=record(__file__), publication_ready=False,
                  limitations=["Current development wheel, not the frozen scientific executor.",
                               "Same-host fixture only, not portability or biological accuracy.",
                               "Does not execute the full phylogenetic pipeline or clear redistribution rights."])
    with (output / "verification.json").open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "python", "wheel", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    verify(args.root.resolve(), args.python.absolute(), args.wheel.resolve(), args.output.absolute())
