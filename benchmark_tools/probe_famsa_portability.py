"""Small alignment portability probe; not benchmark accuracy or timing evidence."""

import argparse
import hashlib
import io
import json
from pathlib import Path
import platform
import subprocess
import sys

import Bio
from Bio import SeqIO


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fixtures():
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    # Hash-derived toy residues avoid random-generator version differences.
    base = "".join(alphabet[b % 20] for i in range(8)
                   for b in hashlib.sha256(f"famsa-probe-{i}".encode()).digest())
    changed = "".join(alphabet[(alphabet.index(c) + 7) % 20] if i % 5 == 0 else c
                      for i, c in enumerate(base))
    return {
        "duplicates": {"a": base, "b": base, "c": changed, "d": changed},
        "indels": {"a": base, "b": base[:80] + base[93:],
                   "c": base[:110] + "ACDEFGHIK" + base[110:],
                   "d": changed[20:], "e": base[:-27]},
        "ambiguous": {"a": base, "b": "X" + base[1:],
                       "c": changed[:90] + "XXX" + changed[93:], "d": base[17:]},
    }


def validate_alignment(text, expected):
    records = list(SeqIO.parse(io.StringIO(text), "fasta"))
    ids = [record.id for record in records]
    if len(ids) != len(set(ids)) or set(ids) != set(expected):
        raise ValueError("Alignment identifiers changed")
    aligned = {record.id: str(record.seq) for record in records}
    if len({len(s) for s in aligned.values()}) != 1:
        raise ValueError("Unequal alignment lengths")
    if any(s.replace("-", "") != expected[name] for name, s in aligned.items()):
        raise ValueError("Input residues changed")
    return aligned


def probe(binary, output):
    binary, output = Path(binary).resolve(), Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    before = sha(binary)
    version = subprocess.run([str(binary), "-h"], capture_output=True, text=True, timeout=30)
    rows = []
    for name, sequences in fixtures().items():
        source = output / f"{name}.fa"
        source.write_text("".join(f">{key}\n{seq}\n" for key, seq in sequences.items()))
        for threads in (1, 4):
            for repeat in (1, 2):
                target = output / f"{name}.t{threads}.r{repeat}.aln"
                argv = [str(binary), "-t", str(threads), str(source), str(target)]
                result = subprocess.run(argv, capture_output=True, text=True, timeout=120)
                row = dict(fixture=name, threads=threads, repeat=repeat, argv=argv,
                           input_sha256=sha(source), exit_code=result.returncode,
                           stdout=result.stdout, stderr=result.stderr, valid=False)
                try:
                    if result.returncode:
                        raise ValueError("FAMSA exited unsuccessfully")
                    row["alignment"] = validate_alignment(target.read_text(), sequences)
                    row["output_sha256"] = sha(target)
                    row["valid"] = True
                except (ValueError, OSError) as exc:
                    row["error"] = str(exc)
                rows.append(row)
    report = dict(machine=platform.machine(), python=sys.version, biopython=Bio.__version__,
                  script_sha256=sha(__file__), binary=str(binary), binary_sha256=before,
                  binary_unchanged=sha(binary) == before,
                  version=dict(exit_code=version.returncode, stdout=version.stdout, stderr=version.stderr),
                  fixtures=fixtures(), rows=rows,
                  limitation="Toy fixtures only; no general equivalence, accuracy or timing claim.")
    (output / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    report = probe(args.binary, args.output)
    sys.exit(0 if report["binary_unchanged"] and all(r["valid"] for r in report["rows"]) else 1)
