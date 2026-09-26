"""Label-free installed-HMM executor; invoke with the pinned Python and -I."""

import argparse
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import sys


SETTINGS = dict(band_width=64, kmer_k=4, evalue_threshold=1e-4,
                threads_per_worker=1, max_candidates_per_query=100)


def fingerprint(path):
    path = Path(path).resolve()
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return dict(path=str(path), bytes=path.stat().st_size, sha256=digest.hexdigest())


def input_inventory(directory):
    paths = sorted(directory.iterdir())
    if not paths or any(p.is_symlink() or not p.is_file() or p.suffix not in {".fasta", ".fa", ".faa"} for p in paths):
        raise ValueError("Expected only regular species FASTA files")
    records, ids, seen = [], {}, set()
    for path in paths:
        records.append(fingerprint(path))
        names, length = [], 0
        with path.open() as stream:
            for line in stream:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if names and not length:
                        raise ValueError("Empty sequence")
                    fields = line[1:].split()
                    if not fields or fields[0] in seen:
                        raise ValueError("Missing or repeated sequence ID")
                    names.append(fields[0])
                    seen.add(fields[0])
                    length = 0
                else:
                    if not names or any(c.isspace() for c in line):
                        raise ValueError("Malformed FASTA sequence")
                    length += len(line)
        if not names or not length:
            raise ValueError("Empty FASTA or sequence")
        ids[path.name] = names
    return records, ids


def write_hits(result, expected_ids, output):
    if result.species_ids != expected_ids:
        raise ValueError("Engine gene IDs differ from FASTA order")
    if set(result.pair_results) != set(itertools.product(expected_ids, repeat=2)):
        raise ValueError("Incomplete species-pair search")
    counts = []
    with output.open("x", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["query_species", "target_species", "query_id", "target_id", "score", "evalue"])
        for key, pair in sorted(result.pair_results.items()):
            if key != (pair.query_species, pair.target_species):
                raise ValueError("Species-pair identity mismatch")
            sizes = {len(pair.query_indices), len(pair.target_indices), len(pair.scores), len(pair.evalues)}
            if len(sizes) != 1:
                raise ValueError("Hit array lengths differ")
            candidates = int(pair.candidate_count)
            if candidates != pair.candidate_count or candidates < len(pair.scores):
                raise ValueError("Invalid candidate count")
            seen = set()
            hits = []
            for qi, ti, score, evalue in zip(pair.query_indices, pair.target_indices, pair.scores, pair.evalues):
                q, t = int(qi), int(ti)
                if (q != qi or t != ti or not 0 <= q < len(expected_ids[key[0]])
                        or not 0 <= t < len(expected_ids[key[1]])):
                    raise ValueError("Invalid local gene index")
                if (q, t) in seen:
                    raise ValueError("Repeated directed hit")
                seen.add((q, t))
                score, evalue = float(score), float(evalue)
                if not math.isfinite(score) or not math.isfinite(evalue) or not 0 <= evalue < SETTINGS["evalue_threshold"]:
                    raise ValueError("Invalid or non-significant native hit")
                hits.append((q, t, score, evalue))
            for q, t, score, evalue in sorted(hits):
                writer.writerow([*key, expected_ids[key[0]][q], expected_ids[key[1]][t], repr(score), repr(evalue)])
            counts.append(dict(query=key[0], target=key[1], candidates=candidates, hits=len(hits)))
    return counts


def run(directory, output):
    if not sys.flags.isolated:
        raise RuntimeError("Run with installed Python -I, not repository imports")
    import orthohmm
    from orthohmm.helpers import SubstitutionMatrix
    from orthohmm.search import engine

    prefix = Path(sys.prefix).resolve()
    for module in (orthohmm, engine):
        if not Path(module.__file__).resolve().is_relative_to(prefix):
            raise RuntimeError("OrthoHMM is not installed inside this interpreter prefix")
    before, ids = input_inventory(directory)
    output.mkdir(parents=True, exist_ok=False)
    result = engine.execute_builtin_search(sorted(ids), str(directory), str(output), 4,
                                          SubstitutionMatrix.blosum62, **SETTINGS)
    counts = write_hits(result, ids, output / "hits.tsv")
    after, after_ids = input_inventory(directory)
    if before != after or ids != after_ids:
        raise ValueError("Input changed during search")
    receipt = dict(status="native_search_completed_pending_independent_readback",
                   inputs=before, hits=fingerprint(output / "hits.tsv"), pairs=counts,
                   source=fingerprint(__file__), engine=fingerprint(engine.__file__),
                   executable=sys.executable, prefix=str(prefix), isolated=True,
                   settings=dict(SETTINGS, cpus=4, matrix="BLOSUM62"),
                   includes_same_species_and_self_hits=True, truth_labels_loaded=False)
    with (output / "receipt.json").open("x") as stream:
        json.dump(receipt, stream, indent=2, sort_keys=True)
        stream.write("\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.input.resolve(), args.output.absolute())
