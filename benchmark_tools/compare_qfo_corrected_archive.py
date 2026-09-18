"""Compare corrected-release canonical FASTAs with the frozen original inputs."""

import argparse
from bisect import bisect_right
from collections import Counter
import gzip
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import sys
import tarfile

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.inventory_swiss_sequences import PREPARED_SHA
from benchmark_tools.resolve_swiss_sequence_aliases import MAPPING_SHA
from benchmark_tools.audit_swiss_missing_input_relations import ALIASES_SHA
from benchmark_tools.qfo_filter_pairs import load_mapping
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check


def missing_by_species(numbers, species, offsets):
    if (len(offsets) != len(species) + 1 or len(set(species)) != len(species)
            or not offsets or offsets[0] != 0
            or any(type(x) is not int for x in offsets)
            or any(a >= b for a, b in zip(offsets, offsets[1:]))):
        raise ValueError("Invalid native species offsets")
    counts = Counter({name: 0 for name in species})
    for number in numbers:
        if type(number) is not int or not 1 <= number <= offsets[-1]:
            raise ValueError("Numeric identity outside species intervals")
        counts[species[bisect_right(offsets, number - 1) - 1]] += 1
    return dict(counts)


def compare(archive, expected, mapping, missing):
    observed, accession_owners, numeric_owners, recovered = {}, {}, {}, {}
    with gzip.open(archive, "rb") as compressed:
        with tarfile.open(fileobj=compressed, mode="r|") as stream:
            for member in stream:
                name = PurePosixPath(member.name).name
                if not name.endswith(".fasta") or name.endswith(("_DNA.fasta", "_additional.fasta")):
                    continue
                if name not in expected or name in observed or not member.isfile():
                    raise ValueError("Unexpected, duplicate or nonregular canonical member")
                with stream.extractfile(member) as handle:
                    data = handle.read()
                if len(data) != member.size:
                    raise ValueError("Truncated member")
                digest = hashlib.sha256(data).hexdigest()
                genes, mapped, unmapped = 0, 0, []
                for entry in SeqIO.parse(io.StringIO(data.decode("utf-8")), "fasta"):
                    parts = entry.id.split("|")
                    if len(parts) != 3 or parts[0] not in {"sp", "tr"} or not all(parts) or not entry.seq:
                        raise ValueError("Invalid canonical FASTA record")
                    accession = parts[1]
                    if accession in accession_owners:
                        raise ValueError("Duplicate canonical accession")
                    accession_owners[accession] = name
                    number = mapping.get(accession)
                    if number is None:
                        unmapped.append(accession)
                    else:
                        if type(number) is not int or number <= 0 or number in numeric_owners:
                            raise ValueError("Invalid or nonunique mapped numeric identity")
                        numeric_owners[number] = accession
                        mapped += 1
                    genes += 1
                    if accession in missing:
                        recovered[accession] = {"file": name, "numeric_protein_id": number,
                                                "length": len(entry.seq), "description": entry.description,
                                                "sequence_sha256": hashlib.sha256(str(entry.seq).encode("ascii")).hexdigest()}
                observed[name] = {"archive_member": member.name, "bytes": len(data), "sha256": digest,
                                  "identical_to_original": (len(data), digest) == (expected[name]["bytes"], expected[name]["sha256"]),
                                  "sequences": genes, "mapped_accessions": mapped, "unmapped_accessions": sorted(unmapped)}
        for _ in iter(lambda: compressed.read(1024 * 1024), b""):
            pass
    if set(observed) != set(expected):
        raise ValueError("Missing canonical proteomes")
    mapped_universe = {number for number in mapping.values() if type(number) is int and number > 0}
    return {"canonical_files": observed, "canonical_proteomes": len(observed),
            "changed_canonical_files": sorted(name for name, value in observed.items() if not value["identical_to_original"]),
            "sequences": len(accession_owners), "unique_mapped_numeric_ids": len(numeric_owners),
            "unmapped_accessions": sum(len(value["unmapped_accessions"]) for value in observed.values()),
            "mapping_numeric_identity_count": len(mapped_universe),
            "mapping_numeric_ids_without_canonical_accession": sorted(mapped_universe - numeric_owners.keys()),
            "recovered_missing_reference_accessions": recovered, "remaining_missing_accessions": sorted(missing - recovered.keys()),
            "gzip_read_to_eof": True}


def audit(archive, prepared_path, mapping_path, alias_path):
    sources = [record(p) for p in (archive, prepared_path, mapping_path, alias_path)]
    if [r["sha256"] for r in sources[1:]] != [PREPARED_SHA, MAPPING_SHA, ALIASES_SHA]:
        raise ValueError("Changed frozen manifests")
    prepared, aliases = json.loads(prepared_path.read_text()), json.loads(alias_path.read_text())
    expected = {Path(r["path"]).name: r for r in prepared["input_fastas"]}
    if len(expected) != 78 or aliases["fasta_inputs"] != prepared["input_fastas"]:
        raise ValueError("Unexpected original input inventory")
    missing = set(aliases["summary"]["missing_genes"])
    with gzip.open(mapping_path, "rt") as stream:
        mapping_data = json.load(stream)
    result = compare(archive, expected, mapping_data["mapping"], missing)
    result["missing_numeric_ids_by_species"] = missing_by_species(
        result["mapping_numeric_ids_without_canonical_accession"], mapping_data["species"], mapping_data["Goff"])
    for identity in sources:
        check(identity)
    return {"status": "qfo_canonical_archive_compared", **result, "inputs": sources,
            "source": record(__file__), "limitations": [
                "Comparison is to frozen original file hashes; no original input or active run changed.",
                "Numeric identity coverage does not independently verify all scorer sequence bytes or annotation versions.",
                "Only canonical FASTA scope compared across releases, not every sidecar or additional file.",
                "No orthology inference, corrected accuracy estimate or retrospective score substitution performed."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("archive", "prepared", "mapping", "aliases", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.archive, args.prepared, args.mapping, args.aliases)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
