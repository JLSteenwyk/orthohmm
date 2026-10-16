"""Distinguish official OrthoFinder 3 groups, later root HOGs and MCL inputs."""

import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmark_tools.audit_ob_orthofinder_provenance import canonical, read_groups, species_ids
from benchmark_tools.orthofinder_mcl_to_orthogroups import iter_mcl_clusters, load_sequence_ids
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.restore_orthofinder_ids import orthofinder_identifier


def padded_hogs(path, mapping, species):
    names = {Path(filename).stem: str(index) for index, filename in species.items()}
    if len(names) != len(species):
        raise ValueError("Ambiguous species column names")
    groups, seen, labels = [], set(), set()
    with path.open() as stream:
        reader = csv.reader(stream, delimiter="\t")
        header = next(reader, [])
        if (header[:3] != ["HOG", "OG", "Gene Tree Parent Clade"]
                or len(header) != len(names) + 3 or set(header[3:]) != set(names)):
            raise ValueError("Unexpected root-HOG columns")
        for row in reader:
            if (len(row) != len(header) or not row[0].startswith("N0.HOG")
                    or row[0] in labels):
                raise ValueError("Malformed or duplicate root-HOG row")
            labels.add(row[0])
            members = []
            for name, cell in zip(header[3:], row[3:]):
                if not cell:
                    continue
                ids = [value.strip() for value in cell.split(",")]
                if any(g not in mapping or g.split("_", 1)[0] != names[name] for g in ids):
                    raise ValueError("Unknown root-HOG gene or wrong species column")
                members.extend(ids)
            if not members or len(members) != len(set(members)) or seen.intersection(members):
                raise ValueError("Empty or overlapping root-HOG membership")
            seen.update(members)
            groups.append([mapping[g] for g in members])
    missing = set(mapping) - seen
    groups.extend([[mapping[g]] for g in missing])
    return groups, len(labels), len(missing)


def audit(root):
    working = root / "WorkingDirectory"
    paths = [working / name for name in ("SequenceIDs.txt", "SpeciesIDs.txt", "N0.ids.tsv",
                                        "clusters_OrthoFinder_I1.2.txt_id_pairs.txt")]
    paths.append(root / "Orthogroups/Orthogroups.txt")
    pins = [record(path) for path in paths]
    mapping = load_sequence_ids(paths[0])
    universe = set(mapping.values())
    if not mapping or len(universe) != len(mapping):
        raise ValueError("Mapped original identifiers must be unique and nonempty")
    species = species_ids(paths[1])
    if {g.split("_", 1)[0] for g in mapping} != {str(i) for i in species}:
        raise ValueError("Sequence mapping and species inventory differ")
    restoration = {orthofinder_identifier(gene): gene for gene in universe}
    if len(restoration) != len(universe):
        raise ValueError("Ambiguous identifier sanitization")
    raw = read_groups(paths[4], True)
    if any(gene not in restoration for group in raw for gene in group):
        raise ValueError("Unknown official group identifier")
    official = canonical(([restoration[g] for g in group] for group in raw), universe)
    restored = sum(restoration[g] != g for group in raw for g in group)
    hogs, hog_rows, singletons = padded_hogs(paths[2], mapping, species)
    latest = canonical(hogs, universe)
    mcl = canonical(([mapping[g] for g in group] for group in iter_mcl_clusters(paths[3])), universe)
    parent = {gene: index for index, group in enumerate(official) for gene in group}
    if any(len({parent[g] for g in group}) != 1 for group in latest):
        raise ValueError("Later root HOG crosses official group boundaries")
    a, b = Counter(official), Counter(latest)
    removed, added = a - b, b - a
    result = dict(genes=len(universe), species=len(species), official_groups=len(official),
        latest_root_hog_rows=hog_rows, latest_root_hog_singleton_padding=singletons,
        latest_root_hog_padded_groups=len(latest), common_groups=sum((a & b).values()),
        official_only_groups=sum(removed.values()), latest_root_hog_only_groups=sum(added.values()),
        genes_in_changed_groups=len({gene for group in removed.keys() | added.keys() for gene in group}),
        identifiers_restored=restored, mcl_groups=len(mcl), official_differs_from_mcl=official != mcl,
        official_equals_latest_root_hog_partition=official == latest,
        latest_root_hogs_refine_official_partition=True, checked=pins)
    for pin in pins:
        check(pin)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", action="append", required=True, help="DATASET=RESULTS_DIRECTORY")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("Preserve prior group-stage audits")
    roots = {}
    for value in args.root:
        name, separator, path = value.partition("=")
        if not separator or not name or not path or name in roots:
            raise ValueError("Require unique named result roots")
        roots[name] = Path(path).resolve()
    helpers = [record(Path(__file__).with_name(name)) for name in (
        "audit_orthofinder_group_stages.py", "audit_ob_orthofinder_provenance.py",
        "orthofinder_mcl_to_orthogroups.py", "restore_orthofinder_ids.py",
        "prepare_ob_candidate_neighborhood.py")]
    rows = {name: audit(root) for name, root in roots.items()}
    for pin in helpers:
        check(pin)
    result = dict(status="orthofinder_group_stage_memberships_checked", rows=rows, helpers=helpers,
        scope="Read-only native partition contrasts. Installed version-specific source and historical run provenance must be checked separately; matching filenames or this partition relation alone do not prove historical execution. No new inference, score selection, current-N0 substitution or changes to retained benchmark scores.",
        defaults_promoted=False, goal_completed=False)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps({name: {k: v for k, v in row.items() if k != "checked"}
                      for name, row in rows.items()}, indent=2))


if __name__ == "__main__":
    main()
