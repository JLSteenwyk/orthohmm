"""Trace reviewed residue-deletion proteins into admitted final groups."""

import argparse
import hashlib
import json
from pathlib import Path

from benchmark_tools.audit_orthomcl_native_groups import HEADER, MEMBER


def read_checked(path, digest):
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != digest:
        raise ValueError("Input digest mismatch")
    return json.loads(data)


def scan(path, digest, targets):
    if not targets:
        raise ValueError("Empty target set")
    found, selected = set(), []
    hashed = hashlib.sha256()
    with path.open("rb") as stream:
        for number, raw in enumerate(stream, 1):
            hashed.update(raw)
            match = HEADER.fullmatch(raw.decode().strip())
            if not match:
                raise ValueError("Malformed group header")
            members = [MEMBER.fullmatch(token) for token in match[4].split()]
            if not all(members):
                raise ValueError("Malformed group member")
            genes = [m[1] for m in members]
            species = [m[2] for m in members]
            if len(genes) != int(match[2]) or len(set(species)) != int(match[3]):
                raise ValueError("Group header count mismatch")
            hits = targets.intersection(genes)
            if not hits:
                continue
            if found.intersection(hits) or len(genes) != len(set(genes)):
                raise ValueError("Duplicate target membership")
            found.update(hits)
            selected.append(dict(line=number, group="ORTHOMCL" + match[1],
                members=[dict(gene=g, species=s) for g, s in zip(genes, species)],
                targets=sorted(hits),
                incident_cross_species_pairs=sum(
                    species[i] != species[j] and (g in targets or genes[j] in targets)
                    for i, g in enumerate(genes) for j in range(i + 1, len(genes)))))
    if hashed.hexdigest() != digest:
        raise ValueError("Group digest mismatch")
    return dict(groups=selected, absent_targets=sorted(targets - found),
        incident_cross_species_pairs=sum(g["incident_cross_species_pairs"] for g in selected))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("trace", "admission", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    for name in ("trace-sha256", "admission-sha256"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args()
    trace = read_checked(args.trace, args.trace_sha256)
    admission = read_checked(args.admission, args.admission_sha256)
    if trace["status"] != "selected_hit_trace_verified_against_supplied_table_identity":
        raise ValueError("Unverified search trace")
    if admission["status"] != "recovered_orthomcl_native_outputs_admitted":
        raise ValueError("Unadmitted native groups")
    targets = {p["gene"] for p in trace["content"]["proteins"]}
    if len(targets) != 7:
        raise ValueError("Require seven reviewed proteins")
    groups = admission["native_groups"]
    result = dict(status="reviewed_residue_final_group_trace_verified",
        trace_sha256=args.trace_sha256, admission_sha256=args.admission_sha256,
        native_groups=groups, source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        content=scan(Path(groups["path"]), groups["sha256"], targets),
        limitations=["Observed final-group trace, not a residue-preserving search counterfactual.",
                     "Zero incident cross-species pairs does not exclude counterfactual clustering changes."])
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")


if __name__ == "__main__":
    main()
