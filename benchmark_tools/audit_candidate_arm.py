"""Check candidate partition/seed/merge consistency without reference labels."""

import csv
import json
from pathlib import Path

from benchmark_tools.audit_historical_profile_ablation import read_partition
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


def audit(arm, seed_record, directory, universe, expanded):
    if arm["seed_partition"] != seed_record or arm["candidate_expansion"] is not expanded:
        raise ValueError("Candidate seed or expansion identity differs")
    working = directory / "orthohmm_working_res"
    partition = working / "orthohmm_edges_clustered.txt"
    if arm["candidate_partition"] != record(partition):
        raise ValueError("Candidate partition changed")
    observed = [record(p) for p in sorted(directory.rglob("*")) if p.is_file()]
    if arm["output_files"] != observed or seed_record in observed:
        raise ValueError("Candidate output inventory differs")
    check(seed_record)
    seeds = read_partition(Path(seed_record["path"]), universe)
    candidates = read_partition(partition, universe)
    seed_by_gene = {g: i for i, group in enumerate(seeds) for g in group}
    candidate_by_gene = {g: i for i, group in enumerate(candidates) for g in group}
    if any(len({candidate_by_gene[g] for g in group}) != 1 for group in seeds):
        raise ValueError("Candidate expansion split a seed family")
    expected_names = {"orthohmm_edges_clustered.txt"}
    merges = 0
    if not expanded:
        if (any(arm["candidate_partition"][k] != seed_record[k] for k in ("sha256", "bytes"))
                or "expansion" in arm or "membership_constraints" in arm):
            raise ValueError("Expansion-off arm differs from exact seed copy")
    else:
        expected_names.update(("phylogeny_candidate_superfamilies.txt", "phylogeny_candidate_seeds.tsv",
                               "phylogeny_candidate_merges.json"))
        details = arm["expansion"]
        if details["profile"] != "satellite_v2" or details["membership_policy"] != "high_confidence_pair":
            raise ValueError("Wrong frozen candidate profile")
        checkpoint = working / "phylogeny_candidate_superfamilies.txt"
        sidecar = working / "phylogeny_candidate_seeds.tsv"
        trace_path = working / "phylogeny_candidate_merges.json"
        for key, path in (("candidate_checkpoint", checkpoint), ("seed_sidecar", sidecar), ("merge_trace_sidecar", trace_path)):
            if details[key] != str(path):
                raise ValueError("Candidate sidecar path differs")
        if arm["membership_constraints"] != record(trace_path):
            raise ValueError("Candidate merge trace changed")
        if checkpoint.read_bytes() != partition.read_bytes():
            raise ValueError("Candidate superfamily copy differs")
        expected_rows = [{"candidate_family": f"Family{i:07d}",
                          "seed_families": ",".join(f"Seed{s:07d}" for s in sorted({seed_by_gene[g] for g in group}))}
                         for i, group in enumerate(candidates)]
        with sidecar.open() as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            if reader.fieldnames != ["candidate_family", "seed_families"] or list(reader) != expected_rows:
                raise ValueError("Candidate seed-family sidecar differs")
        trace = json.loads(trace_path.read_text())
        if not isinstance(trace, list):
            raise ValueError("Merge trace must be a list")
        parents = list(range(len(seeds)))
        def find(i):
            while parents[i] != i:
                parents[i] = parents[parents[i]]
                i = parents[i]
            return i
        for entry in trace:
            sides = []
            for prefix in ("source", "target"):
                genes = entry[prefix + "_genes"]
                if (not isinstance(genes, list) or not genes or any(not isinstance(g, str) for g in genes)
                        or len(set(genes)) != len(genes) or not set(genes) <= universe):
                    raise ValueError("Invalid merge gene list")
                ids = {seed_by_gene[g] for g in genes}
                if set(genes) != set().union(*(seeds[i] for i in ids)):
                    raise ValueError("Merge trace splits an original seed")
                if entry[prefix + "_size"] != len(genes) or entry[prefix + "_seed_families"] != len(ids):
                    raise ValueError("Merge trace size/seed counts differ")
                roots = {find(i) for i in ids}
                if len(roots) != 1:
                    raise ValueError("Merge side was not connected by preceding trace")
                sides.append((set(genes), next(iter(roots))))
            if sides[0][0] & sides[1][0] or sides[0][1] == sides[1][1]:
                raise ValueError("Overlapping or redundant merge")
            if len({candidate_by_gene[g] for g in sides[0][0] | sides[1][0]}) != 1:
                raise ValueError("Merge leaves final candidate")
            parents[sides[0][1]] = sides[1][1]
        reconstructed = {}
        for i in range(len(seeds)):
            reconstructed.setdefault(find(i), set()).add(i)
        expected_components = {frozenset(seed_by_gene[g] for g in group) for group in candidates}
        if {frozenset(ids) for ids in reconstructed.values()} != expected_components:
            raise ValueError("Merge trace does not reproduce candidate partition")
        merges = len(trace)
        if (type(details["merges"]) is not int or details["merges"] != merges
                or details["seed_families"] != len(seeds) or details["candidate_families"] != len(candidates)
                or len(seeds) - len(candidates) != merges):
            raise ValueError("Candidate summary counts differ")
    if {str(p.relative_to(working)) for p in directory.rglob("*") if p.is_file()} != expected_names:
        raise ValueError("Unexpected candidate files")
    for item in [seed_record, *observed]:
        check(item)
    return {"status": "candidate_arm_content_verified", "genes": len(universe), "seed_families": len(seeds),
            "candidate_families": len(candidates), "merges": merges, "accuracy_evaluated": False,
            "checked_records": [seed_record, *observed],
            "limitations": ["Checks recorded merge/partition consistency, not independently recomputed search support.",
                            "Parent execution, frozen profile settings and runtime require separate admission."]}
