"""Prespecified WGD diagnostics on admitted group membership, not orthology F1."""

from collections import Counter


def membership(groups, universe):
    result = {}
    for name, genes in groups.items():
        if not isinstance(name, str) or not name or isinstance(genes, str) or not genes:
            raise ValueError("Empty group ID or membership")
        for gene in genes:
            if not isinstance(gene, str) or not gene or gene.strip() != gene:
                raise ValueError("Invalid group member identifier")
            if gene not in universe or gene in result:
                raise ValueError("Foreign or repeated group member: " + gene)
            result[gene] = name
    return result


def score_pairs(cohort, groups, reference, owners):
    native = membership(groups, owners)
    ref = membership(reference, owners)
    sets = {name: set(genes) for name, genes in groups.items()}
    output = []
    seen = set()
    for pair in cohort:
        anchors = pair["orf_pair"]
        key = tuple(anchors)
        if len(anchors) != 2 or anchors[0] >= anchors[1] or key in seen:
            raise ValueError("Noncanonical or duplicate experimental pair")
        seen.add(key)
        ids = [native.get(gene) for gene in anchors]
        members = [sets.get(group, set()) for group in ids]
        if pair["split_eligible"] and any(owners.get(gene) != "Scerevisiae" for gene in anchors):
            raise ValueError("Eligible anchor absent or wrong species")
        state = ("input_excluded" if not pair["split_eligible"] else
                 "incomplete_assignment" if None in ids else "merged" if ids[0] == ids[1] else "separated")
        row = {**pair, "assignment_state": state, "anchor_groups": ids,
               "anchor_group_sizes": list(map(len, members)),
               "separation_rate": None if state == "input_excluded" else int(state == "separated"),
               "supported_separation_rate": None, "mean_non_scer_coverage": None,
               "homolog_support_by_anchor": None, "coverage_numerator": None,
               "coverage_denominator": None, "foreign_pillar_members": None,
               "unmapped_members": None, "union_size": len(set.union(*members)),
               "pillar_native_group_count": None, "unassigned_pillar_members": None}
        if pair["reference_eligible"]:
            pillar = pair["reference_pillar"]
            if not pair["split_eligible"] or pillar not in reference or any(ref.get(g) != pillar for g in anchors):
                raise ValueError("Inconsistent reference eligibility")
            expected = set(reference[pillar])
            if expected != set(pair["available_pillar_members"]):
                raise ValueError("Prepared reference membership differs")
            non_scer = {g for g in expected if owners[g] != "Scerevisiae"}
            union = set.union(*members)
            support = [len(genes & non_scer) for genes in members]
            covered = len(union & non_scer)
            row.update(supported_separation_rate=int(state == "separated" and all(support)),
                       mean_non_scer_coverage=covered / len(non_scer) if non_scer else None,
                       homolog_support_by_anchor=support, coverage_numerator=covered,
                       coverage_denominator=len(non_scer),
                       foreign_pillar_members=sorted(g for g in union if g in ref and ref[g] != pillar),
                       unmapped_members=sorted(union - ref.keys()),
                       pillar_native_group_count=len({native[g] for g in expected if g in native}),
                       unassigned_pillar_members=sorted(expected - native.keys()))
        output.append(row)
    return output


def summarize(rows):
    eligible = [r for r in rows if r["reference_eligible"]]
    result = {"pairs": len(rows), "input_eligible": sum(r["split_eligible"] for r in rows),
              "reference_eligible": len(eligible),
              "assignment_states": dict(Counter(r["assignment_state"] for r in rows)), "endpoints": {}}
    for endpoint in ("separation_rate", "supported_separation_rate", "mean_non_scer_coverage"):
        values = [r[endpoint] for r in eligible if r[endpoint] is not None]
        result["endpoints"][endpoint] = {"pairs": len(values), "mean": sum(values) / len(values) if values else None}
    return result
