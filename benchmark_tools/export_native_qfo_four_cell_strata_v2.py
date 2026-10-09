"""Recover the fixed-stratum export using an explicit legacy metadata adapter."""

import argparse
import copy
import json
from pathlib import Path

from benchmark_tools import export_native_qfo_four_cell_strata as legacy


AMENDMENT = "NATIVE_QFO_FOUR_CELL_COMPATIBILITY_AMENDMENT_20261009.md"
AMENDMENT_SHA = "0edde329061ad42a9a3fe14d0e2a55373ffb285fa3e17a6715a72b4cac1b0fd4"
COMPONENT_SHA = "024092bb653947ca2e6d1498498fe51273cd3718795abb0806367977073eac10"


def compatibility_view(docs):
    legacy.family_values(docs)
    view = copy.deepcopy(docs)
    changes = []
    for key in ("fixed", "distance"):
        original = docs[key]
        groups = original["bins"]
        expected = ({(s, c, n) for s in groups for c in legacy.CELLS[:3] for n in groups[s]}
                    if key == "fixed" else
                    {("model_distance", c, n) for c in legacy.CELLS[:3] for n in groups})
        seen = set()
        for row, adapted in zip(original["rows"], view[key]["rows"]):
            suite = row["suite"] if key == "fixed" else "model_distance"
            identity = suite, row["cell"], row["stratum"]
            legacy.require(identity in expected and identity not in seen, "Unexpected legacy row identity")
            seen.add(identity)
            normalized = "resolved_native_pairs" if row["cell"] == legacy.CELLS[1] else "group_clique"
            retained = "native_pair" if key == "distance" and row["cell"] == legacy.CELLS[1] else normalized
            legacy.require(row["prediction_semantics"] == retained, "Unexpected legacy prediction vocabulary")
            if retained != normalized:
                changes.append(dict(input=key, suite=suite, cell=row["cell"], stratum=row["stratum"],
                                    retained=retained, normalized=normalized))
                adapted["prediction_semantics"] = normalized
        legacy.require(seen == expected, "Incomplete legacy row inventory")
    legacy.require(len(changes) == 3, "Wrong compatibility mapping count")
    mapping = dict(scope="in_memory_metadata_only", rows=changes, row_count=3,
                   original_input_bytes_changed=False, numerical_values_changed=False)
    return view, mapping


def project(docs):
    view, mapping = compatibility_view(docs)
    return dict(legacy.project(view), compatibility=mapping)


def export(repo, output):
    output = Path(output)
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    docs, refs, checked, protocol = legacy.prepare(repo)
    component = legacy.record(legacy.__file__)
    legacy.require(component["sha256"] == COMPONENT_SHA, "Changed reused exporter")
    amendment = legacy.record(Path(repo) / "benchmark_tools/results" / AMENDMENT)
    legacy.require(amendment["sha256"] == AMENDMENT_SHA, "Changed compatibility amendment")
    checked = [*checked, component, amendment]
    result = dict(schema="native_qfo_four_cell_strata_v2", source=legacy.record(__file__),
                  components=[component], inputs=refs, protocol=protocol, amendment=amendment,
                  checked_inputs=checked, limitations=legacy.LIMITATIONS, **legacy.SCOPE, **project(docs))
    for ref in [*checked, result["source"]]:
        legacy.require(legacy.record(ref["path"]) == ref, "Direct evidence changed before writing")
    output.mkdir(parents=True)
    legacy.write_tsv(output / "scores.tsv", result["rows"], legacy.SCORE_FIELDS)
    legacy.write_tsv(output / "differences.tsv", result["differences"], legacy.DIFF_FIELDS)
    flat = [dict(cell=r["cell"], family=r["family"], **r["counts_without_prior"],
                 **{m:r[m] for m in legacy.METRICS}) for r in result["family_rows"]]
    legacy.write_tsv(output / "family_counts.tsv", flat, legacy.FAMILY_FIELDS)
    (output / "TABLE.md").write_text(legacy.table(result))
    result["outputs"] = [legacy.record(output / name) for name in
                         ("scores.tsv", "differences.tsv", "family_counts.tsv", "TABLE.md")]
    for ref in [*checked, result["source"]]:
        legacy.require(legacy.record(ref["path"]) == ref, "Direct evidence changed while writing")
    with (output / "report.json").open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    result = export(args.repo, args.output)
    print(json.dumps({key:len(result[key]) for key in ("family_rows", "rows", "differences")}))
