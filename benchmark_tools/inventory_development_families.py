"""Inventory explicit reference-family evaluations in an immutable result snapshot."""

import argparse
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess


SNAPSHOT = "f145fff1bd1e425cdd3630ae7abdbde014fd4be0"
ROOTS = ("benchmark_tools/results", "benchmarks/results")
CATALOGS = {
    "OrthoBench": ("benchmark_tools/results/orthobench_phylogeny_partition_20260902.json",
                   "71ce14f0c0df944cb40028a2f63b667875b28e2ac8397831d15c56a421e5826d"),
    "QfO_SwissTrees": ("benchmark_tools/results/qfo_corrected_factorial_complete_20260919/swiss_counts.json",
                       "c9d8bf02ef6f287c56d073fa61f165c982caf834a94fc6e44e6589f03e46eba6"),
}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def file_record(repo, path):
    path = Path(path).resolve()
    require(path.is_relative_to(repo), "Local report outside repository")
    data = path.read_bytes()
    return {"path": str(path.relative_to(repo)), "bytes": len(data), "sha256": digest(data)}


def freeze_local_reports(repo, output, revision=SNAPSHOT):
    repo, output = Path(repo).resolve(), Path(output)
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    snapshot = Snapshot(repo, revision)
    try:
        paths = sorted(p for p in (repo / ROOTS[0]).rglob("*.json")
                       if str(p.relative_to(repo)) not in snapshot.files)
        records = [file_record(repo, p) for p in paths]
    finally:
        snapshot.close()
    require(all(file_record(repo, repo / r["path"]) == r for r in records), "Local file changed during freeze")
    report = {"schema": "local_development_report_manifest_v1", "snapshot_commit": snapshot.revision,
              "local_scan_root": ROOTS[0], "local_reports": records, "raw_files_modified": False}
    with output.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return report


def pointer(parent, field):
    return parent + "/" + str(field).replace("~", "~0").replace("/", "~1")


class Snapshot:
    """Read Git blobs, not mutable worktree files or historical absolute paths."""

    def __init__(self, repo, revision):
        self.repo = Path(repo).resolve()
        self.revision = subprocess.check_output(
            ["git", "rev-parse", "--verify", revision + "^{commit}"], cwd=self.repo, text=True).strip()
        raw = subprocess.check_output(["git", "ls-tree", "-rz", "--full-tree", self.revision,
                                       "--", *ROOTS], cwd=self.repo)
        self.files = {}
        for entry in raw.split(b"\0"):
            if not entry:
                continue
            metadata, name = entry.split(b"\t", 1)
            mode, kind, oid = metadata.decode("ascii").split()
            path = name.decode("utf-8")
            if not path.endswith(".json"):
                continue
            require(kind == "blob" and mode in ("100644", "100755"), "Unsupported result blob")
            self.files[path] = oid
        self.process = subprocess.Popen(["git", "cat-file", "--batch"], cwd=self.repo,
                                        stdin=subprocess.PIPE, stdout=subprocess.PIPE)

    def read(self, path):
        oid = self.files[path]
        self.process.stdin.write((oid + "\n").encode("ascii"))
        self.process.stdin.flush()
        header = self.process.stdout.readline().decode("ascii").split()
        require(len(header) == 3 and header[:2] == [oid, "blob"] and header[2].isdigit(),
                "Unexpected Git blob response")
        size = int(header[2])
        data = self.process.stdout.read(size)
        require(len(data) == size and self.process.stdout.read(1) == b"\n", "Truncated Git blob")
        return json.loads(data), {"path": path, "git_blob": oid, "bytes": size, "sha256": digest(data)}

    def close(self):
        self.process.stdin.close()
        self.process.stdout.close()
        require(self.process.wait(timeout=10) == 0, "Git blob reader failed")


def catalog(ob, swiss):
    names = set(ob["refogs"])
    require(len(names) == 70 and names == set(ob["development"]) | set(ob["validation"])
            and not set(ob["development"]) & set(ob["validation"]), "Wrong original OrthoBench split")
    require(len(ob["development"]) == len(ob["validation"]) == 35, "Wrong partition sizes")
    swiss_names = swiss["families"]
    require(len(swiss_names) == 18 and len(set(swiss_names)) == 18
            and all(isinstance(n, str) and n for n in swiss_names), "Wrong SwissTrees family catalog")
    aliases = {n: ("OrthoBench", n) for n in names}
    aliases.update({n.removesuffix(".txt"): ("OrthoBench", n) for n in names})
    require(not set(aliases) & set(swiss_names), "Reference namespace collision")
    aliases.update({n: ("QfO_SwissTrees", n) for n in swiss_names})
    return aliases, {n: "development" if n in ob["development"] else "validation" for n in names}


def finite_count(value):
    require(type(value) in (int, float) and math.isfinite(value) and value >= 0,
            "Invalid explicit family count")
    return value


def scored_records(records, namespace, aliases):
    require(isinstance(records, list) and bool(records), "Empty/malformed family evaluation")
    converted = []
    for row in records:
        require(isinstance(row, dict), "Malformed scored family record")
        if namespace == "OrthoBench":
            require(row.get("refog") in aliases and aliases[row["refog"]][0] == namespace,
                    "Unknown scored RefOG")
            family = aliases[row["refog"]][1]
            counts = {k: finite_count(row[k]) for k in ("true_positive", "false_positive", "false_negative", "genes")}
        else:
            require(row.get("family") in aliases and aliases[row["family"]][0] == namespace,
                    "Unknown scored SwissTree")
            family = row["family"]
            counts = {k: finite_count(row["counts_without_prior"][k]) for k in ("TP", "FP", "FN", "TN")}
        converted.append({"family": family, "counts": counts})
    require(len({r["family"] for r in converted}) == len(converted), "Duplicate scored family identity")
    return sorted(converted, key=lambda r: r["family"])


def discover(document, aliases):
    blocks, mentions, unresolved = [], [], []

    def visit(node, location):
        if isinstance(node, dict):
            if "refog_records" in node:
                values = node["refog_records"]
                require(isinstance(values, list) and bool(values) and all(isinstance(r, dict) for r in values),
                        "Malformed RefOG record container")
                fields = {"refog", "genes", "true_positive", "false_positive", "false_negative"}
                if all(fields <= r.keys() for r in values):
                    if all(r["refog"] in aliases and aliases[r["refog"]][0] == "OrthoBench" for r in values):
                        records = scored_records(values, "OrthoBench", aliases)
                        blocks.append({"pointer": pointer(location, "refog_records"), "dataset": "OrthoBench", "records": records})
                    else:
                        unresolved.append({"pointer": pointer(location, "refog_records"),
                            "kind": "unmapped_scored_reference_labels", "recorded_labels": [r["refog"] for r in values],
                            "reference_family_mapping_established": False})
                elif all({"refog", "genes", "recovered_pairs", "possible_pairs"} <= r.keys() for r in values):
                    # The historical diagnostic constructs these labels by subset ordinal, not filename.
                    unresolved.append({"pointer": pointer(location, "refog_records"),
                        "kind": "ordinal_candidate_diagnostic", "recorded_labels": [r["refog"] for r in values],
                        "reference_family_mapping_established": False})
                else:
                    raise ValueError("Unrecognized explicit RefOG record schema")
            for name, child in node.items():
                where = pointer(location, name)
                if name in aliases and isinstance(child, (dict, list)):
                    dataset, family = aliases[name]
                    mentions.append({"pointer": where, "dataset": dataset, "family": family,
                                     "kind": "named_family_map"})
                if name in ("families", "refog_names", "development", "validation") and isinstance(child, list):
                    for i, value in enumerate(child):
                        if isinstance(value, str) and value in aliases:
                            dataset, family = aliases[value]
                            mentions.append({"pointer": pointer(where, i), "dataset": dataset, "family": family,
                                             "kind": "declared_family_list"})
                visit(child, where)
        elif isinstance(node, list):
            if node and all(isinstance(v, dict) and "counts_without_prior" in v and "family" in v for v in node):
                records = scored_records(node, "QfO_SwissTrees", aliases)
                blocks.append({"pointer": location, "dataset": "QfO_SwissTrees", "records": records})
            for i, child in enumerate(node):
                visit(child, pointer(location, i))

    visit(document, "")
    return blocks, mentions, unresolved


def collect(repo, revision=SNAPSHOT, local_manifest=None, local_sha=None):
    repo = Path(repo).resolve()
    snapshot = Snapshot(repo, revision)
    files, evaluations, mentions, unresolved, contexts = [], [], [], [], {}
    local_refs, local_pin = [], None
    try:
        if local_manifest is not None:
            local_pin = file_record(repo, local_manifest)
            require(local_pin["sha256"] == local_sha, "Local manifest identity changed")
            plan = json.loads((repo / local_pin["path"]).read_text())
            require(plan["schema"] == "local_development_report_manifest_v1"
                    and plan["snapshot_commit"] == snapshot.revision
                    and plan["local_scan_root"] == ROOTS[0] and plan["raw_files_modified"] is False,
                    "Wrong local report manifest")
            local_refs = plan["local_reports"]
            require(len({r["path"] for r in local_refs}) == len(local_refs)
                    and all(r["path"] not in snapshot.files and Path(r["path"]).is_relative_to(ROOTS[0])
                            and ".." not in Path(r["path"]).parts and r["path"].endswith(".json")
                            for r in local_refs), "Duplicate/nonlocal report identity")
        documents = {}
        for dataset, (path, expected) in CATALOGS.items():
            documents[dataset], ref = snapshot.read(path)
            require(ref["sha256"] == expected, "Family catalog changed")
        aliases, roles = catalog(documents["OrthoBench"], documents["QfO_SwissTrees"])
        def blobs():
            for path in sorted(snapshot.files):
                document, ref = snapshot.read(path)
                yield document, {**ref, "origin": "git_snapshot"}
            for frozen in sorted(local_refs, key=lambda r: r["path"]):
                require(file_record(repo, repo / frozen["path"]) == frozen, "Local report checksum mismatch")
                data = (repo / frozen["path"]).read_bytes()
                require(len(data) == frozen["bytes"] and digest(data) == frozen["sha256"], "Local report changed during read")
                yield json.loads(data), {**frozen, "origin": "frozen_local_only", "git_blob": None}

        for document, ref in blobs():
            path = ref["path"]
            require(isinstance(document, (dict, list)), "Unexpected result document")
            blocks, refs, opaque = discover(document, aliases)
            ref["scored_family_blocks"] = len(blocks)
            ref["explicit_family_mentions"] = len(refs)
            ref["unresolved_family_containers"] = len(opaque)
            files.append(ref)
            if not blocks and not refs and not opaque:
                continue
            context = {k: document[k] for k in ("status", "description", "generated_at", "input_release", "partition",
                                               "command", "source", "git") if isinstance(document, dict) and k in document}
            contexts[path] = context
            for block in blocks:
                payload = json.dumps(block.pop("records"), sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
                records = json.loads(payload)
                evaluations.append({**block, "file": path, "families": [r["family"] for r in records],
                    "count_payload_sha256": digest(payload),
                    "reported_partition": context.get("partition"), "evidence_kind": "explicit_family_counts"})
            mentions.extend({"file": path, **row} for row in refs)
            unresolved.extend({"file": path, **row} for row in opaque)
        require(all(file_record(repo, repo / r["path"]) == r for r in local_refs), "Local report changed during scan")
        if local_pin is not None:
            require(file_record(repo, repo / local_pin["path"]) == local_pin, "Local manifest changed during scan")
    finally:
        snapshot.close()
    family_rows = []
    for dataset in ("OrthoBench", "QfO_SwissTrees"):
        names = sorted({family for namespace, family in aliases.values() if namespace == dataset})
        for family in names:
            relevant = [e for e in evaluations if e["dataset"] == dataset and family in e["families"]]
            roles_seen = Counter(e["reported_partition"] or "not_declared" for e in relevant)
            family_rows.append({"dataset": dataset, "family": family,
                "original_partition": roles.get(family), "publication_exposure": "development_exposed",
                "scored_evidence_blocks": len(relevant), "scored_evidence_files": len({e["file"] for e in relevant}),
                "reported_partition_block_counts": dict(sorted(roles_seen.items())),
                "independent_native_experiments": None, "causal_tuning_influence_established": False})
    summary = {}
    for dataset in ("OrthoBench", "QfO_SwissTrees"):
        blocks = [e for e in evaluations if e["dataset"] == dataset]
        summary[dataset] = {"reference_families": sum(r["dataset"] == dataset for r in family_rows),
            "scored_blocks": len(blocks), "scored_files": len({e["file"] for e in blocks}),
            "family_block_associations": sum(len(e["families"]) for e in blocks),
            "distinct_count_vectors": len({e["count_payload_sha256"] for e in blocks}),
            "families_with_scored_evidence": len({f for e in blocks for f in e["families"]})}
    return {"schema": "development_family_evidence_inventory_v1", "status": "retained_family_evaluations_inventoried_partial",
            "snapshot_commit": snapshot.revision, "scan_roots": list(ROOTS), "json_files_scanned": len(files),
            "snapshot_json_bytes": sum(f["bytes"] for f in files if f["origin"] == "git_snapshot"),
            "local_json_bytes": sum(f["bytes"] for f in files if f["origin"] == "frozen_local_only"),
            "local_report_manifest": local_pin, "files": files, "file_contexts": contexts,
            "evaluations": evaluations, "family_mentions": mentions, "unresolved_family_containers": unresolved,
            "family_rows": family_rows, "summary": summary,
            "source": {"path": "benchmark_tools/inventory_development_families.py",
                       "sha256": digest(Path(__file__).read_bytes()), "bytes": Path(__file__).stat().st_size},
            "native_or_scoring_repeated": False, "causal_tuning_history_complete": False,
            "family_disjoint_validation_established": False, "publication_ready": False,
            "limitations": [
                "Exhaustive scan of tracked JSON result blobs at one immutable Git snapshot plus the separately hash-frozen local-only JSON reports under benchmark_tools/results when supplied. Not all historical commits, local native-result directories, Markdown, raw scores, human inspection history or unrecorded experiments.",
                "Explicit RefOG and SwissTrees count schemas establish evaluated reference families in retained reports, not that every report is an independent native experiment or causally influenced tuning. Repeated vectors are retained and counted separately as evidence blocks; equal vectors do not prove identical execution.",
                "Only the exact 70 RefOG and 18 retained SwissTrees identities are canonicalized. Named/list mentions remain separate from explicit score blocks; plans and family counts without identities are not promoted to evaluations.",
                "OrthoBench's original development/validation assignment is preserved as chronology, not a new independence certificate. Repeatedly evaluated validation families are publication-development-exposed too.",
                "Historical candidate diagnostics synthesize RefOG labels from subset ordinal rather than original filename. Their containers/recorded labels are retained separately as unresolved, not misassigned to canonical families or promoted to official scores.",
                "Scored records outside the exact reference catalog, including retained synthetic runtime fixtures, are inventoried as unmapped labels rather than discarded or relabelled as biological RefOGs.",
                "Original TreeFam-A family reconstruction remains unavailable; pooled relations and newer/derived tree labels cannot replace it. VGNC, GO, EC, FAS, BUSCO, simulations and YGOB use other units; these schemas do not invent a cross-resource family map.",
                "Frozen YGOB provides the separately documented novel-taxon route, not family-disjoint validation. Its admitted development-homology screen and frozen-before-score chronology remain required independent evidence, not recomputed here.",
                "Snapshot source/command text is retained provenance, not launch-input consumption or full runtime attestation. Existing scores, methods, endpoints, timing panel, main PDF and rc3 archive are unchanged.",
            ]}


def association_rows(report):
    return [{"dataset": e["dataset"], "family": family, "evidence_file": e["file"],
             "json_pointer": e["pointer"], "reported_partition": e["reported_partition"],
             "count_payload_sha256": e["count_payload_sha256"], "evidence_kind": e["evidence_kind"],
             "independent_native_experiment": "unestablished"}
            for e in report["evaluations"] for family in e["families"]]


def render(report):
    lines = ["# Retained Development-Family Evidence Inventory", "",
             "Snapshot: `%s`. %d JSON reports; %d committed bytes plus %d separately frozen local bytes." % (
                 report["snapshot_commit"], report["json_files_scanned"], report["snapshot_json_bytes"], report["local_json_bytes"]),
             "Evidence blocks are not independent experiments or proof of causal tuning influence.", "",
             "| Dataset | Families | Scored Blocks | Scored Files | Family-Block Associations | Distinct Count Vectors |",
             "| --- | ---: | ---: | ---: | ---: | ---: |"]
    for name, values in report["summary"].items():
        lines.append("| %s | %d | %d | %d | %d | %d |" % (name, values["reference_families"],
            values["scored_blocks"], values["scored_files"], values["family_block_associations"], values["distinct_count_vectors"]))
    validation = [r for r in report["family_rows"] if r["original_partition"] == "validation"]
    lines += ["", "%d originally validation-labelled RefOGs have explicit score evidence; all remain development-exposed for this publication." %
              sum(r["scored_evidence_blocks"] > 0 for r in validation), "", "## Limits", "",
              *["- " + text for text in report["limitations"]]]
    return "\n".join(lines) + "\n"


def write(repo, output, revision=SNAPSHOT, local_manifest=None, local_sha=None):
    output = Path(output)
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    report = collect(repo, revision, local_manifest, local_sha)
    output.mkdir(parents=True)
    (output / "inventory.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    (output / "inventory.md").write_text(render(report))
    rows = association_rows(report)
    require(bool(rows), "No explicit evaluated reference families")
    with (output / "family_evidence.tsv").open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({k: "NA" if v is None else v for k, v in row.items()} for row in rows)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--revision", default=SNAPSHOT)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--freeze-local", action="store_true")
    parser.add_argument("--local-manifest", type=Path)
    parser.add_argument("--local-manifest-sha256")
    args = parser.parse_args()
    if args.freeze_local:
        result = freeze_local_reports(args.repo, args.output, args.revision)
        print(json.dumps({"frozen_local_reports": len(result["local_reports"])}))
    else:
        result = write(args.repo, args.output, args.revision, args.local_manifest, args.local_manifest_sha256)
        print(json.dumps({"files": result["json_files_scanned"], "summary": result["summary"]}))
