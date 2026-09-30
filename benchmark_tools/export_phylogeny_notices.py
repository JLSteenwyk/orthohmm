"""Export four pinned external-tool notices without code, binaries or clearance."""

import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import tarfile

from benchmark_tools.acquire_publication_fasttree import (
    BASE as FASTTREE_BASE, FILES as FASTTREE_FILES,
    REVISION as FASTTREE_REVISION, identity,
)
from benchmark_tools.build_publication_mafft import ROOT, SHA as MAFFT_SHA, URL
from benchmark_tools.export_bundled_source_notices import PHYLOGENY_SELECTION, verify
from benchmark_tools.run_publication_pipeline import save


def _check(ref):
    current = identity(ref["path"])
    if any(current[key] != ref[key] for key in ("bytes", "sha256")):
        raise ValueError("Changed external-tool evidence: " + ref["path"])


def export(mafft_receipt, mafft_sha, fasttree_receipt, fasttree_sha, output):
    output = Path(output)
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    watched = [identity(mafft_receipt), identity(fasttree_receipt)]
    if [ref["sha256"] for ref in watched] != [mafft_sha, fasttree_sha]:
        raise ValueError("Receipt digest differs from externally supplied pin")
    mafft = json.loads(Path(mafft_receipt).read_bytes())
    fasttree = json.loads(Path(fasttree_receipt).read_bytes())
    if (mafft.get("status") != "core_built_and_installed_phylogeny_fixture_passed"
            or mafft.get("acquisition_url") != URL or mafft.get("extensions_built") is not False
            or mafft["archive"]["sha256"] != MAFFT_SHA):
        raise ValueError("Require pinned core-only MAFFT build evidence")
    if (fasttree.get("status") != "source_notices_acquired_and_installed_binary_matched"
            or fasttree.get("revision") != FASTTREE_REVISION):
        raise ValueError("Require pinned FastTree acquisition evidence")
    archive = mafft["archive"]
    watched.append(archive)
    expected = {}
    for ref in mafft["source_files"]:
        for name in (ROOT + "/license", ROOT + "/license.extensions"):
            if ref["path"].endswith("/" + name):
                if name in expected:
                    raise ValueError("Duplicate MAFFT notice record")
                expected[name] = ref
    if set(expected) != {ROOT + "/license", ROOT + "/license.extensions"}:
        raise ValueError("Missing MAFFT notice record")
    selected = {}
    for name in ("LICENSE", "FastTree.c"):
        matches = [ref for ref in fasttree["acquired_files"]
                   if ref.get("url") == FASTTREE_BASE + name]
        if len(matches) != 1 or matches[0]["sha256"] != FASTTREE_FILES[name]:
            raise ValueError("Missing, duplicate or changed FastTree source identity")
        selected[name] = matches[0]
        watched.append(matches[0])
    for ref in watched:
        _check(ref)

    payloads, seen, total = {}, set(), 0
    # Read archive streams only; never extract archive paths or execute tools.
    with tarfile.open(archive["path"], "r|*") as stream:
        for member in stream:
            name = member.name.rstrip("/")
            path = PurePosixPath(name)
            total += member.size
            if (not path.parts or path.is_absolute() or ".." in path.parts
                    or path.parts[0] != ROOT or name in seen
                    or not (member.isfile() or member.isdir())
                    or member.size < 0 or len(seen) >= 2_000 or total > 20_000_000):
                raise ValueError("Unsafe, duplicate or oversized MAFFT archive")
            seen.add(name)
            if name in expected:
                if not member.isfile() or not 0 < member.size <= 1_000_000:
                    raise ValueError("Invalid MAFFT notice member")
                payload = stream.extractfile(member).read(member.size + 1)
                ref = expected[name]
                if len(payload) != ref["bytes"] or hashlib.sha256(payload).hexdigest() != ref["sha256"]:
                    raise ValueError("MAFFT archive notice differs from build receipt")
                payloads[name] = (payload, archive, {})
    if not set(expected) <= set(payloads):
        raise ValueError("Missing MAFFT archive notice")
    license_payload = Path(selected["LICENSE"]["path"]).read_bytes()
    source = Path(selected["FastTree.c"]["path"]).read_bytes()
    end = source.find(b"*/") + 2
    if not source.startswith(b"/*") or not 2 < end <= 16_384:
        raise ValueError("Missing or oversized FastTree leading notice comment")
    if source[end:end + 1] == b"\n":
        end += 1
    payloads["FastTree/LICENSE"] = (license_payload, selected["LICENSE"], {})
    payloads["FastTree/FastTree.c:leading-comment"] = (
        source[:end], selected["FastTree.c"], {"byte_range": [0, end]})
    for ref in watched:
        _check(ref)

    output.mkdir(parents=True, exist_ok=False)
    files = []
    for member, (payload, ref, extra) in sorted(payloads.items()):
        relative = PHYLOGENY_SELECTION[member]
        target = output / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open("xb") as handle:
            handle.write(payload)
        files.append(dict(member=member, relative_path=relative, source_artifact=ref,
                          bytes=len(payload), sha256=hashlib.sha256(payload).hexdigest(), **extra))
    for ref in watched:
        _check(ref)
    save(output / "SOURCE_NOTICE_INDEX.json", dict(
        scope="selected_external_phylogeny_notices", files=files, inputs=watched,
        source=identity(__file__), redistribution_clearance=False, publication_ready=False,
        limitations=[
            "Four selected exact texts, not complete component attribution or legal clearance.",
            "Extension notice is contextual; optional RNA engines were not built or verified here.",
            "FastTree source-header declaration and separate license file are preserved without reinterpretation.",
            "No source code, binary, source archive or dataset payload exported; no inference repeated.",
            "Build/acquisition receipt bytes are pinned; native fixtures and signed provenance are not rerun.",
            "Verify relocated output using the externally retained index digest, not source paths.",
        ]))
    return verify(output, identity(output / "SOURCE_NOTICE_INDEX.json")["sha256"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mafft-build", type=Path, required=True)
    parser.add_argument("--mafft-build-sha", required=True)
    parser.add_argument("--fasttree-acquisition", type=Path, required=True)
    parser.add_argument("--fasttree-acquisition-sha", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    if args.receipt.exists() or args.receipt.is_symlink():
        raise FileExistsError(args.receipt)
    save(args.receipt, export(args.mafft_build, args.mafft_build_sha,
         args.fasttree_acquisition, args.fasttree_acquisition_sha, args.output.absolute()))
