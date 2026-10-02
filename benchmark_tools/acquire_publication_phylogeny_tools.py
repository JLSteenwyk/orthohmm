"""Acquire frozen MAFFT/FastTree artifacts without installed-tool prerequisites.

Only download and byte verification occur. Source compilation, executable
permission, tool installation and inference are deliberately separate steps.
"""

import argparse
from datetime import datetime, timezone
from pathlib import Path
import sys
from urllib.request import build_opener, HTTPRedirectHandler

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmark_tools import acquire_publication_base as base
from benchmark_tools import acquire_publication_fasttree as fasttree
from benchmark_tools import build_publication_mafft as mafft

PROVIDERS = frozenset({"mafft.cbrc.jp", "raw.githubusercontent.com"})
FASTTREE_SIZES = {
    "FastTree": 1496928,
    "FastTree.c": 395674,
    "LICENSE": 35149,
    "README.md": 662,
    "ChangeLog.txt": 17271,
    "index.html": 49797,
}


class ToolRedirect(HTTPRedirectHandler):
    def redirect_request(self, request, response, code, message, headers, newurl):
        base.provider_url(newurl, PROVIDERS)
        return super().redirect_request(request, response, code, message, headers, newurl)


def network_opener():
    return build_opener(ToolRedirect())


def artifacts():
    result = [dict(role="mafft_source", relative="mafft/mafft-7.525-with-extensions-src.tgz",
                   url=mafft.URL, bytes=758305, sha256=mafft.SHA)]
    result.extend(dict(role="fasttree_binary" if name == "FastTree" else "fasttree_source_notice",
                       relative="fasttree/" + name, url=fasttree.BASE + name,
                       bytes=size, sha256=fasttree.FILES[name])
                  for name, size in FASTTREE_SIZES.items())
    return result


def source_records():
    # Imported local modules are part of the controller, not downloaded inputs.
    paths = {Path(__file__).resolve()}
    paths.update(Path(module.__file__).resolve() for name, module in tuple(sys.modules.items())
                 if name.startswith("benchmark_tools.") and getattr(module, "__file__", None))
    return [base.record(path) for path in sorted(paths)]


def acquire(args):
    output = args.output.absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    if (output.resolve() != output or type(args.timeout) is not int or args.timeout < 1
            or args.acknowledge_historical_runtime is not True):
        raise ValueError("Require a fresh canonical output, positive timeout and historical acknowledgement")
    if output.is_relative_to(Path(__file__).resolve().parent.parent):
        raise ValueError("Acquisition output must be outside the source component")
    rows = artifacts()
    watched = source_records()
    for row in rows:
        base.provider_url(row["url"], PROVIDERS)
    output.mkdir(parents=True, exist_ok=False)
    for name in ("mafft", "fasttree"):
        (output / name).mkdir()
    started = dict(status="acquiring_frozen_phylogeny_artifacts", inputs=watched,
                   artifacts=rows, started_utc=datetime.now(timezone.utc).isoformat(),
                   attempts=1, retry=False, timeout_seconds=args.timeout,
                   historical_runtime=True, downloads_executed=False)
    base.save(output / "started.json", started)
    downloads = []
    scope = dict(attempts=1, retry=False, historical_runtime=True,
                 installed_tools_required=False, extraction_performed=False,
                 compilation_performed=False, installation_performed=False,
                 native_code_executed=False, scientific_inference_executed=False,
                 controlled_timing=False, publication_ready=False,
                 security_clearance=False, redistribution_clearance=False)
    try:
        opener = network_opener()
        for row in rows:
            expected = {key: row[key] for key in ("bytes", "sha256")}
            result = base.download(opener, row["url"], output / row["relative"], row["bytes"],
                                   args.timeout, expected, hosts=PROVIDERS)
            downloads.append(dict(role=row["role"], relative=row["relative"], **result))
        for item in [*watched, *(row["file"] for row in downloads)]:
            if base.record(item["path"]) != item:
                raise ValueError("Source/acquired artifact changed during acquisition")
        if any(Path(row["file"]["path"]).stat().st_mode & 0o777 != 0o644 for row in downloads):
            raise ValueError("Acquired artifacts must remain non-executable regular files")
        report = dict(status="frozen_phylogeny_artifacts_acquired", inputs=watched,
                      downloads=downloads, artifact_files=len(downloads),
                      artifact_bytes=sum(row["file"]["bytes"] for row in downloads),
                      finished_utc=datetime.now(timezone.utc).isoformat(), **scope,
                      fasttree_revision=fasttree.REVISION, mafft_version="7.525",
                      fasttree_version="2.2.0",
                      limitations=[
                          "Exact acquisition, not a complete installed native-tool runtime.",
                          "MAFFT still requires compilation, launcher prefix/helper configuration and tool validation.",
                          "The FastTree binary matches retained bytes; it is not a source-reproducibility claim.",
                          "HTTPS plus retained hashes is not a signed maintainer attestation or security audit.",
                          "Compiler, OS libraries, optional extension engines and selected-artifact rights remain open.",
                          "No scientific/default/lock change, benchmark rerun, timing admission or public binary upload.",
                          "Socket-operation timeout is not a whole-acquisition deadline.",
                      ])
        base.save(output / "complete.json", report)
        return report
    except BaseException as error:
        base.save(output / "failed.json", dict(status="phylogeny_artifact_acquisition_failed",
                  inputs=watched, downloads=downloads, error_type=type(error).__name__,
                  error=str(error), **scope))
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--timeout", type=int, default=60)
    parser.add_argument("--acknowledge-historical-runtime", action="store_true")
    args = parser.parse_args()
    result = acquire(args)
    print(base.json.dumps(dict(status=result["status"], complete=base.record(args.output / "complete.json")),
                          sort_keys=True))


if __name__ == "__main__":
    main()
