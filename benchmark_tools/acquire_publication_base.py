"""Download exact historical base artifacts into a fresh private directory."""

import argparse
import hashlib
import json
from pathlib import Path
import sys
from urllib.parse import unquote, urlsplit
from urllib.request import build_opener, HTTPRedirectHandler, Request

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmark_tools.install_publication_base import PIP_BYTES, PIP_NAME, PIP_SHA, RECEIPT_SHA
from benchmark_tools.run_integrated_publication_workflow import record, save
from benchmark_tools import stage_base_archives as archives

PROVIDERS = {"repo.anaconda.com", "conda.anaconda.org", "pypi.org", "files.pythonhosted.org"}
PIP_METADATA = "https://pypi.org/pypi/pip/26.2.1/json"
METADATA_LIMIT = 2 * 1024 ** 2


def provider_url(url, hosts=PROVIDERS):
    if not isinstance(url, str) or any(c.isspace() for c in url):
        raise ValueError("Require a provider URL string")
    parsed = urlsplit(url)
    if (parsed.scheme != "https" or parsed.hostname not in hosts
            or parsed.username or parsed.password or parsed.port not in {None, 443}
            or parsed.query or parsed.fragment):
        raise ValueError("Require an allowlisted credential-free HTTPS provider URL")
    return url


class ProviderRedirect(HTTPRedirectHandler):
    def redirect_request(self, request, response, code, message, headers, newurl):
        provider_url(newurl)
        return super().redirect_request(request, response, code, message, headers, newurl)


def network_opener():
    return build_opener(ProviderRedirect())


def download(opener, url, target, limit, timeout, expected=None):
    provider_url(url)
    if target.exists() or target.is_symlink() or target.with_name(target.name + ".partial").exists():
        raise FileExistsError(target)
    partial = target.with_name(target.name + ".partial")
    request = Request(url, headers={"User-Agent": "OrthoHMM-publication-acquisition",
                                   "Accept-Encoding": "identity"})
    size = 0
    with opener.open(request, timeout=timeout) as response, partial.open("xb") as stream:
        final_url = provider_url(response.geturl())
        if response.headers.get("Content-Encoding", "identity") != "identity":
            raise ValueError("Unexpected encoded provider response")
        declared = response.headers.get("Content-Length")
        if declared is not None and (int(declared) < 0 or int(declared) > limit
                or expected is not None and int(declared) != expected["bytes"]):
            raise ValueError("Provider response size differs/exceeds bound")
        while True:
            block = response.read(min(1024 ** 2, limit - size + 1))
            if not block:
                break
            if len(block) > limit - size:
                raise ValueError("Provider response exceeds download bound")
            stream.write(block)
            size += len(block)
    actual = archives.identity(partial)
    if declared is not None and actual["bytes"] != int(declared):
        raise ValueError("Provider response ended before its declared size")
    if expected is not None and any(actual[key] != value for key, value in expected.items()):
        raise ValueError("Downloaded artifact differs from frozen identity")
    partial.rename(target)
    target.chmod(0o644)
    return dict(requested_url=url, final_url=final_url, file=record(target), md5=actual["md5"])


def pip_selection(metadata):
    if metadata["info"]["name"].lower() != "pip" or metadata["info"]["version"] != "26.2.1":
        raise ValueError("Wrong pinned pip release metadata")
    matches = [row for row in metadata["urls"] if row["filename"] == PIP_NAME]
    if len(matches) != 1:
        raise ValueError("Require exactly the pinned pip wheel")
    row = matches[0]
    if (row["packagetype"] != "bdist_wheel" or row["yanked"] is not False
            or type(row["size"]) is not int or row["size"] != PIP_BYTES
            or row["digests"]["sha256"] != PIP_SHA):
        raise ValueError("Provider pip metadata differs from frozen wheel")
    url = provider_url(row["url"], {"files.pythonhosted.org"})
    if unquote(urlsplit(url).path.rsplit("/", 1)[-1]) != PIP_NAME:
        raise ValueError("Wrong pip artifact basename")
    return url


def acquire(args):
    output, receipt = args.output.absolute(), args.receipt.absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    if (output.resolve() != output or type(args.timeout) is not int or args.timeout < 1
            or args.acknowledge_historical_runtime is not True):
        raise ValueError("Require a fresh canonical output, positive timeout and historical acknowledgement")
    if receipt.is_symlink() or not receipt.is_file() or receipt.stat().st_size > 1024 ** 2:
        raise ValueError("Require a bounded regular frozen receipt")
    payload = receipt.read_bytes()
    if hashlib.sha256(payload).hexdigest() != RECEIPT_SHA:
        raise ValueError("Changed frozen reconstruction receipt")
    packages = archives.packages_from_receipt(json.loads(payload))
    if len(packages) != 19:
        raise ValueError("Require all 19 frozen base packages")
    watched = [record(receipt), record(__file__), record(archives.__file__)]
    output.mkdir(parents=True, exist_ok=False)
    (output / "archives").mkdir()
    (output / "bootstrap_wheels").mkdir()
    save(output / "started.json", dict(inputs=watched, packages=packages, attempts=1,
        pip_expected=dict(filename=PIP_NAME, bytes=PIP_BYTES, sha256=PIP_SHA),
        pip_metadata_url=PIP_METADATA, timeout=args.timeout, historical_runtime=True))
    downloads = []
    try:
        opener = network_opener()
        for item in packages:
            package, expected = item["package"], item["expected"]
            target = output / "archives" / package["fn"]
            row = download(opener, package["url"], target, expected["bytes"], args.timeout, expected)
            downloads.append(dict(role="base_archive", package=package, **row))
        metadata_path = output / "pip-release.json"
        meta = download(opener, PIP_METADATA, metadata_path, METADATA_LIMIT, args.timeout)
        downloads.append(dict(role="pip_release_metadata", **meta))
        pip_url = pip_selection(json.loads(metadata_path.read_bytes()))
        row = download(opener, pip_url, output / "bootstrap_wheels" / PIP_NAME,
                       PIP_BYTES, args.timeout, dict(bytes=PIP_BYTES, sha256=PIP_SHA))
        downloads.append(dict(role="pip_wheel", **row))
        for item in watched:
            if record(item["path"]) != item:
                raise ValueError("Acquisition input/source changed")
        for item in downloads:
            if record(item["file"]["path"]) != item["file"]:
                raise ValueError("Acquired payload changed")
        result = dict(status="historical_base_artifacts_acquired", inputs=watched,
            downloads=downloads, base_packages=19, artifact_bytes=sum(p["expected"]["bytes"] for p in packages) + PIP_BYTES,
            attempts=1, retry=False, historical_runtime=True, installation_performed=False,
            native_inference_executed=False, publication_ready=False,
            security_clearance=False, redistribution_clearance=False,
            limitations=["Exact historical artifact acquisition, not a patched installation recommendation.",
                "HTTPS and retained hashes, not provider signatures or complete supply-chain authentication.",
                "Conda bootstrap, scientific wheel/tool sets, raw data and OS runtime remain separately supplied.",
                "No solving, installation, native inference, controlled timing or distribution clearance."])
        save(output / "complete.json", result)
        return result
    except BaseException as error:
        save(output / "failed.json", dict(status="base_artifact_acquisition_failed", inputs=watched,
            type=type(error).__name__, error=str(error), downloads=downloads, attempts=1, retry=False,
            installation_performed=False, publication_ready=False))
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--timeout", type=int, default=60)
    parser.add_argument("--acknowledge-historical-runtime", action="store_true")
    print(json.dumps(acquire(parser.parse_args()), indent=2, sort_keys=True))
