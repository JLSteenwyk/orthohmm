"""Compare retained fixed-release downloads with Ensembl BSD CHECKSUMS."""

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import time
from urllib.parse import urlsplit, urljoin
from urllib.request import urlopen
from urllib.error import HTTPError, URLError

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


def parse_checksums(text):
    rows = {}
    for line in text.splitlines():
        if not line.strip():
            continue
        fields = line.split()
        if len(fields) != 3 or not all(v.isdecimal() for v in fields[:2]):
            raise ValueError("Invalid provider checksum row")
        checksum, blocks = map(int, fields[:2])
        name = fields[2]
        if checksum > 65535 or name in rows or Path(name).name != name:
            raise ValueError("Invalid or duplicate provider filename/checksum")
        rows[name] = dict(bsd_checksum=checksum, blocks_1024=blocks)
    if not rows:
        raise ValueError("Empty checksum inventory")
    return rows


def checksum_url(source):
    parsed = urlsplit(source)
    if (parsed.scheme != "https" or parsed.hostname not in {"ftp.ebi.ac.uk", "ftp.ensembl.org"}
            or parsed.query or parsed.fragment or "current_release" in parsed.path):
        raise ValueError("Require fixed official Ensembl HTTPS source")
    return urljoin(source, "CHECKSUMS")


def reacquire(url, expected, destination):
    """Retain bounded fresh bytes separately; never replace benchmark inputs."""
    started = time.monotonic()
    received = 0
    with destination.open("xb") as handle:
        with urlopen(url, timeout=30) as response:
            if urlsplit(response.url).scheme != "https":
                raise ValueError("Refuse non-HTTPS redirect")
            final_url, headers = response.url, dict(response.headers)
            while True:
                data = response.read(min(1024 * 1024, expected["bytes"] + 1 - received))
                if not data:
                    break
                handle.write(data)
                received += len(data)
                if received > expected["bytes"] or time.monotonic() - started > 120:
                    raise ValueError("Fresh download exceeds size or time bound; partial bytes retained")
    actual = record(destination)
    matched = all(actual[k] == expected[k] for k in ("bytes", "sha256"))
    return dict(file=actual, final_url=final_url, headers=headers, exact_retained_bytes=matched)


def audit(source_path, source_sha, output, download_proteomes=False):
    source_ref = record(source_path)
    if source_ref["sha256"] != source_sha:
        raise ValueError("Retained source audit changed")
    source = json.loads(source_path.read_text())
    output.mkdir(parents=True, exist_ok=False)
    binary = record(Path("/usr/bin/sum"))
    rows = []
    for item in source["inputs"]:
        ref = item["files"]["compressed"]
        check(ref)
        row = dict(code=item["code"], retained=ref, source_url=item["source_url"])
        if item["moving_release_url"]:
            row.update(status="unresolved_moving_release", reason="Do not compare historical bytes to a moving release")
        else:
            url = checksum_url(item["source_url"])
            row.update(checksum_url=url, fetched_at=datetime.now(timezone.utc).isoformat())
            try:
                with urlopen(url, timeout=30) as response:
                    data = response.read(1000001)
                    if len(data) > 1000000 or urlsplit(response.url).scheme != "https":
                        raise ValueError("Invalid or oversized provider response")
                    row.update(final_url=response.url, headers=dict(response.headers))
                path = output / (item["code"] + ".CHECKSUMS")
                path.write_bytes(data)
                row["provider_file"] = record(path)
                entries = parse_checksums(data.decode("ascii"))
                name = Path(urlsplit(item["source_url"]).path).name
                if name not in entries:
                    raise ValueError("Expected source filename absent from provider inventory")
                command = [binary["path"], "-r", ref["path"]]
                result = subprocess.run(command, capture_output=True, text=True, check=True)
                values = result.stdout.split(maxsplit=2)
                observed = dict(bsd_checksum=int(values[0]), blocks_1024=int(values[1]))
                row.update(command=command, command_stdout=result.stdout, expected=entries[name],
                           observed=observed, status="matched" if observed == entries[name] else "mismatch")
                if download_proteomes and row["status"] == "matched":
                    row["fresh_download"] = reacquire(item["source_url"], ref, output / (item["code"] + ".fa.gz"))
                    if not row["fresh_download"]["exact_retained_bytes"]:
                        row["status"] = "fresh_download_mismatch"
            except (HTTPError, URLError, TimeoutError, ValueError, UnicodeError) as error:
                row.update(status="unresolved", error_type=type(error).__name__, error=str(error))
        check(ref)
        rows.append(row)
    check(source_ref)
    check(binary)
    result = dict(source=source_ref, auditor=record(Path(__file__)), checksum_binary=binary,
        rows=rows, download_proteomes=download_proteomes,
        scientific_scores_changed=False, redistribution_cleared=False,
        limitations=["BSD sum is a weak 16-bit checksum with rounded block counts, not cryptographic equality.",
                     "Local SHA256 pins alone are not provider verification; optional fresh downloads are compared separately.",
                     "Provider metadata fetched now does not establish its historical contents or per-tool input use.",
                     "No raw dataset replacement, inference, BUSCO rerun or redistribution-rights decision."])
    (output / "report.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--source-sha256", required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--download-proteomes", action="store_true")
    args = parser.parse_args()
    report = audit(args.source, args.source_sha256, args.output_directory, args.download_proteomes)
    print(json.dumps([(row["code"], row["status"]) for row in report["rows"]]))
