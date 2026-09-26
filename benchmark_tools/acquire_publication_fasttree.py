"""Acquire pinned FastTree 2.2.0 source/notices and identify a retained binary.

Downloads remain local and are never executed. This is an acquisition receipt,
not a reproducible-build, signed-release or redistribution-compliance claim.
"""

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

REVISION = "29c5e62fbcd93230ee325f9c6a17b81f00e3c72a"
BASE = f"https://raw.githubusercontent.com/morgannprice/fasttree/{REVISION}/"
FILES = {
    "FastTree": "55a9d997813aae2208bd4c2081bfa690e0ecdba2d6c491805d8689415c43e38e",
    "FastTree.c": "975202a6b74c9996af871404ff043bb2152edcbda539035662514bc12d1f3431",
    "LICENSE": "3972dc9744f6499f0f9b2dbf76696f2ae7ad8af9b23dde66d6af86c9dfb36986",
    "README.md": "7575f311e6098306079988c9a4933ce21c8c7b8266f655d190b49d85dc539a7b",
    "ChangeLog.txt": "cb9388cc08330a90417571e708cdb0e2418bb59a67863bdf7772196e5b19ffe7",
    "index.html": "6b8e5747a9e959127fde76e52636e35dc09105e6817893606f2ee7f55a856326",
}
MAX_BYTES = 4 * 1024 * 1024


def identity(path):
    path = Path(path)
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"Expected nonsymlink regular file: {path}")
    data = path.read_bytes()
    return dict(path=str(path.absolute()), bytes=len(data),
                sha256=hashlib.sha256(data).hexdigest())


def verify(path, expected):
    result = identity(path)
    if result["sha256"] != expected:
        raise ValueError(f"SHA-256 mismatch: {path}")
    return result


def download(url, destination, expected):
    with urlopen(url, timeout=60) as response:
        final_url = response.geturl()
        if not final_url.startswith("https://"):
            raise ValueError("Non-HTTPS download redirect")
        data = response.read(MAX_BYTES + 1)
    if len(data) > MAX_BYTES:
        raise ValueError("Download exceeds size limit")
    # Preserve mismatching bytes for diagnosis, but never treat them as admitted.
    with destination.open("xb") as stream:
        stream.write(data)
    return dict(**verify(destination, expected), url=url, final_url=final_url)


def run(output, installed):
    output = Path(output).absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    installed_record = verify(installed, FILES["FastTree"])
    output.mkdir(parents=True)
    report = dict(
        status="acquiring", observed_utc=datetime.now(timezone.utc).isoformat(),
        upstream="https://github.com/morgannprice/fasttree", revision=REVISION,
        tag_label="v2.2.0", acquisition="immutable revision URLs with pinned SHA-256",
        installed_binary=installed_record, acquired_files=[], script=identity(__file__),
        executed_downloads=False, benchmarks_rerun=False, installed_binary_modified=False,
        reproducible_build_verified=False, redistribution_compliance_verified=False,
        source_notice="FastTree.c declares GPL version 2 or later; preserve its header.",
        license_file="Upstream LICENSE contains the GNU GPL version 3 text; preserve separately.",
        limitations=[
            "Matching upstream binary bytes does not establish its compiler inputs or a reproducible build.",
            "TLS/Git content identity is not a cryptographically signed maintainer attestation.",
            "This acquisition does not include compiler, OS library or all transitive dependency sources.",
            "No binary or upstream source is newly uploaded by this workflow.",
        ],
    )
    try:
        for name, expected in FILES.items():
            report["acquired_files"].append(download(BASE + name, output / name, expected))
        for name, expected in FILES.items():
            verify(output / name, expected)
        if identity(installed) != installed_record:
            raise ValueError("Installed binary changed during acquisition")
        report["installed_matches_upstream_binary"] = True
        report["status"] = "source_notices_acquired_and_installed_binary_matched"
    except Exception as error:
        report.update(status="acquisition_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        with (output / "receipt.json").open("x") as stream:
            json.dump(report, stream, indent=2, sort_keys=True)
            stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--installed", type=Path, required=True)
    args = parser.parse_args()
    run(args.output, args.installed)
