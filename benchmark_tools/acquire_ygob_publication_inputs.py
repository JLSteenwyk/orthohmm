"""Acquire exact YGOB v7 inputs from upstream without redistributing them."""

import argparse
import hashlib
import json
from pathlib import Path
import time
from urllib.request import urlopen

BASE = "http://ygob.ucd.ie/data/v7-Aug2012/"
FILES = {
    "AA.fsa": (68932999, "f4b8389493b567a688650d6b6c41003a9346649e150efabdf7c8b6f15ac8f356"),
    "Pillars.tab": (2778211, "162c1919b762f34006d0be8ea00b090531ea6c67be43d63b6bf683c4bc48212b"),
    "README": (4958, "162fb37cf0f7a50b44fae47f0d2c2f9b126b728660b65d9bee8d3175af8317f6"),
}


def acquire(output):
    output.mkdir(parents=True, exist_ok=False)
    report = dict(status="acquisition_started", files=[], source_url=BASE,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        started_unix_ns=time.time_ns(), redistribution_cleared=False,
        scientific_results_reproduced=False,
        limitations=["HTTP transport is unauthenticated; SHA-256 binds the previously retained input bytes.",
            "Acquisition is not permission to redistribute; raw files remain excluded from the publication bundle.",
            "Preparation, inference and scoring are separate steps; failures are not silently retried."])
    try:
        for name, (size, digest) in FILES.items():
            partial = output / (name + ".partial")
            item = dict(name=name, url=BASE+name, expected_bytes=size, expected_sha256=digest,
                        status="started", downloaded_bytes=0)
            report["files"].append(item)
            checksum = hashlib.sha256()
            with urlopen(item["url"], timeout=30) as response, partial.open("xb") as stream:
                item.update(final_url=response.geturl(), http_status=response.status)
                if response.status != 200:
                    raise ValueError("Unexpected upstream HTTP status")
                while True:
                    chunk = response.read(min(1024**2, size-item["downloaded_bytes"]+1))
                    if not chunk:
                        break
                    stream.write(chunk)
                    checksum.update(chunk)
                    item["downloaded_bytes"] += len(chunk)
                    if item["downloaded_bytes"] > size:
                        raise ValueError("Upstream content exceeds frozen size")
            item["sha256"] = checksum.hexdigest()
            if item["downloaded_bytes"] != size or item["sha256"] != digest:
                raise ValueError("Upstream bytes differ from frozen YGOB input")
            partial.rename(output / name)
            item["status"] = "verified"
        report["status"] = "frozen_ygob_inputs_acquired"
    except BaseException as error:
        report.update(status="acquisition_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        report["finished_unix_ns"] = time.time_ns()
        with (output / "acquisition.json").open("x") as handle:
            json.dump(report, handle, indent=2, sort_keys=True)
            handle.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    acquire(parser.parse_args().output.absolute())
