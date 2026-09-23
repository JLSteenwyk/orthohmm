"""Retain official linked model trees without treating them as benchmark labels."""

import argparse
from datetime import datetime, timezone
from html.parser import HTMLParser
import json
from pathlib import Path
import re
import time
from urllib.request import urlopen

from benchmark_tools.inventory_swisstree_duplications import FAMILIES
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.run_blast_recovery_batch import save_status


class ModelLinks(HTMLParser):
    def __init__(self, identifier):
        super().__init__(convert_charrefs=True)
        if identifier not in FAMILIES:
            raise ValueError("Unknown source family")
        self.pattern = re.compile(rf"/ST/{identifier}/(?:modeltree\.nhx|{identifier}_treemodel\.phyloxml)(?![\w./])")
        self.paths = set()

    def handle_starttag(self, tag, attrs):
        for key, value in attrs:
            if key.lower() in {"href", "onclick", "src"} and value:
                self.handle_data(value)

    def handle_data(self, data):
        self.paths.update(self.pattern.findall(data))


def model_path(page, identifier):
    parser = ModelLinks(identifier)
    parser.feed(page)
    if len(parser.paths) != 1:
        raise ValueError("Require one uniquely linked official model tree")
    return next(iter(parser.paths))


def acquire(output):
    output.mkdir(parents=True, exist_ok=False)
    report = dict(status="acquiring", acquired_utc=datetime.now(timezone.utc).isoformat(),
        source=record(__file__), families={}, benchmark_labels_admitted=False,
        prediction_statistics_evaluated=False, publication_ready=False)
    save_status(output / "manifest.json", report)
    for identifier, family in FAMILIES.items():
        if family is None:
            continue
        directory = output / identifier
        directory.mkdir()
        item = dict(qfo_family=family, downloads=[])
        report["families"][identifier] = item
        try:
            page_url = f"https://swisstree.sib.swiss/cgi-bin/swisst?page={identifier}"
            time.sleep(1)
            with urlopen(page_url, timeout=60) as response:
                content = response.read()
            path = directory / "page.html"
            path.write_bytes(content)
            item["downloads"].append(dict(url=page_url, file=record(path)))
            relative = model_path(content.decode("utf-8"), identifier)
            tree_url = "https://swisstree.sib.swiss" + relative
            time.sleep(1)
            with urlopen(tree_url, timeout=60) as response:
                content = response.read()
            path = directory / Path(relative).name
            path.write_bytes(content)
            item["downloads"].append(dict(url=tree_url, file=record(path)))
            item.update(status="downloaded_pending_tree_and_mapping_validation", tree_format="phyloxml" if path.suffix == ".phyloxml" else "nhx")
        except Exception as error:
            item.update(status="unavailable", error_type=type(error).__name__, error=str(error))
        save_status(output / "manifest.json", report)
    report["status"] = "acquisition_finished_not_admitted"
    save_status(output / "manifest.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = acquire(args.output.absolute())
    print(json.dumps({key: value["status"] for key, value in result["families"].items()}))
