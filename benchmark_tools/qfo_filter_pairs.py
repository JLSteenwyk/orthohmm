#!/usr/bin/env python3
"""Filter pair predictions to identifiers recognized by a QfO mapping."""

import argparse
import gzip
import json
import sys
from pathlib import Path


def load_mapping(path):
    with gzip.open(path, "rt") as handle:
        data = json.load(handle)
    return data["mapping"]


def filter_pairs(input_path, output_path, valid_ids):
    total = 0
    retained = 0
    with input_path.open() as input_handle, output_path.open("w") as output_handle:
        for line_number, line in enumerate(input_handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 2:
                raise ValueError(
                    f"Expected two tab-separated IDs at {input_path}:{line_number}"
                )
            total += 1
            if fields[0] in valid_ids and fields[1] in valid_ids:
                output_handle.write(line)
                retained += 1
    return total, retained


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("mapping", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    valid_ids = load_mapping(args.mapping)
    total, retained = filter_pairs(args.input, args.output, valid_ids)
    print(
        f"Retained {retained} of {total} pairs; removed {total - retained} "
        "pairs containing IDs absent from the QfO mapping",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
