"""Short-lived worker for isolating native Leiden graph state."""

from __future__ import annotations

import json
import os
from pathlib import Path
import sys

import numpy as np

from .externals import _execute_leiden_in_process
from .helpers import IndexedEdges


def main(argv=None) -> int:
    arguments = list(sys.argv[1:] if argv is None else argv)
    if len(arguments) != 1:
        raise SystemExit("usage: python -m orthohmm.leiden_worker PAYLOAD_DIR")
    payload_directory = Path(arguments[0]).resolve()
    metadata = json.loads((payload_directory / "metadata.json").read_text())
    gene_names = (payload_directory / "gene_names.txt").read_text().splitlines()
    edges = IndexedEdges(
        gene_names=gene_names,
        sources=np.load(
            payload_directory / "sources.npy", mmap_mode="r", allow_pickle=False
        ),
        targets=np.load(
            payload_directory / "targets.npy", mmap_mode="r", allow_pickle=False
        ),
        weights=np.load(
            payload_directory / "weights.npy", mmap_mode="r", allow_pickle=False
        ),
    )
    _execute_leiden_in_process(
        metadata["cpm_resolution"],
        metadata["output_directory"],
        edges=edges,
        include_isolates=metadata["include_isolates"],
        seed=metadata["seed"],
    )
    sys.stdout.flush()
    sys.stderr.flush()
    os._exit(0)


if __name__ == "__main__":
    main()
