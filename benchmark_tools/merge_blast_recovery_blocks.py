"""Streaming recovery merge primitives; not a prefix or scientific admission CLI."""

import hashlib
import os
from pathlib import Path
import re


def ordered_blocks(queries, prefix_blocks, replay_blocks, replay_ids, prefix_end):
    """Interleave admitted block inventories in exhaustive original query order."""
    if type(prefix_end) is not int or prefix_end < 0:
        raise ValueError("Invalid prefix boundary")
    if len(replay_ids) != len(set(replay_ids)):
        raise ValueError("Duplicate replay query")
    replay = set(replay_ids)
    prefix_iter, replay_iter = iter(prefix_blocks), iter(replay_blocks)
    retained, new = next(prefix_iter, None), next(replay_iter, None)
    seen, replay_order, end = set(), [], 0
    for gene in queries:
        if not isinstance(gene, str) or not gene or gene in seen:
            raise ValueError("Invalid or duplicated original query")
        seen.add(gene)
        if gene in replay:
            replay_order.append(gene)
            if retained is not None and retained["query"] == gene:
                raise ValueError("Replayed query also present in retained prefix")
            if new is not None and new["query"] == gene:
                yield dict(new, source="replay")
                new = next(replay_iter, None)
        else:
            if retained is None or retained["query"] != gene:
                raise ValueError("Non-replay query lacks ordered retained block")
            if (retained["start"] != end or retained["end"] <= end
                    or retained["end"] > prefix_end or retained.get("final_observed_query") is not False):
                raise ValueError("Discontinuous prefix or final incomplete block included")
            end = retained["end"]
            yield dict(retained, source="prefix")
            retained = next(prefix_iter, None)
    if retained is not None or new is not None or replay_order != replay_ids or end != prefix_end:
        raise ValueError("Unconsumed, reordered or incomplete recovery inventory")


def copy_block(stream, destination, block, digest):
    """Copy one exact query block, refusing malformed rows or changed bytes."""
    start, end, rows = (block[key] for key in ("start", "end", "rows"))
    if (any(type(value) is not int for value in (start, end, rows))
            or start < 0 or end <= start or rows <= 0
            or re.fullmatch(r"[0-9a-f]{64}", block["sha256"]) is None):
        raise ValueError("Invalid byte-range inventory")
    query = block["query"].encode("ascii")
    stream.seek(start)
    position, observed = start, 0
    block_digest = hashlib.sha256()
    while position < end:
        line = stream.readline(min(65536, end-position))
        if (not line or not line.endswith(b"\n") or b"\0" in line
                or len(line.rstrip(b"\n").split(b"\t")) != 12
                or line.split(b"\t", 1)[0] != query):
            raise ValueError("Changed, truncated or misidentified query block")
        destination.write(line)
        block_digest.update(line)
        digest.update(line)
        position += len(line)
        observed += 1
    if observed != rows or block_digest.hexdigest() != block["sha256"]:
        raise ValueError("Query block row count or hash differs")
    return observed, end-start


def merge(blocks, sources, output):
    """Write a new candidate only; callers must separately authorize and admit it.

    Each block has a path key selecting a verified source. Numeric alignment,
    scheduler, no-hit/failure and source provenance gates belong to the caller.
    Failures before rename preserve the partial file. A candidate is never
    an admission marker, including if directory fsync fails after rename.
    """
    output = Path(output)
    output.mkdir(exist_ok=False)
    partial, ready = output / "all.blast.partial", output / "all.blast.candidate"
    streams, identities = {}, {}
    total_rows, total_bytes, count, seen = 0, 0, 0, set()
    digest = hashlib.sha256()
    try:
        for key, path in sources.items():
            path = Path(path)
            stream = path.open("rb")
            streams[key] = stream
            stat = os.fstat(stream.fileno())
            identities[key] = (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)
        with partial.open("xb") as destination:
            for block in blocks:
                if block["query"] in seen or block["path"] not in streams:
                    raise ValueError("Repeated query or unknown source")
                seen.add(block["query"])
                rows, size = copy_block(streams[block["path"]], destination, block, digest)
                total_rows += rows
                total_bytes += size
                count += 1
            if count == 0:
                raise ValueError("Empty merged table")
            for key, stream in streams.items():
                for stat in (os.fstat(stream.fileno()), Path(sources[key]).stat()):
                    current = (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)
                    if current != identities[key]:
                        raise ValueError("Source changed during merge")
            destination.flush()
            os.fsync(destination.fileno())
        os.replace(partial, ready)
        directory = os.open(output, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
    finally:
        for stream in streams.values():
            stream.close()
    return dict(status="merged_candidate_requires_full_admission", query_blocks=count,
        rows=total_rows, bytes=total_bytes, sha256=digest.hexdigest(), path=str(ready),
        search_admitted=False, reuse_authorized=False, publication_ready=False)
