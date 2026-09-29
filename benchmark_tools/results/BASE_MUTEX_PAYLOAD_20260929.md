# Base Mutex Payload Scope

The [exact-archive inspection](base_mutex_payload_20260929.json) narrows the
missing-notice finding for `_libgcc_mutex` 0.1, build `main`. Its retained
3,473-byte archive matches the original notice inventory SHA256. Inspection
of both decompressed tar streams found ten regular metadata members and
**zero payload members**. `info/paths.json` lists no installed paths and
`info/files` is empty. The embedded index has no dependencies. Neither the
index nor `info/about.json` declares a license; the embedded recipes also
contain no license field.

Thus this particular archive does not contain an unidentified executable or
library payload. It is not the archive containing the actual GCC runtime.
Do not transfer this conclusion to `libgcc-ng`, other packages or other
mutex builds. The package metadata/recipe archive still has no license
declaration; no redistribution permission is inferred, and no package has
been removed or replaced in the frozen environment.

This is a content finding from retained bytes, not a legal clearance. The
broader component attribution, source obligations and release review remain
open. No runtime, benchmark output or scientific conclusion changed.

## Replay

Run from the repository root with the existing Python environment containing
`zstandard` (inspection used version 0.19.0). This reads archives without
extracting files or executing package code:

```python
import hashlib
import io
import json
from pathlib import Path
import tarfile
import zipfile
import zstandard

r = json.loads(Path("benchmark_tools/results/base_mutex_payload_20260929.json").read_text())
raw = Path(r["archive"]["path"]).read_bytes()
assert len(raw) == r["archive"]["bytes"]
assert hashlib.sha256(raw).hexdigest() == r["archive"]["sha256"]
with zipfile.ZipFile(io.BytesIO(raw)) as archive:
    assert archive.namelist() == r["zip_members"]
    assert len(archive.namelist()) == len(set(archive.namelist()))
    for stream, expected in r["streams"].items():
        observed = []
        with zstandard.ZstdDecompressor().stream_reader(io.BytesIO(archive.read(stream))) as reader:
            with tarfile.open(fileobj=reader, mode="r|") as members:
                for member in members:
                    assert member.isfile()
                    data = members.extractfile(member).read()
                    assert len(data) == member.size
                    observed.append(dict(name=member.name, bytes=member.size,
                                         sha256=hashlib.sha256(data).hexdigest()))
        assert observed == expected
assert r["streams"]["pkg-_libgcc_mutex-0.1-main.tar.zst"] == []
print("Exact metadata inventory and empty payload verified")
```
