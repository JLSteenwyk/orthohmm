"""Inventory wheel bytes and native linkage without loading native code."""

import argparse
import base64
import csv
from email.parser import BytesParser
import hashlib
import io
import json
from pathlib import Path
import subprocess
import tempfile
import zipfile


def audit(wheel, license_path):
    with zipfile.ZipFile(wheel) as archive:
        names = archive.namelist()
        if len(names) != len(set(names)):
            raise ValueError("Duplicate ZIP entries")
        records = [name for name in names if name.endswith('.dist-info/RECORD')]
        if len(records) != 1:
            raise ValueError("Expected one RECORD")
        prefix = records[0].rsplit('/', 1)[0]
        metadata = BytesParser().parsebytes(archive.read(prefix + '/METADATA'))
        entries = []
        recorded = set()
        for name, digest, size in csv.reader(io.StringIO(archive.read(records[0]).decode())):
            if name in recorded:
                raise ValueError("Duplicate RECORD entry")
            recorded.add(name)
            data = archive.read(name)
            sha = hashlib.sha256(data).hexdigest()
            if name != records[0]:
                expected = 'sha256=' + base64.urlsafe_b64encode(bytes.fromhex(sha)).decode().rstrip('=')
                if digest != expected or int(size) != len(data):
                    raise ValueError('RECORD mismatch: ' + name)
            elif digest or size:
                raise ValueError('RECORD must not hash itself')
            item = dict(path=name, bytes=len(data), sha256=sha)
            if name.endswith('.so'):
                with tempfile.TemporaryDirectory() as directory:
                    native = Path(directory) / 'library.so'
                    native.write_bytes(data)
                    item['readelf_dynamic'] = subprocess.check_output(
                        ['readelf', '-d', str(native)], text=True)
                    symbols = subprocess.check_output(['nm', '-C', str(native)], text=True)
                    item['cuda_runtime_symbol_witnesses'] = [line for line in symbols.splitlines()
                        if any(token in line for token in ('cudaMalloc', 'cudaFree', '__cudaRegisterFatBinary'))]
            entries.append(item)
        if recorded != set(names):
            raise ValueError('RECORD does not cover ZIP inventory')
        embedded = prefix + '/licenses/LICENSE.md'
        license_matches = archive.read(embedded) == license_path.read_bytes()
    return dict(wheel=str(wheel), wheel_sha256=hashlib.sha256(wheel.read_bytes()).hexdigest(),
                wheel_bytes=wheel.stat().st_size, record_verified=True,
                embedded_project_license_matches=license_matches,
                requires_dist=metadata.get_all('Requires-Dist', []), entries=entries,
                redistribution_clearance=False,
                scope='Byte inventory and direct ELF linkage only; not a complete static dependency or legal audit')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('wheel', type=Path)
    parser.add_argument('--license', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(audit(args.wheel, args.license), indent=2))
