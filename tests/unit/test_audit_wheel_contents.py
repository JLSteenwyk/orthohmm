import base64
import hashlib
import zipfile

import pytest

from benchmark_tools.audit_wheel_contents import audit


def fixture_wheel(tmp_path, *, corrupt=False, extra=False):
    prefix = 'example-1.dist-info/'
    files = {prefix + 'METADATA': b'Name: example\nRequires-Dist: numpy\n',
             prefix + 'licenses/LICENSE.md': b'license\n'}
    rows = []
    for name, data in files.items():
        digest = base64.urlsafe_b64encode(hashlib.sha256(data).digest()).decode().rstrip('=')
        rows.append(f'{name},sha256={digest},{len(data)}\n')
    record = prefix + 'RECORD'
    files[record] = (''.join(rows) + f'{record},,\n').encode()
    if corrupt:
        files[prefix + 'METADATA'] += b'changed'
    if extra:
        files['unrecorded.txt'] = b'extra'
    wheel = tmp_path / 'test.whl'
    with zipfile.ZipFile(wheel, 'w') as archive:
        for name, data in files.items():
            archive.writestr(name, data)
    license_path = tmp_path / 'LICENSE.md'
    license_path.write_bytes(b'license\n')
    return wheel, license_path


def test_inventory(tmp_path):
    report = audit(*fixture_wheel(tmp_path))
    assert report['record_verified']
    assert report['embedded_project_license_matches']
    assert report['requires_dist'] == ['numpy']
    assert len(report['entries']) == 3
    assert not report['redistribution_clearance']


@pytest.mark.parametrize('option', ['corrupt', 'extra'])
def test_invalid_inventory(tmp_path, option):
    with pytest.raises(ValueError):
        audit(*fixture_wheel(tmp_path, **{option: True}))


def test_different_license(tmp_path):
    wheel, license_path = fixture_wheel(tmp_path)
    license_path.write_bytes(b'different')
    assert not audit(wheel, license_path)['embedded_project_license_matches']
