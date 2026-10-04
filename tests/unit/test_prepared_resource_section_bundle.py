import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile

import pytest

from benchmark_tools.results import bundle_shared_prepared_component_20261004 as component

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / 'benchmark_tools/results'
SNAPSHOT = RESULTS / 'threadripper_shared_panel_snapshot_20261004_v24'
FIGURE = RESULTS / 'threadripper_shared_resource_figure_20261004_v23'
SECTION = RESULTS / 'resource_manuscript_section_20261004_v24'


def build(tmp_path, monkeypatch, section=SECTION):
    def committed(argv, **kwargs):
        if argv[-1] == 'HEAD':
            return 'a' * 40 + '\n'
        return (ROOT / argv[-1].split(':', 1)[1]).read_bytes()

    with monkeypatch.context() as patch:
        patch.setattr(component.subprocess, 'check_output', committed)
        return component.build(ROOT, SNAPSHOT, FIGURE, tmp_path / 'bundle', section)


def test_section_archive_replays_outside_checkout(tmp_path, monkeypatch):
    result = build(tmp_path, monkeypatch)
    restored = tmp_path / 'restored'
    restored.mkdir()
    with tarfile.open(result['archive']['path']) as handle:
        assert len(handle.getmembers()) == len(component.EXPECTED | component.SECTION_MEMBERS) + 1
        assert all(member.isfile() for member in handle.getmembers())
        handle.extractall(restored, filter='data')
    manifest, payloads = component.verify(restored, result['manifest']['sha256'])
    assert manifest['schema_version'] == 3 and manifest['reviewed_attempts'] == 24
    assert payloads['manuscript/resource_section.md'] == (SECTION / 'resource_section.md').read_bytes()
    common = [str(restored / 'component.py'), 'verify', '--directory', str(restored),
              '--manifest-sha256', result['manifest']['sha256']]
    verified = subprocess.run([sys.executable, '-I', '-S', '-B', *common],
                              cwd=tmp_path, capture_output=True, text=True, check=True)
    assert json.loads(verified.stdout)['publication_ready'] is False
    common[1] = 'replay'
    replay = tmp_path / 'replayed'
    run = subprocess.run([sys.executable, '-I', '-B', *common, '--output', str(replay)],
                         cwd=tmp_path, capture_output=True, text=True, check=True)
    receipt = json.loads(run.stdout)
    assert receipt['exact_resource_section'] is True
    assert receipt['exact_table_files'] == 3 and receipt['equal_png_pixels'] is True
    assert receipt['native_inference_reproduced'] is False
    assert receipt['raw_measurements_revalidated'] is False
    assert (replay / 'resource_section.md').read_bytes() == payloads['manuscript/resource_section.md']


@pytest.mark.parametrize('mutation', ['table', 'source', 'coverage', 'promotion', 'prose'])
def test_mismatched_section_refused_before_output(tmp_path, monkeypatch, mutation):
    copied = tmp_path / 'section'
    shutil.copytree(SECTION, copied)
    path = copied / 'manifest.json'
    receipt = json.loads(path.read_bytes())
    if mutation == 'table':
        receipt['table']['sha256'] = '0' * 64
    elif mutation == 'source':
        receipt['sources'][0]['sha256'] = '0' * 64
    elif mutation == 'coverage':
        receipt['reviewed_attempts'] = 27
    elif mutation == 'promotion':
        receipt['publication_ready'] = True
    else:
        prose = copied / 'resource_section.md'
        prose.write_bytes(prose.read_bytes() + b'Invented superiority.\n')
        receipt['section'].update(component.identity(prose.read_bytes()))
    path.write_text(json.dumps(receipt))
    with pytest.raises(ValueError):
        build(tmp_path, monkeypatch, copied)
    assert not (tmp_path / 'bundle').exists()


def test_reanchored_prose_tampering_is_not_accepted(tmp_path, monkeypatch):
    result = build(tmp_path, monkeypatch)
    directory = Path(result['manifest']['path']).parent
    prose = directory / 'manuscript/resource_section.md'
    prose.write_bytes(prose.read_bytes() + b'Changed conclusion.\n')
    section = directory / 'manuscript/manifest.json'
    receipt = json.loads(section.read_bytes())
    receipt['section'].update(component.identity(prose.read_bytes()))
    section.write_text(json.dumps(receipt))
    manifest = directory / 'bundle.json'
    index = json.loads(manifest.read_bytes())
    for row in index['files']:
        if row['path'] in {'manuscript/resource_section.md', 'manuscript/manifest.json'}:
            row.update(component.identity((directory / row['path']).read_bytes()))
    manifest.write_text(json.dumps(index))
    digest = hashlib.sha256(manifest.read_bytes()).hexdigest()
    with pytest.raises(ValueError, match='prose replay'):
        component.verify(directory, digest)


def test_section_optional_without_changing_old_inventory(tmp_path, monkeypatch):
    result = build(tmp_path, monkeypatch, None)
    directory = Path(result['manifest']['path']).parent
    manifest, payloads = component.verify(directory, result['manifest']['sha256'])
    assert manifest['schema_version'] == 2
    assert set(payloads) == component.EXPECTED
