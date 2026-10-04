import copy
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tarfile

import pytest

from benchmark_tools.results import bundle_shared_prepared_component_20261004 as component

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / 'benchmark_tools/results'
SNAPSHOT = RESULTS / 'threadripper_shared_panel_snapshot_20261004_v22'
FIGURE = RESULTS / 'threadripper_shared_resource_figure_20261004_v21'


@pytest.fixture
def bundle(tmp_path, monkeypatch):
    def committed(argv, **kwargs):
        if argv[-1] == 'HEAD':
            return 'a' * 40 + '\n'
        return (ROOT / argv[-1].split(':', 1)[1]).read_bytes()

    with monkeypatch.context() as patch:
        patch.setattr(component.subprocess, 'check_output', committed)
        result = component.build(ROOT, SNAPSHOT, FIGURE, tmp_path / 'component')
    return Path(result['manifest']['path']).parent, result['manifest']['sha256']


def reanchor(directory, transform):
    path = directory / 'bundle.json'
    manifest = json.loads(path.read_bytes())
    transform(manifest)
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + '\n')
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_actual_prepared_arithmetic_and_archive(bundle, tmp_path):
    directory, digest = bundle
    manifest, payloads = component.verify(directory, digest)
    assert manifest['schema_version'] == 2
    assert manifest['reviewed_attempts'] == 22
    assert manifest['resource_reviewed_attempts'] == 20
    assert manifest['eligible_attempts'] == 19
    assert manifest['excluded_attempts'] == [0, 17, 20]
    assert manifest['pre_native_aborted_indices'] == [17, 20]
    report, _, _ = component.reporting(payloads)
    assert report['runs'][17]['resources'] is None
    assert report['runs'][20]['resources'] is None
    assert report['runs'][0]['resources']['wall_seconds'] > 0
    assert sum(c['eligible_repeats'] == 3 for c in report['cells']) == 2
    archive = directory.with_name(directory.name + '.tar.gz')
    restored = tmp_path / 'restored'
    restored.mkdir()
    with tarfile.open(archive) as handle:
        assert len(handle.getmembers()) == len(component.EXPECTED) + 1
        assert all(m.isfile() for m in handle.getmembers())
        handle.extractall(restored, filter='data')
    command = [sys.executable, '-I', '-S', '-B', str(restored / 'component.py'),
               'verify', '--directory', str(restored), '--manifest-sha256', digest]
    run = subprocess.run(command, capture_output=True, text=True, check=True, cwd=tmp_path)
    assert json.loads(run.stdout)['publication_ready'] is False
    result = component.replay(restored, digest, tmp_path / 'replay')
    assert result['exact_table_files'] == 3 and result['equal_png_pixels'] is True
    assert result['raw_measurements_revalidated'] is False
    assert result['native_inference_reproduced'] is False


@pytest.mark.parametrize('mutation', ['digest', 'bytes', 'extra', 'symlink', 'missing'])
def test_invalid_inventory_refused(bundle, mutation):
    directory, digest = bundle
    path = directory / 'data/attempts.tsv'
    if mutation == 'digest':
        digest = '0' * 64
    elif mutation == 'bytes':
        path.write_bytes(path.read_bytes() + b'changed')
    elif mutation == 'extra':
        (directory / 'unexpected.txt').write_text('extra')
    elif mutation == 'symlink':
        original = path.read_bytes()
        path.unlink()
        target = directory.parent / 'outside.tsv'
        target.write_bytes(original)
        path.symlink_to(target)
    else:
        path.unlink()
    with pytest.raises(ValueError):
        component.verify(directory, digest)


@pytest.mark.parametrize('field,value', [
    ('schema_version', True), ('schema_version', 1),
    ('reviewed_attempts', 27), ('resource_reviewed_attempts', 22),
    ('excluded_attempts', [0]), ('pre_native_aborted_indices', []),
    ('all_planned_attempts_reviewed', True), ('publication_ready', True),
    ('raw_measurements_revalidated', True), ('native_inference_reproduced', True),
])
def test_reanchored_scope_or_coverage_cannot_promote_results(bundle, field, value):
    directory, _ = bundle
    digest = reanchor(directory, lambda m: m.update({field: value}))
    with pytest.raises(ValueError):
        component.verify(directory, digest)


def test_recomputed_panel_rejects_imputed_abort_or_edited_median(bundle):
    directory, digest = bundle
    _, original = component.verify(directory, digest)
    for modify in (
        lambda r: r['runs'][17].update(resources={'wall_seconds': 0}),
        lambda r: r.update(eligible_attempts=20),
        lambda r: r['cells'][0]['resources'].update(wall_seconds={'median': 0}),
    ):
        payloads = copy.copy(original)
        report = json.loads(payloads['data/panel.json'])
        modify(report)
        payloads['data/panel.json'] = json.dumps(report).encode()
        with pytest.raises(ValueError):
            component.reporting(payloads)


def test_output_guards_and_old_helper_unchanged(bundle, tmp_path):
    directory, digest = bundle
    with pytest.raises(FileExistsError):
        component.replay(directory, digest, directory / 'generated')
    occupied = tmp_path / 'occupied'
    occupied.mkdir()
    with pytest.raises(FileExistsError):
        component.replay(directory, digest, occupied)
    with pytest.raises(FileExistsError):
        component.build(ROOT, SNAPSHOT, FIGURE, directory)
    assert component.identity((RESULTS / component.CORE).read_bytes())['sha256'] == component.CORE_SHA
