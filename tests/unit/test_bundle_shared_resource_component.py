from copy import deepcopy
import json
from pathlib import Path
import tarfile

import pytest

from benchmark_tools.results import bundle_shared_resource_component_20261003 as module

ROOT = Path(__file__).resolve().parents[2]


def payloads():
    data = {'sources/' + name: (ROOT / path).read_bytes() for name, path in module.SOURCES.items()}
    data['component.py'] = Path(module.__file__).read_bytes()
    data['LICENSE.md'] = (ROOT / 'LICENSE.md').read_bytes()
    data['README.md'] = b'Reporting only, not native reproduction.\n'
    data['requirements.txt'] = b'matplotlib==3.10.8\n'
    results = ROOT / 'benchmark_tools/results'
    for name in ('panel.json', 'attempts.tsv', 'cells.tsv'):
        data['data/' + name] = (results / 'threadripper_shared_panel_snapshot_20261003_v3' / name).read_bytes()
    for name in ('manifest.json', 'shared_threadripper_resources.png', 'shared_threadripper_resources.pdf', 'shared_threadripper_resources.svg'):
        data['figures/' + name] = (results / 'threadripper_shared_resource_figure_20261003_v2' / name).read_bytes()
    return data


def component(tmp_path):
    directory = tmp_path / 'component'
    data = payloads()
    for name, raw in data.items():
        path = directory / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(raw)
    manifest = dict(schema_version=1, source_commit='0' * 40, scope=module.SCOPE,
        native_inference_reproduced=False, raw_measurements_revalidated=False, publication_ready=False,
        reviewed_attempts=3, eligible_attempts=2, excluded_attempts=[0],
        files=[dict(path=name, **module.identity(raw)) for name, raw in sorted(data.items())])
    pin = rewrite_manifest(directory, manifest)
    return directory, manifest, pin


def rewrite_manifest(directory, manifest):
    path = directory / 'bundle.json'
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + '\n')
    return module.record(path)['sha256']


def test_pure_selection_does_not_execute_source_imports_or_collectors():
    source = 'import forbidden_missing_module\nX = 3\ndef chosen():\n    return X\ndef collector():\n    raise RuntimeError()\n'
    namespace = module.pure(source, {'X', 'chosen'}, {})
    assert namespace['chosen']() == 3
    assert 'collector' not in namespace and 'forbidden_missing_module' not in namespace


@pytest.mark.parametrize('source', ['X=1\nX=2\n', 'Y=1\n'])
def test_changed_pure_entry_points_refused(source):
    with pytest.raises(ValueError): module.pure(source, {'X'}, {})


def test_all_component_files_and_partial_coverage_verify(tmp_path):
    directory, _, pin = component(tmp_path)
    manifest, data = module.verify(directory, pin)
    assert len(data) == 15 and len(manifest['files']) == 15
    report, _, _ = module.reporting(data)
    assert report['reviewed_attempts'] == 3 and report['eligible_attempts'] == 2
    assert report['excluded_attempts'] == [0]
    assert all(cell['resources']['wall_seconds'] is None for cell in report['cells'])


@pytest.mark.parametrize('problem', ['hash', 'missing', 'extra', 'symlink', 'duplicate', 'traversal', 'scope', 'ready', 'coverage', 'schema'])
def test_changed_manifest_payload_or_scope_refused(tmp_path, problem):
    directory, original, pin = component(tmp_path)
    manifest = deepcopy(original)
    if problem == 'hash': (directory / 'data/attempts.tsv').write_bytes(b'changed')
    elif problem == 'missing': (directory / 'LICENSE.md').unlink()
    elif problem == 'extra': (directory / 'extra').write_bytes(b'changed')
    elif problem == 'symlink':
        path = directory / 'LICENSE.md'
        raw = path.read_bytes()
        external = tmp_path / 'external'
        external.write_bytes(raw)
        path.unlink()
        path.symlink_to(external)
    elif problem == 'duplicate': manifest['files'].append(manifest['files'][0])
    elif problem == 'traversal': manifest['files'][0]['path'] = '../outside'
    elif problem == 'scope': manifest['raw_measurements_revalidated'] = True
    elif problem == 'ready': manifest['publication_ready'] = True
    elif problem == 'coverage': manifest['eligible_attempts'] = 3
    else: manifest['schema_version'] = True
    if problem not in {'hash', 'missing', 'extra', 'symlink'}:
        pin = rewrite_manifest(directory, manifest)
    with pytest.raises((ValueError, FileNotFoundError)):
        module.verify(directory, pin)


def test_manifest_requires_external_anchor_and_regular_file(tmp_path):
    directory, _, pin = component(tmp_path)
    with pytest.raises(ValueError): module.verify(directory, '0' * 64)
    path = directory / 'bundle.json'
    original = tmp_path / 'manifest'
    original.write_bytes(path.read_bytes())
    path.unlink()
    path.symlink_to(original)
    with pytest.raises(ValueError): module.verify(directory, pin)


def test_actual_arithmetic_and_figure_replay_is_not_native_admission(tmp_path):
    directory, _, pin = component(tmp_path)
    output = tmp_path / 'replay'
    result = module.replay(directory, pin, output)
    assert result['equal_png_pixels'] is True and result['exact_table_files'] == 3
    assert result['excluded_attempts'] == [0]
    assert not result['raw_measurements_revalidated'] and not result['native_inference_reproduced']
    assert not result['publication_ready']
    assert json.loads((output / 'replay.json').read_text()) == result
    with pytest.raises(FileExistsError): module.replay(directory, pin, output)
    with pytest.raises(FileExistsError): module.replay(directory, pin, directory / 'mutate')


def test_retained_archive_has_exact_regular_payload_and_honest_scope():
    results = ROOT / 'benchmark_tools/results'
    receipt = json.loads((results / 'shared_resource_reporting_component_result_20261003.json').read_text())
    archive = results / 'shared_resource_reporting_component_20261003_v1.tar.gz'
    assert module.record(archive)['sha256'] == receipt['archive']['sha256']
    assert receipt['archive']['sha256'] == 'd14fc99f9ee4580e4b3a57395420f52148b825b2cb2cf0ba7716b57a11e8f6bc'
    with tarfile.open(archive, 'r:gz') as handle:
        members = handle.getmembers()
        assert {member.name for member in members} == module.EXPECTED | {'bundle.json'}
        assert len(members) == 16 and all(member.isfile() for member in members)
        data = {member.name: handle.extractfile(member).read() for member in members}
    assert module.identity(data['bundle.json'])['sha256'] == receipt['manifest']['sha256']
    manifest = json.loads(data['bundle.json'])
    for pin in manifest['files']:
        module.checked(data[pin['path']], pin)
    report = json.loads(data['data/panel.json'])
    assert report['reviewed_attempts'] == 3 and report['eligible_attempts'] == 2
    assert report['excluded_attempts'] == [0]
    observation = receipt['guarded_replay_observation']
    assert observation['original_checkout_canary_denied'] is True
    assert observation['post_canary_original_checkout_open_events'] == 0
    assert observation['project_modules_imported'] is False
    assert observation['os_containment_established'] is False
    assert not receipt['native_inference_reproduced'] and not receipt['raw_measurements_revalidated']
    assert not receipt['publication_ready']
