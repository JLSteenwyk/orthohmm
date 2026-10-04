import json
import shutil
import subprocess
import sys

import pytest

from tests.unit.test_bundle_publication_source import exported, git, rewrite_index
from tests.unit import test_bundle_publication_source as legacy


@pytest.fixture
def source_module():
    return legacy.module


def commit_helpers(repo):
    root = repo / 'benchmark_tools/results'
    root.mkdir(parents=True, exist_ok=True)
    (root / 'fixture_dep.py').write_text('VALUE = 42\n')
    (root / 'fixture_source.py').write_text('from .fixture_dep import VALUE\n')
    (repo / 'tests/unit/test_example.py').write_text(
        'from benchmark_tools.results.fixture_source import VALUE\n'
        'def test_example():\n    assert VALUE == 42\n')
    git(repo, 'add', 'benchmark_tools/results/fixture_dep.py',
        'benchmark_tools/results/fixture_source.py', 'tests/unit/test_example.py')
    git(repo, 'commit', '-qm', 'Result helper fixture')
    # Export must read Git, not mutable working source.
    (root / 'fixture_source.py').write_text('not valid python !')
    return git(repo, 'rev-parse', 'HEAD')


def test_result_helpers_relocated_import_collection_and_verifier(exported, source_module, tmp_path):
    repo, _, original = exported
    revision = commit_helpers(repo)
    output = tmp_path / 'with-helpers'
    original_index = json.loads((tmp_path / 'bundle/SOURCE_INDEX.json').read_text())
    profile = original_index.get('profile', 'source-only')
    result = source_module.build(repo, revision, output, profile, include_result_helpers=True)
    index = json.loads((output / 'SOURCE_INDEX.json').read_text())
    assert index['schema'] == 'publication_source_components_v3'
    assert index['include_result_helpers'] is True
    assert result['result_helpers_included'] is True
    assert result['executable_benchmark_reproduced'] is False
    assert result['publication_ready'] is False
    assert not (output / 'workflow/benchmark_tools/results/private.json').exists()
    assert (output / 'workflow/benchmark_tools/results/fixture_source.py').read_text() == 'from .fixture_dep import VALUE\n'
    assert (output / 'scientific/orthohmm/version.py').read_text() == "VERSION = 'fixture'\n"
    legacy = tmp_path / 'legacy-again'
    source_module.build(repo, revision, legacy, profile)
    assert not (legacy / 'workflow/benchmark_tools/results/fixture_source.py').exists()
    assert json.loads((legacy / 'SOURCE_INDEX.json').read_text())['schema'] == original_index['schema']

    relocated = tmp_path / 'relocated'
    shutil.copytree(output, relocated)
    shutil.rmtree(output)
    shutil.rmtree(repo)
    command = [sys.executable, '-I', '-S', '-B', str(relocated / 'workflow' / source_module.RUNNER),
               'verify', str(relocated), '--manifest-sha256', result['manifest']['sha256']]
    checked = subprocess.run(command, cwd=tmp_path, env={'PATH': '/no-git-here'},
                             capture_output=True, text=True, timeout=10)
    assert checked.returncode == 0, checked.stderr
    assert json.loads(checked.stdout) == result
    workflow = relocated / 'workflow'
    imported = subprocess.run([sys.executable, '-I', '-S', '-B', '-c',
        f'import sys; sys.path.insert(0, {str(workflow)!r}); '
        'from benchmark_tools.results.fixture_source import VALUE; assert VALUE == 42'],
        cwd=tmp_path, env={'PATH': '/no-git-here'}, capture_output=True, text=True, timeout=10)
    assert imported.returncode == 0, imported.stderr
    collected = subprocess.run([sys.executable, '-I', '-B', '-c',
        f'import sys; sys.path.insert(0, {str(workflow)!r}); import pytest; '
        f'raise SystemExit(pytest.main(["--collect-only", "-q", {str(workflow / "tests/unit/test_example.py")!r}]))'],
        cwd=tmp_path, capture_output=True, text=True, timeout=15)
    assert collected.returncode == 0, collected.stderr
    assert '1 test collected' in collected.stdout


@pytest.mark.parametrize('flag', [False, None, 1, 'true'])
def test_current_schema_requires_exact_boolean(exported, source_module, tmp_path, flag):
    repo, _, _ = exported
    revision = commit_helpers(repo)
    output = tmp_path / 'with-helpers'
    source_module.build(repo, revision, output, include_result_helpers=True)
    digest = rewrite_index(output, lambda index: index.update(include_result_helpers=flag))
    with pytest.raises(ValueError, match='explicitly include result helpers'):
        source_module.verify(output, digest)


@pytest.mark.parametrize('exported', ['orthobench-inputs'], indirect=True)
@pytest.mark.parametrize('schema', ['publication_source_components_v1', 'publication_source_components_v2'])
def test_historical_schema_cannot_enable_helpers(exported, source_module, tmp_path, schema):
    repo, _, _ = exported
    revision = commit_helpers(repo)
    output = tmp_path / 'with-helpers'
    source_module.build(repo, revision, output, 'orthobench-inputs', include_result_helpers=True)
    digest = rewrite_index(output, lambda index: index.update(schema=schema))
    with pytest.raises(ValueError, match='Historical source schema'):
        source_module.verify(output, digest)


def test_all_observed_missing_helpers_selected_without_data(source_module):
    helpers = (
        'bundle_shared_prepared_component_20261004', 'bundle_shared_resource_component_20261003',
        'capture_shared_terminal_controller_20261004', 'continue_shared_goal_reaffirmed_20261004',
        'continue_shared_prenative_panel_20261004', 'continue_shared_prepared_panel_20261004',
        'plot_shared_prenative_resources_20261004', 'plot_shared_prepared_resources_20261004',
        'plot_shared_threadripper_resources_20261003', 'prepare_shared_threadripper_launch_20261003',
        'render_shared_resource_section_20261004', 'report_shared_prenative_panel_20261004',
        'report_shared_prepared_panel_20261004', 'report_shared_threadripper_panel_20261003',
        'review_shared_deadline_failure_20261004', 'review_shared_prenative_failure_20261004',
        'review_shared_prenative_panel_20261004', 'review_shared_prepared_panel_20261004',
        'snapshot_shared_prenative_panel_20261004', 'snapshot_shared_prepared_panel_20261004',
    )
    for helper in helpers:
        path = f'benchmark_tools/results/{helper}.py'
        assert source_module.selected(path, 'workflow', 'native-build', True)
        assert not source_module.selected(path, 'workflow', 'native-build')
    for name in ['benchmark_tools/results/private.json', 'benchmark_tools/results/raw.fa',
                 'benchmark_tools/results/nested/helper.py', 'tests/samples/input.fa']:
        assert not source_module.selected(name, 'workflow', 'native-build', True)


def test_builder_refuses_non_boolean_before_output(source_module, tmp_path):
    with pytest.raises(ValueError, match='explicit boolean'):
        source_module.build(tmp_path, 'HEAD', tmp_path / 'never-created', include_result_helpers=1)
    assert not (tmp_path / 'never-created').exists()
