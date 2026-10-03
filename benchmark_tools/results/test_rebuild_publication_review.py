"""Presentation-only replay checks; no scientific workflow is invoked."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import fitz
import pytest


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / 'benchmark_tools/results'
ASSEMBLY = RESULTS / 'publication_main_with_figures_20261002_v3/assembly.json'


def relative(ref):
    return Path('benchmark_tools/results') / ref['path'].split('/benchmark_tools/results/')[1]


@pytest.fixture
def checkout(tmp_path):
    accepted = json.loads(ASSEMBLY.read_text())
    paths = [Path('benchmark_tools/results/rebuild_publication_review.py'),
             ASSEMBLY.relative_to(ROOT),
             relative(accepted['preserved_source_pages'][0]['source']),
             *[relative(item['pdf']) for item in accepted['figures']]]
    # Export only staged Git objects, never retained untracked/raw dependencies.
    for path in paths:
        contents = subprocess.check_output(['git', 'show', ':' + str(path)], cwd=ROOT)
        destination = tmp_path / path
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(contents)
    assert len(list(tmp_path.rglob('*.*'))) == len(paths) == 17
    return tmp_path


def run(checkout, optimized=False):
    command = [sys.executable, '-B']
    if optimized:
        command.append('-O')
    return subprocess.run(command + [str(checkout / 'benchmark_tools/results/rebuild_publication_review.py'),
                          '--output', str(checkout / 'rebuilt')],
                          cwd=checkout, text=True, capture_output=True, check=False)


def signature(link):
    value = {key: item for key, item in link.items() if key not in ('xref', 'id')}
    for key in ('from', 'to'):
        if isinstance(value.get(key), (fitz.Rect, fitz.Point)):
            value[key] = list(value[key])
    return json.dumps(value, sort_keys=True)


@pytest.mark.parametrize('optimized', [False, True])
def test_relocated_presentation_matches_accepted(checkout, optimized):
    completed = run(checkout, optimized)
    assert completed.returncode == 0, completed.stderr
    assert completed.stderr == ''
    report = json.loads((checkout / 'rebuilt/assembly.json').read_text())
    assert report['publication_ready'] is False
    assert len(report['preserved_source_pages']) == 23
    assert len(report['redirected_main_figure_links']) == 6
    assert len(report['restored_original_file_uri_actions']) == 45
    for ref in report['inputs']:
        assert Path(ref['path']).is_relative_to(checkout)
        assert hashlib.sha256(Path(ref['path']).read_bytes()).hexdigest() == ref['sha256']
    accepted = json.loads(ASSEMBLY.read_text())
    with fitz.open(checkout / 'rebuilt/document.pdf') as rebuilt, fitz.open(RESULTS / 'publication_main_with_figures_20261002_v3/document.pdf') as original:
        assert len(rebuilt) == len(original) == 27
        assert not rebuilt.is_repaired
        assert rebuilt.get_toc() == original.get_toc() == accepted['bookmarks']
        for before, after in zip(original, rebuilt):
            assert before.rect == after.rect
            assert before.get_text('words') == after.get_text('words')
            a, b = before.get_pixmap(alpha=False), after.get_pixmap(alpha=False)
            assert (a.width, a.height, a.n, a.samples) == (b.width, b.height, b.n, b.samples)
            assert sorted(map(signature, before.get_links())) == sorted(map(signature, after.get_links()))


@pytest.mark.parametrize('optimized', [False, True])
def test_existing_output_refused_without_writes(checkout, optimized):
    directory = checkout / 'rebuilt'
    directory.mkdir()
    marker = directory / 'untouched'
    marker.write_bytes(b'preserve')
    completed = run(checkout, optimized)
    assert completed.returncode != 0 and 'Refusing existing output directory' in completed.stderr
    assert list(directory.iterdir()) == [marker]
    assert marker.read_bytes() == b'preserve'


@pytest.mark.parametrize('optimized', [False, True])
@pytest.mark.parametrize('receipt', [False, True])
def test_corrupted_input_refused_before_output(checkout, optimized, receipt):
    accepted = json.loads(ASSEMBLY.read_text())
    target = ASSEMBLY.relative_to(ROOT) if receipt else relative(accepted['figures'][0]['pdf'])
    path = checkout / target
    path.write_bytes(path.read_bytes() + b'corrupt')
    completed = run(checkout, optimized)
    assert completed.returncode != 0
    assert 'checksum mismatch' in completed.stderr
    assert not (checkout / 'rebuilt').exists()
