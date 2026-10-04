"""Portable prepared-panel reporting, preserving measured and pre-native failures."""
import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import re
from statistics import median
import subprocess
import tarfile
from types import ModuleType, SimpleNamespace

CORE = 'bundle_shared_resource_component_20261003.py'
CORE_SHA = '98d723bb70910eed22a8ab75ef4e17f0332490ada309121f607fc0ed3d5b5781'
core_bytes = Path(__file__).with_name(CORE).read_bytes()
if hashlib.sha256(core_bytes).hexdigest() != CORE_SHA:
    raise ValueError('Historical reporting helper differs')
core = ModuleType('retained_reporting_core')
core.__file__ = str(Path(__file__).with_name(CORE))
exec(compile(core_bytes, core.__file__, 'exec'), core.__dict__)
identity, record, checked, safe_name, pure = (
    core.identity, core.record, core.checked, core.safe_name, core.pure)
SCOPE = core.SCOPE
SOURCES = {
    'plan.py': 'benchmark_tools/prepare_scaling_inputs.py',
    'scopes.py': 'benchmark_tools/derive_threadripper_resources.py',
    'identity.py': 'benchmark_tools/run_threadripper_scaling.py',
    'prior_tables.py': 'benchmark_tools/results/report_shared_threadripper_panel_20261003.py',
    'tables.py': 'benchmark_tools/results/report_shared_prepared_panel_20261004.py',
    'prior_figure.py': 'benchmark_tools/results/plot_shared_threadripper_resources_20261003.py',
    'figure.py': 'benchmark_tools/results/plot_shared_prepared_resources_20261004.py',
}
EXPECTED = {
    'component.py', CORE, 'LICENSE.md', 'requirements.txt', 'README.md',
    'data/panel.json', 'data/attempts.tsv', 'data/cells.tsv',
    'figures/shared_threadripper_resources.png',
    'figures/shared_threadripper_resources.pdf',
    'figures/shared_threadripper_resources.svg', 'figures/manifest.json',
    *['sources/' + name for name in SOURCES],
}
SECTION_SOURCE = 'benchmark_tools/results/render_shared_resource_section_20261004.py'
SECTION_MEMBERS = {'sources/section.py', 'manuscript/resource_section.md',
                   'manuscript/manifest.json'}


def section_text(payloads):
    report, _, figure = reporting(payloads)
    constants = pure(payloads['sources/prior_figure.py'], {'METHODS'}, {})
    plotter = SimpleNamespace(validate=figure['validate'], **constants)
    renderer = pure(payloads['sources/section.py'], {'METRICS', 'LABELS', 'render'},
                    dict(plotter=plotter))
    return renderer['render'](report).encode()


def check_section(payloads):
    report = json.loads(payloads['data/panel.json'])
    receipt = json.loads(payloads['manuscript/manifest.json'])
    if receipt['status'] != 'reviewed_resource_section_generated':
        raise ValueError('Resource section generation differs')
    checked(payloads['data/panel.json'], receipt['table'])
    checked(payloads['manuscript/resource_section.md'], receipt['section'])
    source_names = {'section.py', 'figure.py', 'tables.py', 'prior_figure.py', 'prior_tables.py'}
    mapping = {Path(SECTION_SOURCE if name == 'section.py' else SOURCES[name]).name: name
               for name in source_names}
    seen = set()
    for pin in receipt['sources']:
        name = mapping.get(Path(pin['path']).name)
        if name is None or name in seen:
            raise ValueError('Resource section source inventory differs')
        checked(payloads['sources/' + name], pin)
        seen.add(name)
    if seen != source_names:
        raise ValueError('Resource section source inventory is incomplete')
    for key, report_key in (
        ('reviewed_attempts', 'reviewed_attempts'),
        ('measured_attempts', 'resource_reviewed_attempts'),
        ('eligible_attempts', 'eligible_attempts'),
        ('excluded_attempts', 'excluded_attempts'),
        ('all_planned_attempts_reviewed', 'all_planned_attempts_reviewed'),
    ):
        if receipt[key] != report[report_key] or type(receipt[key]) is not type(report[report_key]):
            raise ValueError('Resource section coverage differs: ' + key)
    if any(receipt[key] is not False for key in (
            'native_inference_repeated', 'raw_audit_repeated', 'publication_ready')):
        raise ValueError('Resource section scope differs')
    if payloads['manuscript/resource_section.md'] != section_text(payloads):
        raise ValueError('Resource section prose replay differs')


def reporting(payloads, plotting=False):
    plan = pure(payloads['sources/plan.py'], {'METHODS', 'SIZES', 'planned_runs'}, {})
    scopes = pure(payloads['sources/scopes.py'], {'SCOPES'}, {})['SCOPES']
    expect = pure(payloads['sources/identity.py'], {'expect'}, {})['expect']
    executor = SimpleNamespace(record=record, expect=expect,
                               SHARED_SCOPE='shared_host_matched_resources')
    namespace = dict(math=math, median=median, json=json, csv=csv,
                     executor=executor, SCOPES=scopes, **plan)
    prior = pure(payloads['sources/prior_tables.py'],
                 {'METRICS', 'IDENTITY', 'summarize', 'export'}, namespace)
    tables = pure(payloads['sources/tables.py'], {'summarize'},
                  dict(namespace, previous=SimpleNamespace(**prior)))
    figure_constants = pure(payloads['sources/prior_figure.py'],
                            {'METHODS', 'METRICS'}, {})
    figure_namespace = dict(math=math, executor=executor,
        planned_runs=plan['planned_runs'],
        tables=SimpleNamespace(summarize=tables['summarize'], SCOPES=scopes),
        **figure_constants)
    if plotting:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
        figure_namespace.update(plt=plt, Line2D=Line2D)
    figure = pure(payloads['sources/figure.py'], {'validate', 'plot'}, figure_namespace)
    report = json.loads(payloads['data/panel.json'])
    figure['validate'](report)
    return report, prior, figure


def verify(directory, manifest_sha):
    directory = Path(directory).resolve(strict=True)
    path = directory / 'bundle.json'
    if path.is_symlink() or identity(path.read_bytes())['sha256'] != manifest_sha:
        raise ValueError('Externally pinned manifest differs')
    manifest = json.loads(path.read_bytes())
    if (type(manifest['schema_version']) is not int or manifest['schema_version'] not in (2, 3)
            or manifest['scope'] != SCOPE
            or not re.fullmatch('[0-9a-f]{40}', manifest['source_commit'])
            or any(manifest[key] is not False for key in (
                'native_inference_reproduced', 'raw_measurements_revalidated', 'publication_ready'))):
        raise ValueError('Prepared reporting-only scope differs')
    expected = EXPECTED | (SECTION_MEMBERS if manifest['schema_version'] == 3 else set())
    payloads = {}
    for pin in manifest['files']:
        name = safe_name(pin['path'])
        member = directory / name
        if (name not in expected or name in payloads or member.is_symlink()
                or not member.is_file() or not member.resolve().is_relative_to(directory)):
            raise ValueError('Invalid reporting payload')
        payloads[name] = member.read_bytes()
        checked(payloads[name], pin)
    actual = {p.relative_to(directory).as_posix() for p in directory.rglob('*')
              if p.is_file() or p.is_symlink()}
    if set(payloads) != expected or actual != expected | {'bundle.json'}:
        raise ValueError('Missing or extra reporting member')
    if payloads['component.py'] != Path(__file__).read_bytes() or payloads[CORE] != core_bytes:
        raise ValueError('Use exact bundled reporting readers')
    report, _, _ = reporting(payloads)
    checked(payloads['sources/tables.py'], report['source'])
    checked(payloads['sources/prior_tables.py'], report['prior_reporter'])
    for key in ('reviewed_attempts', 'resource_reviewed_attempts', 'eligible_attempts',
                'excluded_attempts', 'pre_native_aborted_indices', 'all_planned_attempts_reviewed'):
        if manifest[key] != report[key] or type(manifest[key]) is not type(report[key]):
            raise ValueError('Manifest coverage differs: ' + key)
    figure = json.loads(payloads['figures/manifest.json'])
    checked(payloads['data/panel.json'], figure['source_results'])
    checked(payloads['sources/figure.py'], figure['plotter'])
    names = set()
    for pin in figure['outputs']:
        name = 'figures/' + Path(pin['path']).name
        if name not in EXPECTED or name in names:
            raise ValueError('Invalid figure output inventory')
        checked(payloads[name], pin)
        names.add(name)
    if names != {f'figures/shared_threadripper_resources.{ext}' for ext in ('png', 'pdf', 'svg')}:
        raise ValueError('Incomplete figure inventory')
    if payloads['requirements.txt'].decode() != 'matplotlib==' + figure['matplotlib'] + '\n':
        raise ValueError('Reporting dependency differs')
    if manifest['schema_version'] == 3:
        check_section(payloads)
    return manifest, payloads


def replay(directory, manifest_sha, output):
    directory = Path(directory).resolve(strict=True)
    output = Path(output).absolute()
    output = output.parent.resolve() / output.name
    if output.exists() or output.is_symlink() or output.is_relative_to(directory):
        raise FileExistsError('Require fresh output outside immutable component')
    manifest, payloads = verify(directory, manifest_sha)
    import matplotlib
    if payloads['requirements.txt'].decode() != 'matplotlib==' + matplotlib.__version__ + '\n':
        raise ValueError('Use recorded reporting dependency')
    report, tables, figure = reporting(payloads, plotting=True)
    output.mkdir(parents=True)
    tables['export'](output / 'tables', report)
    for name in ('panel.json', 'attempts.tsv', 'cells.tsv'):
        if (output / 'tables' / name).read_bytes() != payloads['data/' + name]:
            raise ValueError('Table replay differs')
    plot = figure['plot'](report)
    try:
        plot.savefig(output / 'shared_threadripper_resources.png', dpi=180)
    finally:
        figure['plt'].close(plot)
    import matplotlib.image as mpimg
    import numpy as np
    original = mpimg.imread(io.BytesIO(payloads['figures/shared_threadripper_resources.png']), format='png')
    rendered = mpimg.imread(output / 'shared_threadripper_resources.png')
    if original.shape != rendered.shape or not np.array_equal(original, rendered):
        raise ValueError('Figure PNG pixels differ')
    if manifest['schema_version'] == 3:
        prose = section_text(payloads)
        if prose != payloads['manuscript/resource_section.md']:
            raise ValueError('Resource section prose replay differs')
        (output / 'resource_section.md').write_bytes(prose)
    verify(directory, manifest_sha)
    result = dict(status='prepared_resource_reporting_replayed', scope=SCOPE,
        manifest_sha256=manifest_sha, exact_table_files=3, equal_png_pixels=True,
        reviewed_attempts=manifest['reviewed_attempts'],
        pre_native_aborted_indices=manifest['pre_native_aborted_indices'],
        native_inference_reproduced=False, raw_measurements_revalidated=False,
        publication_ready=False, matplotlib=matplotlib.__version__, numpy=np.__version__)
    if manifest['schema_version'] == 3:
        result['exact_resource_section'] = True
    (output / 'replay.json').write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    return result


def build(repo, snapshot, figure, output, section=None):
    repo, snapshot, figure = [Path(p).resolve(strict=True) for p in (repo, snapshot, figure)]
    output = Path(output).absolute()
    output = output.parent.resolve() / output.name
    archive = output.with_name(output.name + '.tar.gz')
    if any(p.exists() or p.is_symlink() for p in (output, archive)):
        raise FileExistsError(output)
    commit = subprocess.check_output(['git', '-C', str(repo), 'rev-parse', 'HEAD'], text=True).strip()
    sources = {**{'sources/' + k: v for k, v in SOURCES.items()},
        'component.py': str(Path(__file__).resolve().relative_to(repo)),
        CORE: 'benchmark_tools/results/' + CORE, 'LICENSE.md': 'LICENSE.md'}
    if section is not None:
        sources['sources/section.py'] = SECTION_SOURCE
    payloads = {}
    for name, path in sources.items():
        raw = subprocess.check_output(['git', '-C', str(repo), 'show', commit + ':' + path])
        if raw != (repo / path).read_bytes():
            raise ValueError('Uncommitted reporting source: ' + path)
        payloads[name] = raw
    for name in ('panel.json', 'attempts.tsv', 'cells.tsv'):
        payloads['data/' + name] = (snapshot / name).read_bytes()
    for name in ('manifest.json', *['shared_threadripper_resources.' + e for e in ('png', 'pdf', 'svg')]):
        payloads['figures/' + name] = (figure / name).read_bytes()
    figure_manifest = json.loads(payloads['figures/manifest.json'])
    payloads['requirements.txt'] = ('matplotlib==' + figure_manifest['matplotlib'] + '\n').encode()
    payloads['README.md'] = (
        SCOPE + '\n\nUse component.py verify/replay with externally retained bundle.json SHA-256.\n'
        'No original evidence paths are opened; metadata paths are provenance only.\n'
        'Replay preserves excluded measurements, pre-native aborts and incomplete repeats.\n'
        'This is not raw-accounting/native reproduction, isolated speed evidence or publication readiness.\n').encode()
    report, _, _ = reporting(payloads)
    if section is not None:
        section = Path(section).resolve(strict=True)
        for name in ('manifest.json', 'resource_section.md'):
            payloads['manuscript/' + name] = (section / name).read_bytes()
        check_section(payloads)
        payloads['README.md'] += b'Prose replay includes resource_section.md from the same reviewed snapshot.\n'
    manifest = dict(schema_version=2 if section is None else 3, scope=SCOPE, source_commit=commit,
        **{k: report[k] for k in ('reviewed_attempts', 'resource_reviewed_attempts',
            'eligible_attempts', 'excluded_attempts', 'pre_native_aborted_indices',
            'all_planned_attempts_reviewed')},
        native_inference_reproduced=False, raw_measurements_revalidated=False,
        publication_ready=False,
        files=[dict(path=k, **identity(v)) for k, v in sorted(payloads.items())])
    output.mkdir(parents=True)
    for name, raw in payloads.items():
        target = output / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(raw)
    (output / 'bundle.json').write_text(json.dumps(manifest, indent=2, sort_keys=True) + '\n')
    digest = record(output / 'bundle.json')['sha256']
    verify(output, digest)
    with archive.open('xb') as stream, tarfile.open(fileobj=stream, mode='w:gz') as handle:
        for name in sorted(set(payloads) | {'bundle.json'}):
            raw = (output / name).read_bytes()
            member = tarfile.TarInfo(name)
            member.size, member.mode = len(raw), 0o644
            handle.addfile(member, io.BytesIO(raw))
    return dict(status='prepared_resource_reporting_component_built',
                manifest=record(output / 'bundle.json'), archive=record(archive),
                files=len(payloads), scope=SCOPE, publication_ready=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    build_parser = sub.add_parser('build')
    for name in ('repo', 'snapshot', 'figure', 'output'):
        build_parser.add_argument('--' + name, type=Path, required=True)
    build_parser.add_argument('--section', type=Path)
    for action in ('verify', 'replay'):
        child = sub.add_parser(action)
        child.add_argument('--directory', type=Path, required=True)
        child.add_argument('--manifest-sha256', required=True)
        if action == 'replay':
            child.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.action == 'build':
        result = build(args.repo, args.snapshot, args.figure, args.output, args.section)
    elif args.action == 'replay':
        result = replay(args.directory, args.manifest_sha256, args.output)
    else:
        manifest, _ = verify(args.directory, args.manifest_sha256)
        result = dict(status='prepared_resource_reporting_component_verified',
                      files=len(manifest['files']), scope=SCOPE, publication_ready=False)
    print(json.dumps(result, indent=2))
