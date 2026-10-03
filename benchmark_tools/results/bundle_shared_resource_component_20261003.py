"""Archive and replay shared-host resource reporting, not native measurements."""
import argparse
import ast
import csv
import hashlib
import io
import json
import math
from pathlib import Path, PurePosixPath
import re
from statistics import median
import subprocess
import sys
import tarfile
from types import SimpleNamespace

SCOPE = 'shared_host_resource_reporting_not_raw_measurement_reproduction'
SOURCES = {
    'plan.py': 'benchmark_tools/prepare_scaling_inputs.py',
    'scopes.py': 'benchmark_tools/derive_threadripper_resources.py',
    'tables.py': 'benchmark_tools/results/report_shared_threadripper_panel_20261003.py',
    'figure.py': 'benchmark_tools/results/plot_shared_threadripper_resources_20261003.py',
}
EXPECTED = {'component.py', 'LICENSE.md', 'requirements.txt', 'README.md',
    'data/panel.json', 'data/attempts.tsv', 'data/cells.tsv',
    'figures/shared_threadripper_resources.png', 'figures/shared_threadripper_resources.pdf',
    'figures/shared_threadripper_resources.svg', 'figures/manifest.json',
    *['sources/' + name for name in SOURCES]}


def identity(data):
    return dict(bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def record(path):
    path = Path(path).resolve()
    return dict(path=str(path), **identity(path.read_bytes()))


def safe_name(name):
    path = PurePosixPath(name)
    if not name or path.is_absolute() or '..' in path.parts or str(path) != name or '\\' in name:
        raise ValueError('Unsafe component member')
    return name


def checked(raw, pin):
    if identity(raw) != {key: pin[key] for key in ('bytes', 'sha256')}:
        raise ValueError('Component bytes differ')


def pure(source, names, namespace):
    nodes = []
    found = set()
    for node in ast.parse(source).body:
        name = node.name if isinstance(node, ast.FunctionDef) else (
            node.targets[0].id if isinstance(node, ast.Assign) and len(node.targets) == 1
            and isinstance(node.targets[0], ast.Name) else None)
        if name in names:
            if name in found:
                raise ValueError('Duplicate pure reporting entry point')
            found.add(name)
            nodes.append(node)
    if found != set(names):
        raise ValueError('Missing pure reporting entry point')
    exec(compile(ast.Module(body=nodes, type_ignores=[]), '<bundled-reporting-source>', 'exec'), namespace)
    return namespace


def reporting(payloads, plotting=False):
    plan = pure(payloads['sources/plan.py'], {'METHODS', 'SIZES', 'planned_runs'}, {})
    scopes = pure(payloads['sources/scopes.py'], {'SCOPES'}, {})['SCOPES']
    namespace = dict(math=math, median=median, json=json, csv=csv,
        executor=SimpleNamespace(record=record), **{key: plan[key] for key in ('METHODS', 'SIZES', 'planned_runs')},
        SCOPES=scopes)
    tables = pure(payloads['sources/tables.py'], {'METRICS', 'IDENTITY', 'summarize', 'export'}, namespace)
    import_namespace = dict(math=math, planned_runs=plan['planned_runs'],
        tables=SimpleNamespace(summarize=tables['summarize'], SCOPES=scopes),
        executor=SimpleNamespace(SHARED_SCOPE='shared_host_matched_resources'))
    if plotting:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
        import_namespace.update(plt=plt, Line2D=Line2D)
    figure = pure(payloads['sources/figure.py'], {'METHODS', 'METRICS', 'validate', 'plot'}, import_namespace)
    report = json.loads(payloads['data/panel.json'])
    figure['validate'](report)
    return report, tables, figure


def verify(directory, manifest_sha):
    directory = Path(directory).resolve()
    if (directory / 'bundle.json').is_symlink():
        raise ValueError('Manifest must be a regular local file')
    raw = (directory / 'bundle.json').read_bytes()
    if identity(raw)['sha256'] != manifest_sha:
        raise ValueError('Externally pinned manifest differs')
    manifest = json.loads(raw)
    if (type(manifest['schema_version']) is not int or manifest['schema_version'] != 1
            or not re.fullmatch('[0-9a-f]{40}', manifest['source_commit'])
            or manifest['scope'] != SCOPE or manifest['native_inference_reproduced'] is not False
            or manifest['raw_measurements_revalidated'] is not False or manifest['publication_ready'] is not False):
        raise ValueError('Wrong reporting-only scope')
    payloads = {}
    for pin in manifest['files']:
        name = safe_name(pin['path'])
        path = directory / name
        if name not in EXPECTED or name in payloads or path.is_symlink() or not path.is_file():
            raise ValueError('Wrong component inventory or member type')
        payloads[name] = path.read_bytes()
        checked(payloads[name], pin)
    actual = {path.relative_to(directory).as_posix() for path in directory.rglob('*') if path.is_file() or path.is_symlink()}
    if set(payloads) != EXPECTED or actual != EXPECTED | {'bundle.json'}:
        raise ValueError('Missing or extra component member')
    if payloads['component.py'] != Path(__file__).read_bytes():
        raise ValueError('Run the exact bundled reporting reader')
    report, _, _ = reporting(payloads)
    checked(payloads['sources/tables.py'], report['source'])
    if (manifest['reviewed_attempts'], manifest['eligible_attempts'], manifest['excluded_attempts']) != (
            report['reviewed_attempts'], report['eligible_attempts'], report['excluded_attempts']):
        raise ValueError('Manifest coverage differs from reviewed snapshot')
    figure = json.loads(payloads['figures/manifest.json'])
    checked(payloads['data/panel.json'], figure['source_results'])
    checked(payloads['sources/figure.py'], figure['plotter'])
    for pin in figure['outputs']:
        checked(payloads['figures/' + Path(pin['path']).name], pin)
    if payloads['requirements.txt'].decode() != 'matplotlib==' + figure['matplotlib'] + '\n':
        raise ValueError('Reporting dependency version differs')
    return manifest, payloads


def replay(directory, manifest_sha, output):
    directory, output = Path(directory).resolve(), Path(output)
    output = output.parent.resolve() / output.name
    if output.exists() or output.is_symlink() or output.is_relative_to(directory):
        raise FileExistsError('Require a fresh output outside the component')
    manifest, payloads = verify(directory, manifest_sha)
    import matplotlib
    if payloads['requirements.txt'].decode() != 'matplotlib==' + matplotlib.__version__ + '\n':
        raise ValueError('Use the recorded reporting environment')
    report, tables, figure = reporting(payloads, plotting=True)
    output.mkdir(parents=True)
    tables['export'](output / 'tables', report)
    for name in ('panel.json', 'attempts.tsv', 'cells.tsv'):
        if (output / 'tables' / name).read_bytes() != payloads['data/' + name]:
            raise ValueError('Resource table replay differs')
    plot = figure['plot'](report)
    try:
        for extension in ('png', 'pdf', 'svg'):
            plot.savefig(output / ('shared_threadripper_resources.' + extension), dpi=180)
    finally:
        figure['plt'].close(plot)
    import matplotlib.image as mpimg
    import numpy as np
    expected = mpimg.imread(io.BytesIO(payloads['figures/shared_threadripper_resources.png']), format='png')
    actual = mpimg.imread(output / 'shared_threadripper_resources.png')
    if expected.shape != actual.shape or not np.array_equal(expected, actual):
        raise ValueError('Resource figure pixels differ')
    verify(directory, manifest_sha)
    result = dict(status='shared_resource_reporting_replayed', scope=SCOPE,
        bundle_manifest_sha256=manifest_sha, reviewed_attempts=manifest['reviewed_attempts'],
        eligible_attempts=manifest['eligible_attempts'], excluded_attempts=manifest['excluded_attempts'],
        exact_table_files=3, equal_png_pixels=True, raw_measurements_revalidated=False,
        native_inference_reproduced=False, publication_ready=False,
        matplotlib=matplotlib.__version__, numpy=np.__version__)
    (output / 'replay.json').write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    return result


def build(repo, snapshot, figure, output):
    repo, snapshot, figure = map(lambda path: Path(path).resolve(), (repo, snapshot, figure))
    output = Path(output)
    output = output.parent.resolve() / output.name
    archive = output.parent / (output.name + '.tar.gz')
    if output.exists() or output.is_symlink() or archive.exists() or archive.is_symlink():
        raise FileExistsError(output)
    commit = subprocess.check_output(['git', '-C', str(repo), 'rev-parse', 'HEAD'], text=True).strip()
    source = Path(__file__).resolve()
    source_names = dict(SOURCES, **{'component.py': str(source.relative_to(repo)), 'LICENSE.md': 'LICENSE.md'})
    payloads = {}
    for target, name in source_names.items():
        raw = subprocess.check_output(['git', '-C', str(repo), 'show', commit + ':' + name])
        if raw != (repo / name).read_bytes():
            raise ValueError('Reporting source differs from committed bytes')
        payloads[target if target in {'component.py', 'LICENSE.md'} else 'sources/' + target] = raw
    for name in ('panel.json', 'attempts.tsv', 'cells.tsv'):
        payloads['data/' + name] = (snapshot / name).read_bytes()
    for name in ('manifest.json', 'shared_threadripper_resources.png', 'shared_threadripper_resources.pdf', 'shared_threadripper_resources.svg'):
        payloads['figures/' + name] = (figure / name).read_bytes()
    figure_manifest = json.loads(payloads['figures/manifest.json'])
    payloads['requirements.txt'] = ('matplotlib==' + figure_manifest['matplotlib'] + '\n').encode()
    payloads['README.md'] = (SCOPE + '\n\nReplay with component.py and the externally supplied bundle.json SHA-256.\n'
        'Original absolute evidence paths are provenance only and are not read.\n'
        'Tables and PNG pixels reproduce; raw cgroup/runtime/native evidence is not revalidated.\n'
        'Excluded attempts and incomplete repeats remain explicit. No isolated speedup or publication readiness.\n').encode()
    report, _, _ = reporting(payloads)
    manifest = dict(schema_version=1, scope=SCOPE, source_commit=commit,
        reviewed_attempts=report['reviewed_attempts'], eligible_attempts=report['eligible_attempts'],
        excluded_attempts=report['excluded_attempts'], native_inference_reproduced=False,
        raw_measurements_revalidated=False, publication_ready=False,
        files=[dict(path=name, **identity(raw)) for name, raw in sorted(payloads.items())])
    output.mkdir(parents=True)
    for name, raw in payloads.items():
        path = output / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(raw)
    manifest_path = output / 'bundle.json'
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + '\n')
    digest = record(manifest_path)['sha256']
    verify(output, digest)
    with archive.open('xb') as stream, tarfile.open(fileobj=stream, mode='w:gz') as handle:
        for name in sorted(EXPECTED | {'bundle.json'}):
            raw = (output / name).read_bytes()
            member = tarfile.TarInfo(name)
            member.size, member.mode = len(raw), 0o644
            handle.addfile(member, io.BytesIO(raw))
    return dict(status='shared_resource_reporting_component_built', scope=SCOPE,
        manifest=record(manifest_path), archive=record(archive), files=len(payloads),
        native_inference_reproduced=False, raw_measurements_revalidated=False, publication_ready=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('build', 'verify', 'replay'))
    parser.add_argument('--repo', type=Path)
    parser.add_argument('--snapshot', type=Path)
    parser.add_argument('--figure', type=Path)
    parser.add_argument('--directory', type=Path)
    parser.add_argument('--manifest-sha256')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    required = ('repo', 'snapshot', 'figure', 'output') if args.action == 'build' else ('directory', 'manifest_sha256')
    if args.action == 'replay':
        required += ('output',)
    if any(getattr(args, key) is None for key in required):
        parser.error('Missing required ' + args.action + ' arguments: ' + ', '.join(required))
    if args.action == 'build':
        result = build(args.repo, args.snapshot, args.figure, args.output)
    elif args.action == 'verify':
        manifest, _ = verify(args.directory, args.manifest_sha256)
        result = dict(status='shared_resource_reporting_component_verified', files=len(manifest['files']), scope=SCOPE,
            native_inference_reproduced=False, raw_measurements_revalidated=False, publication_ready=False)
    else:
        result = replay(args.directory, args.manifest_sha256, args.output)
    print(json.dumps(result, indent=2))
