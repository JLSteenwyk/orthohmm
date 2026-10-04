"""Generate a resource manuscript section from the independently reviewed panel."""
import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from benchmark_tools.results import plot_shared_prepared_resources_20261004 as plotter

executor = plotter.executor
METRICS = (
    ('wall_seconds', 'Native command wall (seconds)', 1.0, 4),
    ('cpu_seconds', 'Task-subtree CPU bracket (CPU-seconds)', 1.0, 4),
    ('peak_memory_bytes', 'Native-step lifetime peak (GiB)', 1024**3, 3),
)
LABELS = {method: label for method, label, _, _ in plotter.METHODS}


def render(report):
    plotter.validate(report)
    complete = report['all_planned_attempts_reviewed']
    title = '### Shared-Host Resource Panel' + ('' if complete else ' (Interim)')
    measured = report['resource_reviewed_attempts']
    eligible = report['eligible_attempts']
    cells = sum(c['eligible_repeats'] == 3 for c in report['cells'])
    lines = [title, '',
        f"The snapshot retains {report['reviewed_attempts']} of 27 reviewed attempts, "
        f'{measured} with measured native resources and {eligible} eligible observations. '
        f'{cells} method/size cells have three eligible repeats.', '',
        'These are shared-host matched-resource observations with 32 native CPUs and a '
        '128-GiB limit per run, not isolated comparative timing. Contention distortion '
        'is unknown and potentially method dependent. No background overhead is subtracted.', '',
        'Tables report median [minimum, maximum] only for cells with three eligible repeats. '
        'An unavailable summary is not zero or a failed native inference. Ranges are not '
        'confidence intervals. Eligibility counts retain measured exclusions and pre-native '
        'aborts rather than selecting the fastest attempts.', '']
    for key, heading, divisor, digits in METRICS:
        lines += [f'#### {heading}', '',
            '| Method | Proteomes | Eligible/planned | Median [minimum, maximum] |',
            '| --- | --- | --- | --- |']
        for cell in report['cells']:
            summary = cell['resources'][key]
            value = 'Unavailable'
            if summary is not None:
                numbers = [f"{summary[k]/divisor:.{digits}f}" for k in ('median', 'minimum', 'maximum')]
                value = f'{numbers[0]} [{numbers[1]}, {numbers[2]}]'
            lines.append(f"| {LABELS[cell['method']]} | {cell['proteomes']} | "
                f"{cell['eligible_repeats']}/{cell['planned_repeats']} | {value} |")
        lines.append('')
    excluded = ', '.join(str(i) for i in report['excluded_attempts']) or 'none'
    aborted = ', '.join(str(i) for i in report['pre_native_aborted_indices']) or 'none'
    pending = ', '.join(str(r['index']) for r in report['runs'] if r['status'] == 'not_yet_reviewed') or 'none'
    foreign = [r['whole_run_maximum_foreign_average_cores'] for r in report['runs']
        if r['resources'] is not None]
    lines += [f'Excluded attempt indices: {excluded}. Pre-native abort indices: {aborted}; '
        'these have no native resource endpoints, not zero measurements. '
        f'Unreviewed indices: {pending}.', '',
        f'Across measured attempts, maximum observed foreign CPU demand ranges from '
        f'{min(foreign):.4f} to {max(foreign):.4f} core equivalents. These are '
        'process-interval observations, not reservations or estimates of causal slowdown.', '',
        'CPU includes the native task-subtree wrapper bracket. Peak memory includes the '
        'native-step launcher and is not pure algorithm RSS. Preparation, conversion and '
        'scoring are outside the native command timer. This single nested taxon series '
        'co-varies proteome count and taxon composition, not taxon-invariant scaling. '
        'Native output validation does not provide new prediction-accuracy evidence. '
        'Historical DGX and earlier shared-host timings are not pooled.', '']
    return '\n'.join(lines)


def build(table, sha256, output):
    output = Path(output)
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    table_ref = executor.record(Path(table).resolve())
    report = executor.read_frozen(Path(table_ref['path']), sha256)
    prose = render(report)
    sources = [executor.record(path) for path in (__file__, plotter.__file__, plotter.tables.__file__,
        plotter.previous.__file__, plotter.tables.previous.__file__)]
    for pin in [table_ref, *report['evidence'], *sources]:
        executor.check(pin)
    output.mkdir()
    section = output / 'resource_section.md'
    with section.open('x') as stream:
        stream.write(prose)
    receipt = dict(status='reviewed_resource_section_generated', table=table_ref,
        sources=sources, section=executor.record(section),
        reviewed_attempts=report['reviewed_attempts'], measured_attempts=report['resource_reviewed_attempts'],
        eligible_attempts=report['eligible_attempts'], excluded_attempts=report['excluded_attempts'],
        all_planned_attempts_reviewed=report['all_planned_attempts_reviewed'],
        native_inference_repeated=False, raw_audit_repeated=False, publication_ready=False,
        scope='Prose and tables from reviewed resource evidence, not new accuracy or isolated speed')
    with (output / 'manifest.json').open('x') as stream:
        json.dump(receipt, stream, indent=2, sort_keys=True)
        stream.write('\n')
    return receipt


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--table', type=Path, required=True)
    parser.add_argument('--sha256', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(build(args.table, args.sha256, args.output), indent=2))
