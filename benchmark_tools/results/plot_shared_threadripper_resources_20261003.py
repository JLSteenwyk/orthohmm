"""Plot reviewed shared-host attempts without hiding exclusions or missing repeats."""
import argparse
import json
import math
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from benchmark_tools.results import report_shared_threadripper_panel_20261003 as tables
from benchmark_tools.prepare_scaling_inputs import planned_runs
from benchmark_tools import run_threadripper_scaling as executor

METHODS = (
    ('orthohmm_high_sensitivity', 'OrthoHMM high sensitivity', '#007d83', 'o'),
    ('orthohmm_satellite_v2', 'OrthoHMM inferred phylogeny', '#a66b0b', 's'),
    ('orthofinder_3_1_5_full', 'OrthoFinder 3.1.5 full', '#755297', '^'))
METRICS = (
    ('wall_seconds', 60., 'A  Native elapsed time', 'Minutes'),
    ('cpu_seconds', 3600., 'B  Native CPU use', 'CPU-hours'),
    ('peak_memory_bytes', 1024.**3, 'C  Native-step memory peak', 'GiB (not process RSS)'))


def validate(report):
    expected = planned_runs()
    if len(report['runs']) != len(expected):
        raise ValueError('Changed planned attempt inventory')
    reviews = []
    for row, identity in zip(report['runs'], expected):
        if any(type(row[key]) is not type(value) or row[key] != value for key, value in identity.items()):
            raise ValueError('Changed frozen figure identity')
        if row['status'] == 'not_yet_reviewed':
            continue
        eligible = row['comparative_timing_eligible']
        if type(eligible) is not bool:
            raise ValueError('Reviewed row lacks explicit eligibility')
        reviews.append(dict(identity, job_id=row['job_id'], execution_scope=executor.SHARED_SCOPE,
            uncontended_timing=False, scientific_timings_admitted=False,
            resource_scopes=tables.SCOPES, primary_resources_replayed=True,
            shared_host_resources_reviewed=eligible, original_environment_protocol_passed=eligible,
            status='shared_attempt_independently_reviewed' if eligible else 'shared_attempt_reviewed_with_environment_failure',
            scheduler_state='COMPLETED' if eligible else 'FAILED', scheduler_exit_code='0:0' if eligible else '1:0',
            resources=row['resources'], preflight_foreign_average_cores=row['preflight_foreign_average_cores'],
            whole_run_maximum_foreign_average_cores=row['whole_run_maximum_foreign_average_cores']))
        for key in ('preflight_foreign_average_cores', 'whole_run_maximum_foreign_average_cores'):
            value = row[key]
            if type(value) not in (int, float) or not math.isfinite(value) or value < 0:
                raise ValueError('Invalid contention annotation')
    rebuilt = tables.summarize(expected, reviews)
    for key, value in rebuilt.items():
        if report.get(key) != value or type(report[key]) is not type(value):
            raise ValueError('Snapshot differs from recomputed table: ' + key)
    return [row for row in report['runs'] if row['status'] != 'not_yet_reviewed']


def plot(report):
    rows = validate(report)
    fig, axes = plt.subplots(1, 3, figsize=(15, 9))
    fig.subplots_adjust(left=.065, right=.985, bottom=.44, top=.76, wspace=.32)
    title = 'Threadripper shared-host resource observations'
    if not report['all_planned_attempts_reviewed']:
        title += ' (partial)'
    fig.suptitle(title, x=.045, ha='left', y=.975, fontsize=18)
    fig.text(.045, .924, f"{report['reviewed_attempts']}/27 attempts reviewed | "
        f"{report['eligible_attempts']} eligible | {len(report['excluded_attempts'])} excluded | "
        '32 native CPUs and 128 GiB per run', fontsize=11)
    for ax, (metric, divisor, title, ylabel) in zip(axes, METRICS):
        for method_index, (method, _, color, marker) in enumerate(METHODS):
            for position, size in enumerate((4, 8, 12)):
                x = position + (method_index - 1) * .22
                selected = [row for row in rows if row['method'] == method and row['proteomes'] == size]
                for row in selected:
                    eligible = row['comparative_timing_eligible']
                    ax.scatter([x + (row['repeat'] - 1) * .04], [row['resources'][metric] / divisor],
                        s=34, marker=marker if eligible else 'x', color=color if eligible else '#777777',
                        linewidths=1.4, zorder=4)
                cell = next(cell for cell in report['cells'] if cell['method'] == method and cell['proteomes'] == size)
                summary = cell['resources'][metric]
                if summary is not None:
                    ax.plot([x, x], [summary['minimum'] / divisor, summary['maximum'] / divisor],
                        color=color, linewidth=1.5, zorder=2)
                    ax.plot([x - .045, x + .045], [summary['median'] / divisor] * 2,
                        color=color, linewidth=2.5, zorder=3)
        ax.set_title(title, loc='left', fontsize=12, pad=12)
        ax.set_xticks(range(3), ['4', '8', '12'])
        ax.set_xlabel('Complete proteomes')
        ax.set_ylabel(ylabel)
        ax.set_xlim(-.5, 2.5)
        maximum = max((row['resources'][metric] / divisor for row in rows), default=1.)
        ax.set_ylim(0, maximum * 1.10)
        ax.grid(axis='y', color='#dddddd', linewidth=.6)
        ax.set_axisbelow(True)
        ax.spines[['top', 'right']].set_visible(False)
    handles = [Line2D([], [], color=color, marker=marker, linestyle='none', label=label)
        for _, label, color, marker in METHODS]
    handles.append(Line2D([], [], color='#777777', marker='x', linestyle='none', label='Excluded raw attempt'))
    fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(.04, .875),
        ncol=4, frameon=False, fontsize=10)
    notes = [
        'Symbols: individual reviewed attempts. Gray crosses: excluded raw values, never pooled into eligible summaries.',
        'Median ticks and range bars require three eligible repeats; observed ranges are not confidence intervals.',
        'Unreviewed attempts are absent, not zero. Shared-host contention is unknown and potentially method dependent.',
        'No uncontended timing, causal speedup, fastest-repeat selection, overhead subtraction or historical-time pooling.',
        'CPU: native task-subtree bracket including wrapper. Peak: native-step lifetime including launcher; not pure algorithm RSS.',
        'One nested taxon series: 73,266 / 165,168 / 251,378 proteins; taxon composition and proteome count co-vary.',
    ]
    for y, note in zip((.37, .335, .30, .265, .23, .195), notes):
        fig.text(.045, y, note, fontsize=10)
    for number, (method, label, _, _) in enumerate(METHODS):
        cells = [cell for cell in report['cells'] if cell['method'] == method]
        coverage = '; '.join(f"{cell['proteomes']} proteomes: {cell['eligible_repeats']}/3 eligible, "
            f"{len(cell['excluded_indices'])} excluded" for cell in cells)
        fig.text(.045, .14 - number * .035, label + ' | ' + coverage, fontsize=10)
    return fig


def export(report_path, sha256, output):
    report_path, output = Path(report_path).resolve(), Path(output).resolve()
    if output.exists():
        raise FileExistsError(output)
    report_ref = executor.record(report_path)
    report = executor.read_frozen(report_path, sha256)
    refs = [report_ref, executor.record(__file__), executor.record(tables.__file__),
        report['source'], report['plan'], *report['evidence']]
    for ref in refs:
        executor.check(ref)
    figure = plot(report)
    output.mkdir(parents=True)
    try:
        outputs = []
        for extension in ('png', 'pdf', 'svg'):
            path = output / ('shared_threadripper_resources.' + extension)
            figure.savefig(path, dpi=180)
            outputs.append(executor.record(path))
        for ref in refs:
            executor.check(ref)
        manifest = dict(source_results=report_ref, plotter=refs[1], inputs=refs,
            outputs=outputs, matplotlib=matplotlib.__version__,
            reviewed_attempts=report['reviewed_attempts'], eligible_attempts=report['eligible_attempts'],
            excluded_indices=report['excluded_attempts'], all_planned_attempts_reviewed=report['all_planned_attempts_reviewed'],
            primary_scopes=report['primary_scopes'], execution_scope=executor.SHARED_SCOPE,
            uncontended_timing=False, scientific_timings_admitted=False, publication_ready=False)
        with (output / 'manifest.json').open('x') as stream:
            json.dump(manifest, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write('\n')
        return manifest
    finally:
        plt.close(figure)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, required=True)
    parser.add_argument('--sha256', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = export(args.results, args.sha256, args.output)
    print(json.dumps({key: result[key] for key in ('reviewed_attempts', 'eligible_attempts', 'excluded_indices', 'outputs')}, indent=2))
