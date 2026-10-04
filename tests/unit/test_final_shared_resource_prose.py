import json
from pathlib import Path
from statistics import median

from benchmark_tools.results import render_shared_resource_section_20261004 as section

RESULTS = Path(__file__).resolve().parents[2] / 'benchmark_tools/results'
TABLE = RESULTS / 'threadripper_shared_panel_snapshot_20261004_v27/panel.json'
PROSE = RESULTS / 'threadripper_shared_resource_section_20261004_v27'


def test_final_panel_retains_previous_rows_and_missing_repeats():
    report = json.loads(TABLE.read_text())
    previous = json.loads((RESULTS / 'threadripper_shared_panel_snapshot_20261004_v26/panel.json').read_text())
    assert report['runs'][:26] == previous['runs'][:26]
    assert (report['reviewed_attempts'], report['resource_reviewed_attempts'],
            report['eligible_attempts']) == (27, 25, 24)
    assert report['all_planned_attempts_reviewed'] is True
    assert report['all_cells_have_three_eligible_repeats'] is False
    assert report['excluded_attempts'] == [0, 17, 20]
    assert report['pre_native_aborted_indices'] == [17, 20]
    assert not any(row['status'] == 'not_yet_reviewed' for row in report['runs'])
    assert report['publication_ready'] is False
    assert sum(cell['eligible_repeats'] == 3 for cell in report['cells']) == 6


def test_final_prose_values_match_eligible_raw_rows_and_pinned_sources():
    report = json.loads(TABLE.read_text())
    receipt = json.loads((PROSE / 'manifest.json').read_text())
    text = (PROSE / 'resource_section.md').read_text()
    assert receipt['table'] == section.executor.record(TABLE)
    for pin in [receipt['table'], receipt['section'], *receipt['sources']]:
        section.executor.check(pin)
    assert receipt['all_planned_attempts_reviewed'] is True
    assert receipt['native_inference_repeated'] is False
    assert receipt['raw_audit_repeated'] is False
    assert text == section.render(report)
    assert '(Interim)' not in text
    assert 'Unreviewed indices: none' in text
    assert 'unknown and potentially method dependent' in text
    assert 'Historical DGX and earlier shared-host timings are not pooled' in text
    rows = [line.split('|')[1:-1] for line in text.splitlines()
            if line.startswith('| Ortho')]
    assert len(rows) == 27
    unavailable = 0
    for metric_index, (metric, _, divisor, digits) in enumerate(section.METRICS):
        for cell_index, cell in enumerate(report['cells']):
            row = [value.strip() for value in rows[metric_index * 9 + cell_index]]
            eligible = [attempt for attempt in report['runs']
                        if attempt['comparative_timing_eligible'] is True
                        and (attempt['method'], attempt['proteomes'])
                        == (cell['method'], cell['proteomes'])]
            assert row[:3] == [section.LABELS[cell['method']],
                               str(cell['proteomes']), f'{len(eligible)}/3']
            if len(eligible) != 3:
                assert row[3] == 'Unavailable'
                unavailable += 1
                continue
            values = [attempt['resources'][metric] for attempt in eligible]
            numbers = [f'{value / divisor:.{digits}f}'
                       for value in (median(values), min(values), max(values))]
            assert row[3] == f'{numbers[0]} [{numbers[1]}, {numbers[2]}]'
    assert unavailable == 9
