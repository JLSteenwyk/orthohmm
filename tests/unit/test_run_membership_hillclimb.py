import json
import os
from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmark_tools import run_membership_hillclimb as runner


def original():
    return ['python', 'old/replay.py', '--fasta-directory', 'inputs',
            '--candidate-clusters', 'candidates', '--membership-constraints', 'constraints',
            '--output-directory', 'prior', '--json', 'metrics', '--cpu', '32',
            '--species-tree-mode', 'infer', '--checkpoint-source', 'other_partition']


def test_command_uses_matching_candidate_checkpoint_and_explicit_policy():
    prior = original()
    command = runner.derive_command(prior, Path('snapshot'), Path('new'), 'high_confidence_path', 8)
    assert command[1] == 'snapshot/benchmark_tools/replay_phylogeny.py'
    assert runner.argument(command, '--checkpoint-source') == 'prior'
    assert runner.argument(command, '--output-directory') == 'new/output'
    assert runner.argument(command, '--satellite-membership-rule') == 'high_confidence_path'
    assert runner.argument(command, '--cpu') == '8'
    assert prior == original()


@pytest.mark.parametrize('argv', [original() + ['--cpu', '4'], ['--cpu']])
def test_argument_requires_unique_flag(argv):
    with pytest.raises((ValueError, IndexError)):
        runner.argument(argv, '--cpu')


def test_command_rejects_supplied_species_tree():
    prior = original()
    prior[prior.index('--species-tree-mode') + 1] = 'supplied'
    with pytest.raises(ValueError, match='inferred'):
        runner.derive_command(prior, Path('snapshot'), Path('out'), 'high_confidence_pair', 8)


def test_partition_requires_full_unique_gene_universe(tmp_path):
    before, after = tmp_path / 'before', tmp_path / 'after'
    before.write_text('a b\nc\n')
    after.write_text('a\nb c\n')
    assert runner.validate_partition(before, after)['all_candidate_genes_preserved']
    for text in ('a b\n', 'a b\na c\n', 'a b\nc d\n'):
        after.write_text(text)
        with pytest.raises(ValueError):
            runner.validate_partition(before, after)


def test_golden_comparison_fails_closed(tmp_path):
    path = tmp_path / 'pairs.tsv'
    path.write_text('original\n')
    golden = {'pairs.tsv': runner.record(path)}
    assert runner.compare_golden(tmp_path, golden) == {'pairs.tsv': True}
    path.write_text('changed\n')
    with pytest.raises(ValueError, match='exact retained-output'):
        runner.compare_golden(tmp_path, golden)


def mini_plan(tmp_path, monkeypatch):
    source = tmp_path / 'source'
    tools = source / 'benchmark_tools'
    tools.mkdir(parents=True)
    script = tools / 'run_membership_hillclimb.py'
    script.write_text('frozen runner identity\n')
    monkeypatch.setattr(runner, '__file__', str(script))
    experiment = tmp_path / 'experiment.json'
    experiment.write_text('{}\n')
    candidates = tmp_path / 'candidates'
    candidates.write_text('a b\nc\n')
    golden_dir = tmp_path / 'golden'
    golden_dir.mkdir()
    golden = {}
    for name in (*runner.PRODUCTS, 'species_tree.rooted.nwk'):
        path = golden_dir / name
        path.write_text('golden\n')
        golden[name] = runner.record(path)
    output = tmp_path / 'results'
    output.mkdir()
    cases = []
    for index, policy in enumerate(('high_confidence_pair', 'high_confidence_path')):
        arm = output / 'OrthoBench' / str(index)
        cases.append({'label': str(index), 'index': index, 'dataset': 'OrthoBench',
                      'cpm_resolution': .1, 'membership_rule': policy, 'status': 'prepared',
                      'output': str(arm), 'command': ['fake', '--output-directory', str(arm / 'output'),
                                                    '--candidate-clusters', str(candidates)],
                      'inputs': [runner.record(candidates)], 'golden': golden, 'scoring_inputs': []})
    cases.append({'label': 'missing', 'dataset': 'OrthoBench', 'cpm_resolution': .08,
                  'membership_rule': 'high_confidence_pair', 'status': 'unavailable', 'reason': 'Missing checkpoint'})
    plan = {'repo': str(tmp_path), 'source_root': str(source), 'source_revision': 'test',
            'plan': runner.record(experiment), 'sources': [runner.record(script)],
            'output': str(output), 'scope': 'test', 'cases': cases}
    plan['tool_path'] = os.environ['PATH']
    path = output / 'execution_plan.json'
    runner.write_json(path, plan)
    return path


def fake_replay(command, **kwargs):
    output = Path(runner.argument(command, '--output-directory'))
    phylogeny = output / 'orthohmm_phylogeny'
    working = output / 'orthohmm_working_res'
    phylogeny.mkdir(parents=True)
    working.mkdir()
    for name in (*runner.PRODUCTS, 'species_tree.rooted.nwk'):
        (phylogeny / name).write_text('golden\n')
    (working / 'orthohmm_edges_clustered.txt').write_text('a\nb c\n')
    runner.write_json(phylogeny / 'reconciliation_summary.json', {'ortholog_pairs': 1})
    runner.write_json(phylogeny / 'provenance_manifest.json', {'membership_reconciliation': {}})
    runner.write_json(output.parent / 'metrics.json', {'status': 'complete'})
    return SimpleNamespace(returncode=0)


def test_runner_scores_after_validating_control_and_retains_missing_cells(tmp_path, monkeypatch):
    path = mini_plan(tmp_path, monkeypatch)
    monkeypatch.setattr(runner.subprocess, 'run', fake_replay)
    monkeypatch.setattr(runner, 'score_orthobench', lambda *_: {'scores': {'OrthoBench F1': .75}})
    report = runner.run(path, 'OrthoBench')
    assert [r['status'] for r in report['arms']] == [
        'scored_development_result', 'scored_development_result', 'unavailable']
    assert report['arms'][0]['golden_matches']
    assert all(r['inputs_reverified_unchanged'] for r in report['arms'][:2])
    assert not report['defaults_promoted']
    assert not report['all_metric_screening_complete']


def test_runner_never_interprets_path_arm_after_control_mismatch(tmp_path, monkeypatch):
    path = mini_plan(tmp_path, monkeypatch)
    calls = []
    def changed(command, **kwargs):
        calls.append(command)
        result = fake_replay(command, **kwargs)
        output = Path(runner.argument(command, '--output-directory'))
        (output / 'orthohmm_phylogeny' / runner.PRODUCTS[0]).write_text('wrong\n')
        return result
    monkeypatch.setattr(runner.subprocess, 'run', changed)
    monkeypatch.setattr(runner, 'score_orthobench', lambda *_: pytest.fail('Scoring must not run'))
    report = runner.run(path, 'OrthoBench')
    assert len(calls) == 1
    assert report['arms'][0]['status'] == report['arms'][1]['status'] == 'failed'
    assert 'exact retained-output' in report['arms'][0]['error']
    assert 'Matching direct control failed' in report['arms'][1]['error']


def test_failed_scoring_keeps_native_evidence_and_does_not_abort_other_arms(tmp_path, monkeypatch):
    path = mini_plan(tmp_path, monkeypatch)
    monkeypatch.setattr(runner.subprocess, 'run', fake_replay)
    def failed(*args):
        raise RuntimeError('scorer failed')
    monkeypatch.setattr(runner, 'score_orthobench', failed)
    report = runner.run(path, 'OrthoBench')
    assert report['status'] == 'finished_with_explicit_gaps'
    for row in report['arms'][:2]:
        assert row['status'] == 'failed'
        assert row['error'] == 'scorer failed'
        assert row['outputs']
        assert row['native_summary']['ortholog_pairs'] == 1
