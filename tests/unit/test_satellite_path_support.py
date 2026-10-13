from dataclasses import replace

import pytest

from benchmark_tools.replay_phylogeny import build_parser, main
from orthohmm.phylogeny import PhylogenyConfig, PhylogenyConfigurationError
from orthohmm.phylogeny_pipeline import (
    FamilyOutcome,
    _family_tree_input_hash,
    apply_satellite_membership_constraints,
)


def family(pairs=(('a', 'b'), ('b', 'c')), groups=(('a', 'b', 'c'),), confidence=None):
    return FamilyOutcome(
        family_id='F0', genes=('a', 'b', 'c'), root_groups=groups,
        ortholog_pairs=pairs,
        ortholog_pair_confidence=confidence or tuple((*pair, 'high') for pair in pairs),
        reconciliation=None, checkpoint_hit=False,
    )


CONSTRAINT = [{'source_genes': ['a'], 'target_genes': ['c']}]


def test_path_support_retains_existing_predictions_not_transitive_pairs():
    original = family()
    direct, direct_audit = apply_satellite_membership_constraints([original], CONSTRAINT)
    path, audit = apply_satellite_membership_constraints(
        [original], CONSTRAINT, policy='high_confidence_path'
    )
    assert direct[0].ortholog_pairs == (('b', 'c'),)
    assert direct_audit['detached_constraints'] == 1
    assert path == [original]
    assert ('a', 'c') not in path[0].ortholog_pairs
    assert audit['path_supported_constraints'] == 1
    assert audit['ortholog_pairs_removed'] == 0


@pytest.mark.parametrize('confidence', ['medium', 'low'])
def test_path_support_rejects_weaker_links(confidence):
    original = family(confidence=(('a', 'b', 'high'), ('b', 'c', confidence)))
    path, audit = apply_satellite_membership_constraints(
        [original], CONSTRAINT, policy='high_confidence_path'
    )
    assert audit['detached_constraints'] == 1
    assert audit['path_supported_constraints'] == 0
    assert ('a', 'b') not in path[0].ortholog_pairs


def test_support_path_cannot_cross_native_root_groups():
    original = family(groups=(('a', 'b'), ('c',)))
    _, audit = apply_satellite_membership_constraints(
        [original], CONSTRAINT, policy='high_confidence_path'
    )
    assert audit['path_supported_constraints'] == 0
    assert audit['detached_constraints'] == 1


def test_unpredicted_confidence_row_cannot_supply_indirect_support():
    original = family(pairs=(('a', 'b'),),
                      confidence=(('a', 'b', 'high'), ('b', 'c', 'high')))
    _, audit = apply_satellite_membership_constraints(
        [original], CONSTRAINT, policy='high_confidence_path'
    )
    assert audit['path_supported_constraints'] == 0


def test_disconnected_satellites_do_not_support_each_other():
    original = family(pairs=(('a', 'b'),))
    constraints = [
        {'source_genes': ['a'], 'target_genes': ['c']},
        {'source_genes': ['b'], 'target_genes': ['c']},
    ]
    forward, audit = apply_satellite_membership_constraints(
        [original], constraints, policy='high_confidence_path'
    )
    backward, _ = apply_satellite_membership_constraints(
        [original], list(reversed(constraints)), policy='high_confidence_path'
    )
    assert audit['detached_constraints'] == 2
    assert set(forward[0].root_groups) == set(backward[0].root_groups)
    assert forward[0].ortholog_pairs == backward[0].ortholog_pairs == ()


def test_direct_support_and_default_audit_remain_unchanged():
    original = family(pairs=(('a', 'c'),))
    implicit = apply_satellite_membership_constraints([original], CONSTRAINT)
    explicit = apply_satellite_membership_constraints(
        [original], CONSTRAINT, policy='high_confidence_pair'
    )
    assert implicit == explicit
    assert implicit[0] == [original]
    assert 'path_supported_constraints' not in implicit[1]
    path, audit = apply_satellite_membership_constraints(
        [original], CONSTRAINT, policy='high_confidence_path'
    )
    assert path == [original]
    assert audit['path_supported_constraints'] == 0


def test_invalid_policy_fails_without_modifying_inputs():
    original = family()
    with pytest.raises(PhylogenyConfigurationError, match='Unsupported'):
        apply_satellite_membership_constraints([original], CONSTRAINT, policy='relaxed')
    assert original == family()


def test_membership_policy_does_not_invalidate_raw_tree_checkpoint():
    config = PhylogenyConfig()
    sequences = {'a': 'MAAA', 'b': 'MBBB', 'c': 'MCCC'}
    changed = replace(config, satellite_membership_rule='high_confidence_path')
    assert _family_tree_input_hash('F0', ('a', 'b', 'c'), sequences, config) == (
        _family_tree_input_hash('F0', ('a', 'b', 'c'), sequences, changed)
    )


def test_replay_parser_defaults_to_direct_support_and_accepts_path():
    argv = ['--fasta-directory', 'in', '--candidate-clusters', 'groups',
            '--output-directory', 'out', '--json', 'report']
    assert build_parser().parse_args(argv).satellite_membership_rule == 'high_confidence_pair'
    assert build_parser().parse_args(argv + [
        '--satellite-membership-rule', 'high_confidence_path'
    ]).satellite_membership_rule == 'high_confidence_path'


def test_replay_rejects_path_policy_without_constraints_before_writing(tmp_path):
    with pytest.raises(SystemExit, match='requires --membership-constraints'):
        main(['--fasta-directory', str(tmp_path),
              '--candidate-clusters', str(tmp_path / 'groups.txt'),
              '--output-directory', str(tmp_path / 'out'),
              '--json', str(tmp_path / 'report.json'),
              '--satellite-membership-rule', 'high_confidence_path'])
    assert not (tmp_path / 'out').exists()


@pytest.mark.parametrize('extra', [[], ['--phylogeny', 'reconcile', '--species_tree_mode', 'infer']])
def test_production_rejects_path_without_satellite_v2(tmp_path, extra):
    from orthohmm.args_processing import process_args
    from orthohmm.parser import create_parser
    args = create_parser().parse_args([
        str(tmp_path), '--accuracy_profile', 'high_sensitivity',
        '--satellite_membership_rule', 'high_confidence_path', *extra])
    with pytest.raises(SystemExit) as error:
        process_args(args)
    assert error.value.code == 1


def test_production_passes_opt_in_rule_without_changing_default(tmp_path, monkeypatch):
    import orthohmm.args_processing as processing
    from orthohmm.parser import create_parser
    monkeypatch.setattr(processing, 'validate_phylogeny_config', lambda config, files: config)
    default = processing.process_args(create_parser().parse_args([str(tmp_path)]))
    assert default['satellite_membership_rule'] == 'high_confidence_pair'
    args = create_parser().parse_args([
        str(tmp_path), '--accuracy_profile', 'high_sensitivity',
        '--phylogeny', 'reconcile', '--species_tree_mode', 'infer',
        '--phylogeny_candidates', 'satellite_v2',
        '--satellite_membership_rule', 'high_confidence_path'])
    configured = processing.process_args(args)
    assert configured['satellite_membership_rule'] == 'high_confidence_path'
