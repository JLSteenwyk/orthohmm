import pytest

from benchmark_tools.run_canonical_ob_phylogeny import replay_arguments,run
from benchmark_tools.replay_phylogeny import build_parser


def plan():
    return dict(inputs="/inputs",candidate={"path":"/candidate"},constraints={"path":"/constraints"},
                inference="/output/inference",output="/output",cpu=32,aligner="/mafft",tree_builder="/FastTree")


def test_exact_frozen_reconciliation_settings():
    parsed=build_parser().parse_args(replay_arguments(plan()))
    assert parsed.species_tree_mode=="infer"
    assert parsed.species_tree_rooting=="min_variance"
    assert parsed.root_rule=="species_overlap" and parsed.pair_rule=="positive_paralogy"
    assert parsed.cpu==32 and str(parsed.membership_constraints)=="/constraints"
    assert str(parsed.output_directory)=="/output/inference"
    assert parsed.aligner=="/mafft" and parsed.tree_builder=="/FastTree"


def test_fresh_stage_never_exposes_optional_scoring_or_reuse():
    p=plan()
    p.update(checkpoint_source="/old",official_benchmark="/labels",species_tree="/supplied")
    parsed=build_parser().parse_args(replay_arguments(p))
    assert parsed.checkpoint_source is None and parsed.official_benchmark is None
    assert parsed.species_tree is None and not parsed.unconstrained_membership


@pytest.mark.parametrize("sha",[None,"", "0"*64])
def test_unpinned_or_changed_plan_cannot_run(tmp_path,sha):
    path=tmp_path/"plan.json"
    path.write_text("{}")
    with pytest.raises(ValueError,match="unpinned"):
        run(path,sha)
    assert sorted(p.name for p in tmp_path.iterdir())==["plan.json"]
