from benchmark_tools.summarize_three_kingdoms_parity import (
    METHODS,
    format_duration,
    parse_run_metadata,
    parse_score,
)


def test_parse_score(tmp_path):
    score = tmp_path / "score.txt"
    score.write_text(
        "=== tool ===\n"
        "  reference OGs            :        255\n"
        "  reference genes          :      2,035\n"
        "  ref-genes in prediction  :      1,820 (89.4%)\n"
        "  predicted OGs            :     28,019\n"
        "  TP gene pairs            :      5,565\n"
        "  FP gene pairs            :          0\n"
        "  FN gene pairs            :      1,787\n"
        "  precision                :     1.0000\n"
        "  recall                   :     0.7569\n"
        "  F-score                  :     0.8617\n"
    )

    result = parse_score(score)

    assert result["reference_orthogroups"] == 255
    assert result["reference_genes_in_prediction"] == 1820
    assert result["predicted_orthogroups"] == 28019
    assert result["precision"] == 1.0
    assert result["f_score"] == 0.8617
    assert result["reference_gene_coverage"] == 1820 / 2035


def test_parse_run_metadata_and_format_duration(tmp_path):
    metadata = tmp_path / "run_metadata.tsv"
    metadata.write_text("slurm_job_id\t20892\nexit_code\t0\nstate\tcomplete\n")

    assert parse_run_metadata(metadata) == {
        "slurm_job_id": 20892,
        "exit_code": 0,
        "state": "complete",
    }
    assert format_duration(6332) == "1:45:32"
    assert format_duration(None) == "n/a"


def test_fastoma_variant_describes_scored_output():
    fastoma = next(method for method in METHODS if method.key == "fastoma_0_3_5")

    assert fastoma.variant == "final orthologous groups; supplied species tree"


def test_orthomcl_variant_discloses_compatible_conversion():
    orthomcl = next(method for method in METHODS if method.key == "orthomcl_1_4")

    assert "exact-compatible conversion and pair parallelism" in orthomcl.variant
    assert orthomcl.runtime_kind == (
        "measured stage sum; checkpoint-recovered downstream stage"
    )
