import pytest
import matplotlib.pyplot as plt

from benchmark_tools import plot_species_tree_robustness as plotter


def report():
    native = {"status": "native_group_output_verified", "native_outputs_validated": True, "accuracy_evaluated": False}
    return {"baseline": plotter.BASELINE, "replicates": 20000, "seed": 20260918, "alpha": .05,
            "families": list(range(70)), "multiplicity": "Bonferroni tail adjustment over 18 reported contrasts/metrics",
            "native_validation": {"status": "native_panel_validated_unscored", "accuracy_evaluated": False,
                "control_validation": {"status": "equivalent", "native_validation": native},
                "variants": {label: {"native_validation": native} for label in plotter.PERTURBATIONS}},
            "point_estimates_percent": {label: {metric: 70. for metric in plotter.METRICS}
                for label in (plotter.BASELINE, *plotter.PERTURBATIONS)},
            "comparisons": {label: {"versus": plotter.BASELINE, "family_f1_wins": 0, "family_f1_ties": 70,
                "family_f1_losses": 0, "metrics": {metric: {"difference_percentage_points": 0.,
                    "paired_percentile_ci": [-1., 1.], "bonferroni_percentile_ci": [-2., 2.]}
                    for metric in plotter.METRICS}} for label in plotter.PERTURBATIONS}}


@pytest.mark.parametrize("change", ["missing", "baseline", "interval", "point", "nonfinite", "multiplicity"])
def test_invalid_figure_evidence_rejected(change):
    data = report()
    row = data["comparisons"][plotter.PERTURBATIONS[0]]
    if change == "missing":
        data["comparisons"].pop(plotter.PERTURBATIONS[0])
    elif change == "baseline":
        row["versus"] = "other"
    elif change == "interval":
        row["metrics"]["f_score"]["bonferroni_percentile_ci"] = [0., 0.]
    elif change == "point":
        row["metrics"]["f_score"]["difference_percentage_points"] = 1.
    elif change == "nonfinite":
        data["point_estimates_percent"][plotter.BASELINE]["f_score"] = float("nan")
    elif change == "multiplicity":
        data["multiplicity"] = "3 endpoints"
    with pytest.raises(ValueError):
        plotter.validate(data)


def test_all_scores_and_eighteen_effects_rendered():
    fig = plotter.plot(report())
    try:
        fig.canvas.draw()
        assert len(fig.axes) == 4
        assert len(fig.axes[0].texts) == 21
        assert all(len(ax.lines) == 20 for ax in fig.axes[1:])
        assert "18 endpoints" in " ".join(text.get_text() for text in fig.texts)
    finally:
        plt.close(fig)
