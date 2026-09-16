from benchmark_tools.render_simulation_results import render
from benchmark_tools.summarize_simulation_panel import CONDITIONS, METHODS, SEEDS, summarize


def test_all_failures_render_as_na_with_counts_and_no_paired_estimates():
    rows = [{"condition": c, "method": m, "seed": s, "status": "failed", "reason": "fixture"}
            for c in CONDITIONS for m in METHODS for s in SEEDS]
    text = render(summarize(rows))
    assert "| baseline | OrthoFinder full | 0 | 10 | 0 | NA | NA | NA |" in text
    assert "| baseline | OrthoHMM high sensitivity | 0 | NA | NA | NA |" in text
    assert "conditional on success" in text and "not pooled gene pairs" in text
