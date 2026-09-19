from copy import deepcopy

import pytest

from benchmark_tools import reproduce_corrected_swiss_strata as module
from benchmark_tools.bootstrap_corrected_swiss_strata import bootstrap
from benchmark_tools.run_corrected_swiss_strata import assemble
from tests.unit.test_run_corrected_swiss_strata import fixture


@pytest.fixture(scope="module")
def result():
    factorial, comparator, strata = fixture()
    for m, source in enumerate((factorial["cells"][4], factorial["cells"][7], comparator)):
        for i, row in enumerate(source["families"]):
            tp, fp = (i + m) % 4, (i + 2 * m) % 4
            row["counts_without_prior"] = dict(TP=tp, FN=3 - tp, FP=fp, TN=3 - fp)
    counts, membership = assemble(factorial, comparator, strata)
    value = bootstrap(counts, membership)
    value.update(status="corrected_swiss_primary_stratified_intervals", scientific_inputs_admitted=True,
        uncertainty_admitted=True, reconstructed_counts=dict(factorial=factorial, orthofinder=comparator))
    return value


def test_independent_family_sums_reproduce_all_primary_endpoints(result):
    assert module.verify(result) == 27


@pytest.mark.parametrize("problem", ["seed", "point", "interval", "wins", "interaction", "cell", "method", "missing"])
def test_modified_result_rejected(result, problem):
    value = deepcopy(result)
    if problem == "seed":
        value["seed"] += 1
    elif problem == "point":
        value["bins"]["lower"]["comparisons"][0]["metrics"]["F1"]["difference"] += .01
    elif problem == "interval":
        value["bins"]["higher"]["comparisons"][0]["metrics"]["F1"]["bonferroni_percentile_ci"][0] += .01
    elif problem == "wins":
        value["bins"]["lower"]["comparisons"][0]["metrics"]["F1"]["family_wins"] += 1
    elif problem == "interaction":
        value["interactions"][0]["metrics"]["PPV"]["difference"] += .01
    elif problem == "cell":
        value["reconstructed_counts"]["factorial"]["cells"][4]["cell"] = "p1c0r0"
    elif problem == "method":
        value["reconstructed_counts"]["orthofinder"]["method"] = "orthofinder_sequence_only"
    else:
        value["bins"]["missing"]["comparisons"][0]["metrics"]["F1"]["difference"] = 0
    with pytest.raises(ValueError):
        module.verify(value)
