"""Conditional native FAS intervals with two unknown finite-population means."""

import math
from numbers import Integral

from benchmark_tools.native_fas_sampling_interval import count_interval, mixture_weight, require


def bounded_mean_interval(size, mean, error):
    require(isinstance(size, Integral) and not isinstance(size, bool) and size >= 0,
            "Nonnegative integer sample size required")
    require(math.isfinite(error) and 0 < error < 1, "Invalid mean error")
    if size == 0:
        require(mean is None, "Empty sample has no mean")
        return [0., 1.]
    require(mean is not None and math.isfinite(mean) and 0 <= mean <= 1,
            "Invalid sample mean")
    width = math.sqrt((math.log(2) - math.log(error)) / (2 * size))
    return [max(0., mean - width), min(1., mean + width)]


def method_interval(population, draws, precomputed_population, precomputed_draws,
                    precomputed_sample_mean, returned, return_sample_mean, component_error):
    """Project three covered components; per-method failure <=3*error."""
    require(all(isinstance(x, Integral) and not isinstance(x, bool)
                for x in (precomputed_population, precomputed_draws))
            and 1 <= precomputed_draws <= precomputed_population, "Invalid precomputed counts")
    counts = count_interval(population, draws, returned, component_error)
    old = bounded_mean_interval(precomputed_draws, precomputed_sample_mean, component_error)
    if precomputed_draws == precomputed_population:
        old = [precomputed_sample_mean, precomputed_sample_mean]
    new = bounded_mean_interval(returned, return_sample_mean, component_error)
    weights = [mixture_weight(population, g, draws, precomputed_draws) for g in counts]
    corners = [w * a + (1 - w) * b for w in weights for a in old for b in new]
    return dict(success_count_bounds=list(counts), precomputed_mean_bounds=old,
                return_mean_bounds=new, expected_native_mean_bounds=[min(corners), max(corners)],
                component_error=component_error, method_error_bound=3 * component_error,
                target="expected_native_post_attrition_ratio_under_fixed_design")
