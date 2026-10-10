"""Conditional fixed-population intervals for the expected native FAS ratio.

These numerical kernels do not establish their design assumptions or admit
historical benchmark uncertainty. See the frozen native sampling protocol.
"""

import math
from numbers import Integral

from scipy.stats import hypergeom


def require(condition, message):
    if not condition:
        raise ValueError(message)


def validate_counts(population, draws, returned):
    require(all(isinstance(x, Integral) and not isinstance(x, bool)
                for x in (population, draws, returned)), "Counts must be integers")
    require(0 <= returned <= draws <= population and draws <= 9000, "Invalid native count support")


def count_interval(population, draws, returned, error):
    """Invert both inclusive count tails on the integer success-count support."""
    validate_counts(population, draws, returned)
    require(math.isfinite(error) and 0 < error < 1 and error / 2 > 0, "Invalid count error")
    if draws == 0:
        return 0, population
    threshold = error / 2
    support_low, support_high = returned, population - draws + returned

    def tail(successes, upper):
        value = (hypergeom.sf(returned - 1, population, successes, draws) if upper else
                 hypergeom.cdf(returned, population, successes, draws))
        require(math.isfinite(value) and 0 <= value <= 1, "Invalid hypergeometric tail")
        return value

    low, high = support_low, support_high
    while low < high:
        middle = (low + high) // 2
        if tail(middle, True) >= threshold:
            high = middle
        else:
            low = middle + 1
    first = low
    low, high = support_low, support_high
    while low < high:
        middle = (low + high + 1) // 2
        if tail(middle, False) >= threshold:
            low = middle
        else:
            high = middle - 1
    require(first <= low and tail(first, True) >= threshold and tail(low, False) >= threshold,
            "Empty count confidence set")
    return first, low


def mixture_weight(population, successes, draws, precomputed_draws):
    """E[k/(k+R)], retaining the hypergeometric random denominator."""
    validate_counts(population, draws, 0)
    require(isinstance(successes, Integral) and not isinstance(successes, bool) and
            0 <= successes <= population, "Invalid success count")
    require(isinstance(precomputed_draws, Integral) and not isinstance(precomputed_draws, bool) and
            precomputed_draws >= 1, "Positive integer k required")
    if draws == 0 or successes == 0:
        return 1.
    lower, upper = max(0, draws - population + successes), min(draws, successes)
    probabilities = hypergeom.pmf(range(lower, upper + 1), population, successes, draws)
    require(all(math.isfinite(p) and p >= 0 for p in probabilities), "Invalid hypergeometric mass")
    mass = math.fsum(probabilities)
    require(abs(mass - 1) <= 1e-8, "Hypergeometric mass not normalized")
    # Normalize only numerical mass drift; no tail is dropped or approximated.
    return math.fsum(p * precomputed_draws / (precomputed_draws + r)
                     for p, r in zip(probabilities, range(lower, upper + 1))) / mass


def expected_native_mean(population, successes, draws, precomputed_draws, precomputed_mean, return_mean):
    require(all(math.isfinite(x) and 0 <= x <= 1 for x in (precomputed_mean, return_mean)),
            "Means must be in [0,1]")
    weight = mixture_weight(population, successes, draws, precomputed_draws)
    return weight * precomputed_mean + (1 - weight) * return_mean


def method_interval(population, draws, precomputed_draws, precomputed_mean,
                    returned, return_sample_mean, component_error):
    """Two-component rectangle projection; failure probability <=2*error."""
    validate_counts(population, draws, returned)
    count_low, count_high = count_interval(population, draws, returned, component_error)
    if returned == 0:
        require(return_sample_mean is None, "Zero returns have no sample mean")
        mean_low, mean_high = 0., 1.
    else:
        require(return_sample_mean is not None and math.isfinite(return_sample_mean) and
                0 <= return_sample_mean <= 1, "Invalid returned sample mean")
        # G=1 is trivial: its positive-return sample observes the entire mean.
        # The same enclosing interval covers it without invoking an N>=2 theorem.
        width = math.sqrt((math.log(2) - math.log(component_error)) / (2 * returned))
        mean_low, mean_high = max(0., return_sample_mean - width), min(1., return_sample_mean + width)
    weights = [mixture_weight(population, g, draws, precomputed_draws) for g in (count_low, count_high)]
    require(math.isfinite(precomputed_mean) and 0 <= precomputed_mean <= 1, "Invalid precomputed population mean")
    corners = [weight * precomputed_mean + (1 - weight) * mean
               for weight in weights for mean in (mean_low, mean_high)]
    return dict(success_count_bounds=[count_low, count_high], return_mean_bounds=[mean_low, mean_high],
                expected_native_mean_bounds=[min(corners), max(corners)],
                component_error=component_error, method_error_bound=2 * component_error,
                target="expected_native_post_attrition_ratio_under_fixed_design")


def difference_interval(left, right):
    """Project a jointly covered rectangle; no cross-method independence needed."""
    a, b = left["expected_native_mean_bounds"], right["expected_native_mean_bounds"]
    require(len(a) == len(b) == 2 and all(math.isfinite(x) for x in a + b) and
            0 <= a[0] <= a[1] <= 1 and 0 <= b[0] <= b[1] <= 1, "Invalid method bounds")
    return [a[0] - b[1], a[1] - b[0]]
