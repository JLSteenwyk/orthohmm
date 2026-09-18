# Interval Native Recipe Packaging Correction

The original integration attempt21807 failed at worker import before native
inference: the recipe omitted`benchmark_tools/__init__.py`. After the launcher
changed to the frozen native working directory, a regular package in the old
checkout took precedence over the recipe's namespace package. The worker
could not import`benchmark_tools.measure_native_counter_step` there.

Preserve the original recipe, manifest, failed scheduler records, wrapper
verification and step logs. This is a packaging failure, not a biological
method failure or an adverse timing result. No completed scientific run is
discarded. Do not edit or resume the failed output directories.

Create a fresh`interval_native_recipe_v2` including the existing tracked
package initializer, and relocate outputs/cache prefixes to
`interval_native_smoke_v2`. All scientific commands, environments, input
identities, measurement code, thresholds and resources remain unchanged.
Add a regression test executing the copied worker import from a checkout
containing a competing regular package, verifying that imports resolve to
the copied recipe. Verify every transferred file before a new sequential
three-task array. This addendum replaces only the original protocol's output
prefix and recipe packaging; all scientific limits and retention rules stand.
