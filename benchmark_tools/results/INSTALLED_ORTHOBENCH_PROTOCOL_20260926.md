# Full Installed-Package OrthoBench Reproduction

## Scope and Justification

The setup-overlay wheel and relocated/rebuilt external tools passed small
fixtures, but those do not reproduce a publication dataset. Run one fresh,
end-to-end OrthoBench inference with the separately installed frozen-source
wheel, rebuilt MAFFT core and identified FastTree 2.2.0 executable. This is a
new artifact/environment reproduction check, not independent biological
validation, a replacement for historical predictions, or controlled timing.
Existing expensive checkpoints remain untouched. No parameter selection is
performed from the new result.

## Fixed Inputs and Configuration

Use the 12 byte-checked proteomes and 251,378 unique identifiers in the pinned
factorial preparation manifest. Copy only FASTAs into a fresh private directory.
Do not copy search, graph, candidate, tree or reconciliation checkpoints.
The plan binds the installed-source audit, all recorded installed inputs,
MAFFT's 34 helpers/launcher, FastTree, the historical comparator partition and
reference files. The inference subprocess receives no reference-label path.

Invoke the installed Python with `-I -m orthohmm`, builtin search, Leiden,
`high_sensitivity`, `satellite_v2` candidates, inferred species tree with
`min_variance` rooting, `species_overlap` gene-tree root rule and
`positive_paralogy` pair rule. Use 32 CPUs and one thread per search worker;
BLAS/OpenMP numerical thread limits are one. `MAFFT_BINARIES` points explicitly
to the private rebuilt helper directory. Strip PYTHONPATH/PYTHONHOME from the
child environment; use `/usr/bin:/bin` PATH and absolute external-tool paths.

The CPU-baseline wheel and dependency versions differ from the historical
runtime. The test must record discrepancies, not assume full equivalence
from frozen scientific source. No automatic retry, resume, new default or
benchmark-score substitution is allowed. An existing inference directory
causes refusal. Errors preserve logs and a failed execution report.

## Execution and Validation Gates

Submit one Slurm job on bizon, partition gpu, 32 CPUs, 128 GiB RAM and 24 hours.
The child has a 23-hour timeout. The wrapper pins the executor Git commit and
plan SHA-256; the executor verifies input/tool/package records before and after.
The installed module path must be under its intended venv at preparation.
The unchanged historical sources and scores are retained for comparison.

`/usr/bin/time -v` records descriptive process/child resource observations.
Its maximum-RSS field is not a summed concurrent-process peak or cgroup memory
measurement. No timing number from this shared host enters the 27-run controlled
scaling panel. No DGX access or unrelated process/service change is authorized.

Successful native exit is only `native_completed_pending_independent_scientific_readback`.
Subsequent readback must check scheduler completion, coverage exactly once,
fresh inferred trees/reconciliation, reference-file identities, aggregate and
per-family OrthoBench scores, and label-invariant partition agreement with
the historical P1/C1/R1 partition. Preserve disagreement and investigate its
stage before any claim of full-dataset reproduction. Native success alone
does not admit scientific scores or close the publication goal.

```sh
/usr/bin/python3 -S -m benchmark_tools.run_installed_orthobench --repo . \
  --prepare benchmarks/work/publication_installed_orthobench_20260926
python -m pytest -q tests/unit/test_run_installed_orthobench.py
# Submit the tracked .slurm wrapper with a clean pinned executor worktree,
# its commit, the absolute plan path and the plan SHA-256.
```

Preparation is complete; job submission and live state are recorded separately
in the progress ledger. This document does not assert completion of inference.
