# Regression Audit At 970b496

## Unit Suite

`python -m pytest -q tests/unit`:5259passed,9skipped in92.13seconds.
This is unit-suite evidence, not an all-tests success claim.

The initial unrestricted`python -m pytest -q` traversed retained worktrees
because no test-discovery configuration is present. It was interrupted
before any tests ran (51.37seconds). A subsequent`pytest -q tests` included
legacy integration tests and failed at MCL execution; it was interrupted
after five failures. The default shell MCL is the OrthoMCL-specific02-063
build. Main-worktree sample outputs were regenerated during this invocation;
those already-dirty generated files were neither staged nor reverted.
No scientific Slurm job was stopped or restarted.

## Isolated Integration Evidence

Created a detached worktree at970b496:
`benchmarks/work/publication_regression_970b496`.
Selected MCL14-137 from OrthoFinder's bundled bin directory only for the
test command. One initial restricted-PATH invocation completed5failed/3passed
in2.11seconds but did not preserve the usual dependency PATH, so it is not
the retained dependency-complete run. The subsequent invocation prepended
that MCL directory to the original PATH:

```sh
env PATH="/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/orthofinder_3.1.5/lib/python3.12/site-packages/orthofinder/bin:$PATH" \
  python -m pytest -q tests/integration --tb=short --junitxml=ABSOLUTE_XML_PATH
```

It completed5failed/3passed in65.75seconds. Retained[JUnit evidence](integration_regression_970b496_20260918.xml),
SHA-256`708297a059fd917c76c0e50845052a78979d222c397e0161dbfe465e8f9cc91f`.
The first failures compare gene-count text with different species-column
and group ordering. Two other assertions compare individual FASTA files.
An independent normalized partition comparison found exactly matching
gene-membership sets for both simple(38genes) and long-name(1155genes)
fixtures. That observation does not validate every failed output assertion.

Required follow-up: audit gene-count columns by species and group membership,
audit FASTA identities/sequences independent of presentation order, and
isolate integration outputs in temporary directories. Preserve meaningful
regression checks rather than merely updating expected files. Do not claim
the complete test suite passes until these failures are resolved.
