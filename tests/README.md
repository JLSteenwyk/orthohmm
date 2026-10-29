# Application Tests

This suite tests OrthoHMM, its CLI, installation, and native search kernels.
The benchmark, manuscript, and historical replay suites are preserved on
[research-archive-20261029](https://github.com/JLSteenwyk/orthohmm/tree/research-archive-20261029),
at pre-split commit `ebfa635ca786796287c9ff2a3caba1c2618f016e`.
No research checkout, archived dataset, or old Git revision is needed here.

## Run Tests

From a fresh checkout in an activated Python 3.10+ environment:

```sh
python -m pip install . -r tests/requirements.txt
make kernels
make test
make test.coverage
```

The checkout-kernel helper requires a fresh output directory and refuses
inherited shared libraries. For another fresh build, use a fresh worktree or
remove only your generated `orthohmm/search/csrc/*.so` and kernel build directory.
Normal installed-package inference does not require this test-only build step.

The default suite collects every application test, including installation and
CLI checks. No private benchmark markers or research-specific exclusions are
used. Tests of optional HMMER/MCL or MAFFT/FastTree workflows report capability
skips if those executables are unavailable. Linux CI installs these tools;
macOS CI exercises the portable core and legacy HMMER/MCL workflow.

## Coverage

Coverage measures `orthohmm/`, including Python child processes and workers
that end with `os._exit`. The application gate is 85% line coverage. Unit and
integration data are retained separately and combined into `coverage.xml`,
`coverage.json`, and `.coverage`. Tests verify that combination preserves both
inputs, missing application source is an error, and child working-directory
changes do not include test-harness code in reports.

Python coverage does not measure C/CUDA instructions. Independent numerical
oracles compare Python reference bodies and compiled search kernels. GPU and
native capability skips are not counted as successful native execution.

Coverage now describes the application repository only. Removing research
workflows changes the repository-wide denominator; any resulting percentage
change is not itself additional test coverage.
