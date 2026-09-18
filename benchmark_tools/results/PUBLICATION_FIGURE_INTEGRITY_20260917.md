# Retained Figure Integrity Audit

`audit_publication_figures.py` checked 14 explicitly retained figure manifests:
accuracy/QfO endpoints, method diagram, OrthoBench factorial, YGOB validation,
sequence-search controls, OrthoBench feature strata, parameter neighborhood,
species-tree perturbations, corrected variable/fixed-length simulations,
simulation tree robustness, WGD application, SwissTrees comparator intervals
and SwissTrees domain strata.

All49 output records and all83 file-reference occurrences match their saved
byte counts and SHA-256 identities. This includes PNG/PDF/SVG files, the WGD
example table, directly recorded numerical inputs, source scripts and other
manifest evidence. Manifests and inspected files were rechecked after reading.
The machine-readable report is `publication_figure_integrity_20260917.json`.

The explicit inventory excludes superseded and defective-runtime figures.
The retained corrected fixed-length simulation panel is still a failure
diagnostic, not a successful OrthoFinder accuracy comparison. QfO factorial
and controlled scaling figures cannot be added until their results are
complete and admitted.

## Recoverable Detached Dependency

One referenced file is not tracked at its local location in the main checkout:
`benchmarks/work/publication_method_native_v2/benchmark_tools/replay_high_sensitivity.py`.
It is evidence for the method diagram. The file matches the recorded SHA-256
`852c3e4fc1a53de7e6046aa78324da587376c0df6a0db1cd8265af65c75bea0f`.

The detached worktree HEAD was checked as
`7f3a9e40dd7e79f842cc2c11fb8b548f9a802806`. Reading
`git show 7f3a9e4:benchmark_tools/replay_high_sensitivity.py` produces the same
SHA-256. Thus these source bytes can be recovered from the retained Git
revision; the final portable bundle must export them explicitly rather than
requiring this machine's worktree path. This does not verify every compiled
runtime dependency or make the full method portable.

## Scope and Reproduction

Nine tests check matching files, tracked-versus-local classification, changed
or missing artifacts, duplicate outputs, absent formats, output paths outside
their figure directory, conflicting hashes and unsupported record schemas.
They pass. Actual audit completed successfully without modifying any figure
or historical manifest.

```bash
python -m pytest -q tests/unit/test_audit_publication_figures.py
python benchmark_tools/audit_publication_figures.py --repo . --output /tmp/new-figure-integrity-report.json
```

Use a fresh output path. This audit checks retained byte identity, not plotted
numerical correctness, rendering, statistical validity, raw-data lineage or
copyright permissions. Existing per-figure scientific and visual checks remain
separate. Git tracking is not proof that local bytes equal a public commit;
transitive raw-data dependencies are not traversed. No complete archive or
publication-readiness claim follows from these checks.
