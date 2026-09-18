# Relocated Domain-Stratified Reproduction

Extended the existing SwissTrees reproduction runner with optional
`--include-domain-strata`. Its original comparator analysis remains enabled.
The additional workflow exports the annotation inventory, frozen protocol,
analysis, expected results/table, plotter and import dependencies.

## Executed Check

Exported 22 files from committed analysis revision
`87a9a843035d987343915b98a3490a50a446c519` into
`/tmp/orthohmm-swiss-domain-reproduction-20260917/source`, outside the working
checkout. Used the existing patched isolated analysis environment: Python
3.10.13, NumPy2.2.6, Matplotlib3.10.8 and Pillow12.3.0, with all11 packages
recorded. This was not a newly installed environment or another machine.
The runner itself is recorded separately by its exact file hash in the report.

All four subprocesses succeeded in isolated Python mode, with user-site imports
disabled, Python/loader overrides removed, numerical threads set to one, and
a fresh plotting configuration directory:

1. Comparator bootstrap: scientific JSON exact, Markdown byte-identical.
2. Comparator figure: PNG/PDF/SVG and manifest generated.
3. Domain-stratified bootstrap: scientific JSON exact, Markdown byte-identical.
4. Domain-stratified figure: PNG/PDF/SVG and manifest generated.

Only relocated provenance paths may differ. Input, source and ordered helper
hashes/byte counts must remain identical. The comparison includes every
interval, family/bin membership, point estimate, random seed, numerical-library
version, secondary descriptive result and limitation. Plotting occurs only
after the complete corresponding statistical content reproduces.
All exported bytes are rechecked after execution. The relocated domain PNG
was visually inspected for legibility and overlap; no bitwise figure identity
is claimed because renderer metadata can vary.

Eight new comparison tests reject altered domain intervals, bin membership,
versions, input order, helper/source bytes and extra result fields while
allowing path relocation. Combined reproduction/domain-analysis/plot tests:
27passed. Actual relocated execution supplements those mutation tests.

Evidence: `swiss_domain_relocated_reproduction_20260917.json`, SHA-256
`22a9b42691870d0965b9899b2c614f94d05cc9ea51d5594ef618ef8be8128ba7`.

## Reproduce

Use the current patched hash-pinned analysis lock and a new export/report path:

```bash
python benchmark_tools/reproduce_swiss_comparators.py --revision 87a9a84 --python benchmarks/work/swiss_analysis_env_20260917/bin/python --include-domain-strata --output /tmp/swiss-domain-new-export --report /tmp/swiss-domain-new-report.json
```

For a new installation, follow the environment setup in
[the original reproduction note](SWISS_RELOCATED_REPRODUCTION_20260917.md)
using the current lock, not the historical Pillow12.2.0 pin.

## Remaining Scope

This reproduces statistics from retained count and annotation inventories.
It does not regenerate Pfam annotations, rerun inference or QfO scoring,
establish biological independence, validate another platform, settle third-party
redistribution rights, or complete the publication archive and versioned release.
Historical absolute paths inside source reports remain provenance only; the
analysis processes use exported data, not the original raw installations.
