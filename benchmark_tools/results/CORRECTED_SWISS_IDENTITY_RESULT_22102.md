# Corrected SwissTrees Identity Analysis

The [frozen protocol](CORRECTED_SWISS_IDENTITY_PROTOCOL_20260923.md) and
preparer were committed at 21cf7bb before alignment execution. Job 22102
completed 0:0 in 29 seconds with all 18 families and 563 reference proteins.
No alignment failed. Its times are shared-node preparation costs, not
controlled inference timings. No DGX access was required.

The [independent feature check](corrected_swiss_identity_admission_22102.json)
re-reads corrected FASTAs and every alignment, verifies exact accession and
ungapped-residue identities, commands, hashes, per-family status, native
manifests and successful accounting. A scalar implementation recomputes
10,765 unordered pair identities independently of the preparer's NumPy
implementation; all family means agree within 1e-12. It shares the Biopython
FASTA parser and does not establish an independent alignment method.
The frozen and scalar median family identities both equal
0.42912207396202784. Nine families fall in each bin; none is missing.

The [32-row descriptive table](swiss_identity_strata_20260923/scores.md)
retains eight method slots across all/lower/higher/missing bins. Seven methods
have admitted corrected counts; OrthoMCL remains missing. Each row recomputes
macro P/R and harmonic F1 from selected family counts with the existing
raw/2+1 convention, not pooled pairs or average family F1. Full-family
statistics reproduce the admitted overall point estimates. The manifest
retains bin membership, native output semantics and unrounded differences.

| Configuration | Lower Identity F1 (%) | Higher Identity F1 (%) |
|---|---:|---:|
| OrthoHMM sensitive | 64.023 | 72.435 |
| OrthoHMM phylogenetic | 73.973 | 91.914 |
| OrthoFinder full | 76.705 | 92.600 |
| OrthoFinder sequence-only diagnostic | 65.457 | 72.561 |
| SonicParanoid | 73.372 | 86.055 |
| Proteinortho | 60.488 | 81.853 |
| FastOMA, supplied tree | 69.618 | 85.377 |
| OrthoMCL | NA | NA |

All seven admitted methods have lower F1 in the lower-identity bin.
Phylogenetic OrthoHMM minus full OrthoFinder is -2.732 percentage points in
the lower bin and -0.685 in the higher bin. These are descriptive patterns,
not significance, equivalence, interaction or causal claims. Identity is
alignment-dependent, not calibrated evolutionary divergence; bins can differ
in domains, family size, taxa and composition. This does not resolve missing
independent fragment or duplication-history annotations.

All 64 focused verifier/preparer/helper/export tests passed. Recheck using
fresh destinations and the retained local raw evidence:

```sh
python -m benchmark_tools.verify_corrected_swiss_identity --root . \
  --output /tmp/swiss-identity-verification-new.json
python -m benchmark_tools.export_swiss_identity_strata \
  --counts benchmark_tools/results/qfo_fastoma_swiss_uncertainty_22098.json \
  --admission benchmark_tools/results/corrected_swiss_identity_admission_22102.json \
  --output /tmp/swiss-identity-scores-new
```

The exporter pins the retained admission bytes
cb04162af62fbd58fcfe8f02cbb78bc53bcf49ac20a487aae89911bc3a4e7b2d;
it does not silently substitute a later verification report. Alignment
preparation remains reproducible via the retained batch, but executing that
batch again would fail rather than overwrite the completed output. Raw
alignments remain outside Git. These commands do not re-run native orthology
inference or make this a complete portable publication archive.
