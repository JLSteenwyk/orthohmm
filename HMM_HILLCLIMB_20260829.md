# HMM-Centered Accuracy Hill-Climb, 2026-08-29

## Outcome

The primary outperformance target was not reached. The strongest generalizable
candidate is selective weakest-member jackknife calibration for profile HMMs
whose seed cluster contains at most one gene per represented species. It
improved cached-search OrthoBench replay from F=70.4 to F=70.5, but remains 2.2
points below OrthoFinder 3.1.5 at F=72.7. It is retained as an experimental
opt-in and is not the default.

No candidate passed the predeclared OrthoBench gate, so no candidate-specific
QfO or Three Kingdoms run was made. This prevents final-benchmark tuning and
avoids attributing the prior high-sensitivity QfO and Three Kingdoms scores to
the new candidate.

## Cross-Benchmark Comparison

`-` means that exact tool/version/method was not assessed on that benchmark.
Historical OrthoBench tool scores use the recorded exploration outputs; the
primary OrthoFinder row is the fresh 3.1.5 result.

| method | OrthoBench F | P | R | exact | QfO mean | Three Kingdoms F | P | R |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OrthoFinder 3.1.5 DIAMOND | 72.7 | 66.1 | 80.9 | 19 | 0.782071 | - | - | - |
| OrthoHMM selective HMM jackknife | 70.5 | 78.7 | 63.8 | 12 | - | - | - | - |
| OrthoHMM high sensitivity | 70.4 | 78.9 | 63.5 | 13 | 0.682548 | 0.8263 | 1.0000 | 0.7040 |
| FastOMA, historical | 50.8 | 52.0 | 49.7 | 21 | - | - | - | - |
| SonicParanoid, historical | 46.8 | 36.5 | 64.9 | 19 | - | 0.8859 | 1.0000 | 0.7952 |
| ProteinOrtho, historical | 45.1 | 97.1 | 29.3 | 6 | - | 0.9305 | 1.0000 | 0.8700 |
| OrthoFinder 2.5.5 DIAMOND | - | - | - | - | - | 0.8552 | 1.0000 | 0.7470 |
| OrthoHMM prior CPM-M4 experiment | - | - | - | - | - | 0.9070 | 1.0000 | 0.8298 |

Against OrthoFinder 3.1.5 on OrthoBench, selective jackknife is +12.6 points
in precision, -17.1 points in recall, -2.2 points in F-score, and -7 exact
RefOGs. The limiting problem is recall, not precision.

## QfO Metrics

These are fresh official assessments of the existing high-sensitivity method,
not the selective-jackknife candidate. The candidate was gated out before QfO.

| metric | OrthoHMM high sensitivity | OrthoFinder 3.1.5 | delta |
| --- | ---: | ---: | ---: |
| VGNC F | 0.667448 | 0.988229 | -0.320781 |
| SwissTrees F | 0.673851 | 0.859160 | -0.185310 |
| TreeFam-A F | 0.575870 | 0.742887 | -0.167017 |
| EC | 0.931045 | 0.941982 | -0.010937 |
| GO | 0.472282 | 0.468166 | +0.004117 |
| FAS | 0.774794 | 0.692001 | +0.082794 |
| mean | 0.682548 | 0.782071 | -0.099522 |

OrthoHMM wins GO and FAS but trails substantially on all three reference-tree
F-scores. It produced 390,817 groups versus OrthoFinder's 46,593, consistent
with severe fragmentation rather than a globally permissive false-positive
problem.

## Experiments

| HMM-centered mechanism | OrthoBench F | decision |
| --- | ---: | --- |
| baseline high sensitivity | 70.4 | comparator |
| selective single-copy weakest-member calibration | 70.5 | retain opt-in |
| HMM-supported component refinement | 70.5 | neutral; reverted |
| selective absent-species second pass | 70.4 | rejected and reverted |
| direct profile-only fallback | 70.3 | rejected and reverted |
| reciprocal split-cluster repair | 70.3 | rejected and reverted |
| pair-seeded profiles | 69.9 | rejected and reverted |
| profile-profile repair, similarity 0.40 | 69.9 | rejected and reverted |
| per-target-residue profile calibration | 69.7 | rejected and reverted |
| global weakest-member calibration | 69.1 | rejected |
| profile-profile repair, similarity 0.60 | 70.5 | inactive and reverted |
| global second profile pass | 70.2 | previously rejected |
| species-balanced profile emissions | gate failed | rejected before OrthoBench |

Minimum two- and three-species profile gates scored F=69.9 and F=69.8. Global
BLOSUM45, a wider uncalibrated BLOSUM62 band, reduced-alphabet fallback,
unconditional iteration, and scalar RBNH changes had already failed earlier
screens and were not repeated.

The independent synthetic holdouts were useful for rejecting unsafe methods,
but several methods with perfect negative rejection still regressed complete
proteomes. Their generated families did not adequately model multidomain
proteins, fragments, and recent paralogs. Final OrthoBench values were not used
to select another threshold after these failures.

## Resources

| run | CPUs | wall time | peak RSS | note |
| --- | ---: | ---: | ---: | --- |
| OrthoHMM high sensitivity, OrthoBench | 32 | 2,185.48 s | 11.90 GiB | fresh production; co-scheduled |
| OrthoFinder 3.1.5, OrthoBench | 32 search / 8 analysis | 3,421.07 s | 4.57 GiB | fresh production |
| OrthoHMM selective jackknife, OrthoBench | 32 | 347.15 s | 5.56 GiB | cached-search replay only |
| OrthoHMM high sensitivity, QfO | 32 | 69,438 s | 16.70 GiB | fresh production |
| OrthoFinder 3.1.5, QfO | 32 | 80,953 s | 6.10 GiB | fresh production |
| OrthoHMM high sensitivity, Three Kingdoms | 32 | 6,556.43 s | 11.12 GiB | fresh production |

The replay is not a full-runtime comparator. OrthoHMM's separate QfO run was
14.2% faster than OrthoFinder's but used 2.74x peak memory; this is descriptive,
not a controlled co-scheduled speed claim.

## Reproducibility

Key checksums:

| artifact | SHA-256 |
| --- | --- |
| selective-jackknife final clusters | `55f9e3a2912bf3b3aabb1371a34f8dcae91fed6144cbd7273aea47253fa1263a` |
| baseline replay final clusters | `8ee100f4370fef6e3e6a3e0fd36112adcc3aefeb60197331115667c519ca121f` |
| cached normalized hits | `78a5af40ea2683a69e1baefefb0966c3549cffe71b24bda03b9ba4a6e29e65e1` |
| OrthoFinder 3.1.5 OrthoBench groups | `e0bf8d21b7ddc4494f0df1224dc2cd57d964cb92507eb40bd1f058cbc6f47af7` |
| OrthoFinder 3.1.5 QfO pairs | `5581c9b81314e7d074b04564b4ce7d0ac66d909c0e512e88ba30e14e05dba6a2` |
| OrthoHMM high-sensitivity QfO pairs | `3abde0f64063ce8ce0432d358faa340556df95a8b5d4aa6d30b0e19deea7cb0f` |
| OrthoHMM Three Kingdoms groups | `28ba45202afc7f51e719d59405674023c9b284a2d5ea18575b03c7b949e13db1` |

Relevant commits:

| phase | commit |
| --- | --- |
| reciprocal profile repair implementation/result | `528bdd5`, `6a60b77` |
| accepted selective jackknife implementation/result | `6d48b35`, `90d591a` |
| length calibration implementation/result | `ba1fed4`, `6a597e6` |
| profile-profile implementation/results | `2c9024f`, `a803b4e`, `e8f66eb`, `eb6efa8` |
| component refinement implementation/result | `63fc686`, `368ba57` |
| direct fallback implementation/result | `da1c881`, `b4fd0b5` |
| pair-seeded profiles implementation/result | `c1615a5`, `3b8a0f3` |
| species-completion implementation/fix/result | `6161532`, `54d81d0`, `9ad6b5f` |
| failed-variant production cleanup | `23008d8` |

Verification after cleanup: `274 passed` in the unit suite. Scoped Ruff checks
pass for all retained modified implementation, harness, and test files. The
repository-wide Ruff scan still has 37 pre-existing findings. No scheduler jobs
remain active.

## Blocker and Recommendation

Current Viterbi kernels return scores only. They do not expose alignment start,
end, traceback, or model/query coverage. The strongest remaining hypothesis is
domain-coverage-aware HMM evidence, potentially followed by multiple profiles
for multimodal or paralog-rich families. Implementing and calibrating that
requires both:

1. score-and-coverage output from the CPU HMM kernels without violating the
   current memory bound; and
2. an independently labeled development corpus containing multidomain
   proteins, fragments, domain shuffles, and recent paralogs.

Without those two pieces, further scalar threshold hill-climbing would tune the
final benchmark and violate the no-leakage constraint. Keep `standard` as the
default. Retain high sensitivity and selective single-copy jackknife as opt-in
research modes. Do not run another full QfO assessment until a coverage-aware
candidate beats F=72.7 on an untouched OrthoBench production run.
