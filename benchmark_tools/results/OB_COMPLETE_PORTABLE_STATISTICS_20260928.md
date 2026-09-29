# Portable OrthoBench Statistical Replay

The complete eight-method comparison now has a standalone statistical replay.
It requires just `check_ob_complete_uncertainty.py`, the
[compact input](ob_complete_portable_statistics_20260928.json), Python and NumPy.
No project import, repository checkout, absolute historical path, network,
reference membership file or prediction file is needed in portable mode.

The input contains 560 rows: 70 named RefOGs for each of eight methods,
with gene count and TP/FP/FN sufficient statistics. It includes the frozen
analysis settings and expected results, not sequences or gene memberships.
Export checked the retained local input records using the existing loader;
source result names, sizes and hashes are included without local paths.
The input is 114,026 bytes, SHA-256
`cdb30ff987100843f3127a39a75eb3649123a136229b8b41e6ed156209aff30e`.

Copying the two files into `/tmp/orthohmm-ob-complete-statistics-20260928`
and running the copied script with Python `-I -B`, cwd `/tmp`, succeeded.
The [replay receipt](ob_complete_portable_replay_20260928.json) binds both
copied files and records Python 3.10.13 and NumPy 2.2.6. All 21 contrasts,
nominal/adjusted intervals and seven family-outcome triples agree; maximum
error is 2.14e-14 percentage points. The original local-input mode also
passes after refactoring. Fifty-four focused tests pass, including an actual
isolated subprocess outside the checkout and checksum/schema rejection.

With the checker and count file in the current directory, use a fresh output:

```sh
python -I -B check_ob_complete_uncertainty.py \
  --portable ob_complete_portable_statistics_20260928.json \
  --sha256 cdb30ff987100843f3127a39a75eb3649123a136229b8b41e6ed156209aff30e \
  --output replay
```

The output explicitly labels its scope as derived-count replay, not raw-input
verification. It records the actual Python/NumPy versions; these dependencies
are not bundled or installed by the command. This is same-host relocation,
not a new-machine installation, scoring reproduction, inference reproduction,
complete study archive or deposition. Shared NumPy RNG/quantile code and
conditional family-exchangeability assumptions remain. The full-reference
weighted statistic, 100,000 draws, seed, 21-endpoint correction, thresholds
and all retained scientific conclusions are unchanged.
