# Simulation Arithmetic Replay Result

Status: bounded retained-count replay passed, not publication-ready. Base source
`afdcafb6e35278a4224b7e5e08ab8d943779a364`; no timing reservation, workload poll,
DGX access or unrelated process/service action. Native inference, raw scoring
and scientific defaults/results are unchanged.

See the [standalone command and limits](../SIMULATION_ARITHMETIC_REPLAY.md) and
[machine-readable review](simulation_arithmetic_replay_review_20261002.json).
The review is 11,273 bytes, SHA256
`1ec2e5394ea6a63a9c12cce2e6439c84f41af79eb442160bf3cec0a70f63febf`.
Both original reports and the producer's summary source are byte-identical to
the base commit. The independent verifier is 14,573 bytes, SHA256
`f475cb322346a23ae2afeb19a898a9a79fd786645e5db24ce7b3d07221b014b6`.

## Executed Scope

| Panel | Planned outcomes | Complete | Failed | Estimated contrasts | Unavailable contrasts |
| --- | ---: | ---: | ---: | ---: | ---: |
| Fixed length | 280 | 134 | 146 | 0 | 14 |
| Variable length | 280 | 267 | 13 | 14 | 0 |

All method outcome lists, failure fractions and 168 mean cells reproduce,
including unavailable values. Variable-length comparisons reproduce 42 paired
metric effects and 112 scalar bounds within 1e-12 using the original 20,000
draws per contrast. Independent rational metrics and explicit seed summation
replace the producer's floating ratios/matrix product, not its RNG/quantiles.
No unavailable fixed-length OrthoFinder comparison is turned into a tool win.
Failures remain explicit; no pooling across gene pairs or length panels.

The first direct local replay used the earlier 14,382-byte verifier, SHA256
`782de5186b9e09404fb4002c81d01b3b0be658d42c87415b97118fe8f60bf117`.
Its retained receipt is not relabeled as execution of the final strict-shape
source. Before final replay, numeric checks were tightened to reject boolean,
nonfinite and broadcast-shaped observations.

Actual final reproduction runs from three copied files in a fresh `/tmp`
directory, using private Python 3.12.3/NumPy 2.2.6 with `-I -B`. All payloads
are byte-verified before/after. Checkout, `/proc` and `/sys` canaries are
rejected; zero forbidden opens follow. No project modules are imported.
The temporary copy is removed; logs and result receipts remain retained.

Two earlier harness attempts exit one: the guard catches lazy runtime reads
for a standard-library helper, then NumPy random. Both failed receipts/logs
are preserved. The final harness preloads these dependencies before the guard
and preserves the private interpreter path instead of resolving its symlink.
Its third child exits zero at 13:47:59 UTC. The guard is not relaxed to allow
checkout/raw-data access. Warmed runtime and Python-event checks are not OS
containment, hermetic dependencies or cross-host portability. This is a copied
tree test, not a newly extracted archive or public upload.

## Tests and Remote Status

Initial two-module regression: 53 passes in 2.08s. Final panel: **55 passes in
2.53s**, zero failures/errors/skips, comprising 42 new replay cases and 13
existing producer cases. JUnit is 7,491 bytes, SHA256
`b016155eb42f3dc0fffff032b40a97f5b2213414ee31ae3ce95ececc88fdb25a`.
Panels overlap, not additive. Checks include both real reports, corrupted
counts/coverage/truth/design/intervals, explicit failures, diagnostic-parent
hierarchy, unavailable and single-pair cases, toy mean-versus-pooled distinction,
independent paired draws, no-overwrite and a fresh copied standalone CLI.

At 13:51:08 UTC, all eight jobs of prior-source public CI run `37013500430`
are successful at exact `afdcafb6`. No job execution counts are inferred or
downloaded here. This confirms the public profile's status, not full private
raw/native success or execution of this new simulation code. Its post-push
CI outcome must be observed separately.

The replay does not re-admit native outputs/generator/truth history, establish
seed exchangeability or remove success-conditioned selection bias. It does
not validate tree-perturbation/matched-recall simulations or prove superiority
on arbitrary datasets. Controlled timing, other-QfO uncertainty, complete
runtime/rights/release/deposition and final manuscript reconciliation remain
open. The publication goal remains active.
