# Frozen Checkout Runtime Failure

## Finding

OrthoBench replay job **21088** finished native inference with exit zero but
failed its partition-equivalence gate (Slurm FAILED, exit 1, elapsed 00:03:49).
Both non-profile stages match the historical partitions byte-for-byte.
Both profile-enabled stages instead equal their non-profile counterparts.
The native replay records zero profiles built, zero candidates and zero edges.

The inference checkout `benchmarks/work/publication_method_7f3a9e4` contains
native C sources but no compiled `pair_align.so`. The center-star MSA loader
requires that exact library. Calling it directly raises OSError; the profile
worker catches all exceptions and returns None. This makes a missing runtime
dependency look like successful inference without any usable profiles.

[Isolated runtime probe](profile_runtime_missing_20260916.json) records the
source path/hash and actual loader failure. It tests a fixed synthetic cluster,
not reference labels. The failure occurs before profile construction. No
frozen source, binaries or active job inputs were modified.

The same probe [passes in the development checkout](profile_runtime_development_control_20260916.json),
with the identical `msa_profile.py` SHA256 and a recorded native library hash,
producing a 20-position profile. This is a diagnostic positive control, not a
replacement frozen runtime or proof of historical end-to-end equivalence.

## Affected Evidence

- Replay v1 is **not equivalent** and cannot feed the publication factorial.
  Preserve its preflight, logs and verification under
  `benchmarks/results/publication_ob_replay_check_v1`.
- Both simulation method manifests name the same incomplete checkout.
  All 140 fixed-length OrthoHMM metrics records report zero profiles.
  At the live variable-length audit, all 96 available OrthoHMM metrics records
  (48 per mode) also report zero profiles. These counts include failed runs;
  they are not counts of admitted accuracy results.
- Fixed-length OrthoHMM scores are retained as defective-runtime diagnostics,
  **not valid estimates of the intended high-sensitivity or satellite method**.
  Species-tree failures may depend on the missing profile stage; previous
  causal interpretations need reevaluation after a correctly built rerun.
- Variable-length OrthoHMM outcomes must not be used for publication comparison.
  Array 21010 remains unmodified to finish reusable OrthoFinder outputs.
  No variable-length accuracy was inspected to select a correction.
- OrthoFinder's independently demonstrated constant-length normalization failure
  is not explained away by this OrthoHMM defect. Nor does the defect establish
  any comparative accuracy advantage.
- YGOB job 20917 was still PENDING, with no runtime. Cancelled only that pending
  job; accounting confirms CANCELLED by 1000, elapsed 00:00:00. Its recovered
  [batch script](ygob_cancelled_batch_20260916.sh) is retained. No YGOB accuracy
  has been inspected; inference must be resubmitted with validated binaries.
- Historical main-benchmark results are not automatically invalidated: they
  used other execution contexts and showed effective profile expansion. Their
  native binary provenance still needs explicit audit rather than assumption.

## Correction And Remaining Gates

1. Create a separate checkout of the same frozen source revision. Build native
   kernels there with recorded compiler, flags, source and binary hashes.
   Never add binaries under the checkout used by active simulations.
2. Require an isolated exact-interpreter/exact-checkout synthetic profile smoke
   test. The new replay launcher applies this before creating output or running
   expensive inference. Extend this gate and binary manifests to the simulation
   and YGOB workflows before rerunning them.
3. Repeat the label-blind historical replay in a new output directory. A working
   smoke test alone does not prove historical partition equivalence.
4. Only after equivalence, execute the guarded factorial preparation and cells.
   Do not queue it behind failed job 21088.
5. Amend execution manifests prospectively for corrected-runtime simulation
   and YGOB runs, retaining scientific settings and all prior evidence. Reuse
   validated simulated inputs and valid comparator outputs where provenance
   permits; do not silently relabel old results as corrected results.
6. Fix the broad profile exception handling in a separately tested development
   change. Do not silently alter the frozen scientific implementation mid-run.

The source-only manifests and prior native-output validators did not establish
that a requested computational stage actually executed. This is an unresolved
publication gate, not a completed runtime repair.
