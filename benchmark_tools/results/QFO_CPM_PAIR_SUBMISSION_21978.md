# CPM Native Pair Conversion

Submitted September 19, 2026 at 10:33 EDT as serial array `21978_0..1`.
Both tasks were confirmed PENDING on dependencies. Each requests two CPUs,
64 GiB and four hours on `bizon`, without requeue. Dependencies are
`afterany:21976,aftercorr:21976`: require terminal whole-array accounting
and successful matching native admission. Index 0 is `cpm_low`, index 1
is `cpm_high`. No CPM dataset output or accuracy is admitted yet.

Executor commit `ba4005be49a050c538306f2d7e64a0c4ac387eeb` is pushed and
frozen in `benchmarks/work/publication_qfo_cpm_pairs_v1`.

| File | SHA-256 |
| --- | --- |
| `benchmark_tools/prepare_qfo_cpm_pairs.py` | `53124b59b10ee111151d2b4332f643d69306c9c2915736dbcca4a848b6e3a197` |
| `benchmark_tools/results/qfo_cpm_pairs_batch_20260919.sh` | `fac98ed71630b2938380cd7433c4316baa2ff45da1fe1f0682d4e278582a3972` |

```bash
sbatch --parsable benchmark_tools/results/qfo_cpm_pairs_batch_20260919.sh \
  /mnt/ca1e2e99-718e-417c-9ba6-62421455971a/ORTHOHMM/orthohmm/benchmarks/work/publication_qfo_cpm_pairs_v1 \
  ba4005be49a050c538306f2d7e64a0c4ac387eeb
```

The converter pins the native-admitter revision and source, checks the
matching CPM context/candidate admission and unscored status, rehashes
inputs and reruns the frozen native validator. Its fresh report must equal
the admitted report exactly. Native pair fields/species ownership are
validated by the shared streaming converter; no RootHOG clique expansion
is substituted. Frozen reference mapping must retain every pair. Mapping
losses and all other post-preflight failures preserve partial outputs,
counts and a failure report. Existing output directories are not reused.

Outputs will be in `benchmarks/results/qfo_cpm_pairs_v1/<arm>/`, including
`pairs.tsv`, `pairs.qfo.tsv`, conversion counts, fresh admission evidence,
and `results.json`. Successful status is
`cpm_native_pairs_prepared_unscored`; accuracy and publication flags remain
false. Participants are `ohmm_qfo_parameter_cpm_low` and
`ohmm_qfo_parameter_cpm_high`. Assessment and score admission are still
required, followed by the prespecified SwissTrees parameter comparison.

112 focused tests pass in 0.83 seconds: new conversion orchestration,
parameter conversion, CPM native admission, native factorial conversion,
and QfO native-pair conversion. Tests exercise actual pair/mapping files
with mocked scheduler/frozen-admitter calls, including mapping loss,
changed inputs, mismatched fresh reports, source/revision mismatch and
unfinished admission rejection. A fixture retry initially retained its
injected changed-input mock; the fixture was corrected before the passing
run. No production behavior was relaxed. Batch `bash -n` and staged
`git diff --check` passed. Production conversion remains pending.
