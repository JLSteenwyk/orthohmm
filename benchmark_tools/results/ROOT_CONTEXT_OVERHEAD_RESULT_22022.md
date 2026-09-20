# Paired Root-Context Overhead Panel 22022

## Outcome

All 18 prescribed native commands validated, including raw collector replay,
provenance, output semantics and canonical equality to the same-method
22021 baseline. All nine pairs and all three complete method medians passed
the frozen engineering budgets: at most +10% for every pair and +5% for each
three-pair median. No task, pair or CPU flag was excluded.

Signed change is `root_context_wall / lineage_wall - 1`, not a comparison
between native methods. Negative values remain negative and do not establish
that extra monitoring makes a method faster. The contrast is supplementary
root context versus periodic lineage collection, not total instrumentation
overhead. Fixed ordering, drift, caches and background activity limit causal
interpretation. These budgets are not confidence intervals or scientific
timing admission criteria.

| Pair | Method | Lineage seconds | Root-context seconds | Signed change |
| --- | --- | ---: | ---: | ---: |
| 0 | OrthoHMM high sensitivity | 556.150353 | 550.980204 | -0.929631% |
| 1 | OrthoHMM satellite_v2 | 808.428589 | 819.284502 | +1.342841% |
| 2 | OrthoFinder 3.1.5 full | 619.908470 | 610.873433 | -1.457479% |
| 3 | OrthoHMM satellite_v2 | 804.030631 | 805.017766 | +0.122773% |
| 4 | OrthoFinder 3.1.5 full | 613.073274 | 615.229356 | +0.351684% |
| 5 | OrthoHMM high sensitivity | 544.793355 | 545.980729 | +0.217950% |
| 6 | OrthoFinder 3.1.5 full | 621.957447 | 619.816909 | -0.344162% |
| 7 | OrthoHMM high sensitivity | 558.196609 | 545.611659 | -2.254573% |
| 8 | OrthoHMM satellite_v2 | 793.032790 | 804.116018 | +1.397575% |

Complete method medians are -0.929631%, +1.342841% and -0.344162% for
high sensitivity, satellite_v2 and full OrthoFinder, respectively. The
machine-readable audit retains full precision, native output types/counts,
GNU-time accounting and separate native-step cgroup memory measurements.
All commands used the frozen four-proteome, 73,266-protein input.

## Environment And Flags

Slurm job 22022 completed `0:0` in 03:21:00, from September 19 21:16:17
to September 20 00:37:17, 2026, America/New_York. It used one exclusive
20-CPU, 96-GiB allocation, a five-hour limit, and no requeue/restart.
The bounded submitting SSH session passed comparison to terminal scheduler
evidence. The local recorder retained 2,405 observations without errors.
No additional DGX SSH connections or transfers occurred during the panel.

All 11,846 observed intervals remain, including 344 original and 60 narrow
flags. All raw pressure and supplementary context data remain in hashed
reports. No residual was used to assign a cause or correct a timing.

The bounded user journal covers 01:16:07 through 04:37:27 UTC on September
20. Its printed timezone is UTC-07:00. User-manager startup was recorded just
before the printed scheduler start, with no shutdown reported in the window.
The unrelated `samwise-daemon-samwise.service` scheduled 2,297 restarts within
the printed job window after missing-working-directory failures; two further
entries occur in the post-job margin. No service or persistent login setting
was changed. This is not background-free execution or causal attribution of
the CPU flags. Non-CPU isolation and prospective timing eligibility remain
unresolved; no scientific comparison is admitted by this panel.

## Reproduction And Provenance

Execution source is `272a94d`, frozen before submission in `ee3a025`.
Recipe SHA-256 is
`ef3c4d0e083e31273a098cb342fa89b797911a84f7fe643ca256b27829c40b82`.
The protocol and plan are unchanged from
[the prospective design](ROOT_CONTEXT_OVERHEAD_PROTOCOL_20260919.md).

| Retained artifact | SHA-256 |
| --- | --- |
| `root_context_overhead_audit_22022_20260920.json.gz` | `6da8441923e51a404fb2e4fec4c7302593fc1a1bc7e338b8996b9d2d2f24641a` |
| `root_context_overhead_receipts_22022.tar.gz` | `5af248d568d14be85bdf6ab686b712d222537a15988f2f33c124c4a7a00d5097` |
| `root_context_overhead_session_audit_22022_20260920.json` | `53e1f423f1b415d48ae052825ae0eb57f4d46b01e055402e829358ed54a96668` |

The 1.6-GB raw archive remains outside Git in `benchmarks/work/` and on the
DGX: `root_context_overhead_archive_22022.tar.gz`, SHA-256
`1623f2a8c186904002efea72b32b3efcb79f6ef25202f91d1e5c21951c4ad441`.
Reuse the pinned `root_context_native_inputs_22021.tar.gz` input archive,
SHA-256 `91830e6740313a74214a447db0aff379fe3913fd5a3e025af9e6c200d786dc9f`.
The journal is included in the receipt archive, SHA-256
`e4049363be9f6e5647b995bb5b0d4bfc06f37a2d3fb91a158f412205f452e964`.
These local archives do not imply external deposit or redistribution rights.

Extract the raw and input archives together into a fresh directory. Extract
the receipts under `benchmarks/work/`. Retain the baseline audit's referenced
raw evidence, which is checked independently. From the repository root:

```bash
python -B -m benchmark_tools.audit_root_context_overhead \
  --archive benchmarks/work/root_context_overhead_archive_22022 \
  --results benchmark_tools/results \
  --recipe benchmark_tools/results/dgx_root_context_overhead_recipe_20260919.json \
  --recipe-sha ef3c4d0e083e31273a098cb342fa89b797911a84f7fe643ca256b27829c40b82 \
  --scheduler benchmarks/work/root_context_overhead_scheduler_22022/scheduler_22022.txt \
  --job 22022 \
  --prior-audit benchmark_tools/results/root_context_native_audit_22021_20260919.json.gz \
  --session-directory benchmarks/work/root_context_overhead_submission_v1 \
  --scheduler-timezone America/New_York \
  --output benchmarks/work/root_context_overhead_reaudit_22022.json
```

The fresh relocated archive independently validates all 18 tasks without
panel issues. Native results, all flags, supplementary context, accounting,
observation bounds, paired comparisons and waiting-session results exactly
match the original audit. Location-dependent evidence records are excluded
from this result comparison, but were checked by each audit. The exact
excluded fields and both uncompressed audit hashes are retained in
`root_context_overhead_relocation_22022_20260920.json`. No native command was
rerun for this check. Both full raw audits remain in `benchmarks/work/`.
