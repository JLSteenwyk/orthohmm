# DGX Periodic CPU Flags

## Evidence

`summarize_pressure_panel_flags.py` describes the already audited reports;
it does not replace the raw-counter replay or admit scientific timings.
The input compressed audit SHA-256 is
`d656ae37dcb64d617a132c82391745e33bfd5f078f6c1df7186156f0022ea006`.
All eight report files were checked against their audit hashes before and
after reading. Identical duplicated evidence records are permitted;
conflicting records or multiple report paths fail validation.

Output: `dgx_pressure_flags_21889_20260919_v2.json`, SHA-256
`d560fa5d42458cd7914a7ce93dfc87f64ae2d19784a0f26d1c13a8a6ac255781`.
All 18 task statuses are retained, including failed tasks 5 and 12.
The earlier uncommitted JSON without overhang summaries is exploratory,
not the retained result. Thirteen focused tests pass, including an
end-to-end compressed audit fixture and evidence-tampering rejection.

## Results

All flagged intervals have only `excess_unassigned_cpu`; none has a
negative-accounting or steal-time flag. Every whole-command screen passes.
Values below are rounded; residuals are signed average cores, not identified
foreign CPU usage. Overhang and outside-frontier columns are medians.

| Task | Method | Flagged / intervals | Interval residual median | Whole-command residual | Read overhang (ms) | Outside-frontier CPU (ms/interval) |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 1 | High sensitivity | 358 / 548 | 0.399903 | 0.053815 | 18.884 | 0.887 |
| 3 | Satellite v2 | 527 / 815 | 0.379417 | 0.071294 | 13.651 | 1.108 |
| 6 | Satellite v2 | 493 / 820 | 0.313992 | 0.068725 | 12.052 | 1.083 |
| 8 | OrthoFinder full | 177 / 617 | 0.077034 | 0.037788 | 9.525 | 1.195 |
| 10 | High sensitivity | 380 / 547 | 0.417606 | 0.054306 | 18.812 | 0.840 |
| 13 | OrthoFinder full | 175 / 619 | 0.076466 | 0.037265 | 9.467 | 1.334 |
| 15 | High sensitivity | 333 / 546 | 0.291920 | 0.047933 | 13.441 | 0.879 |
| 17 | Satellite v2 | 554 / 790 | 0.398846 | 0.074009 | 18.397 | 1.050 |

Native CPU PSI `some` totals span 1.58-1.78 seconds for high sensitivity,
23.70-23.84 seconds for satellite and 2.84-2.88 seconds for OrthoFinder.
These are pressure durations for the native step and descendants, not
CPU-seconds, a wall-time correction, or a measure of external interference.
Full CPU, I/O and memory pressure totals are retained in JSON.

## Interpretation and Next Action

`probe_interval_cpu.interval` subtracts native-step CPU between native
reads from host CPU measured across wider outer brackets. The host windows
overlap across adjacent intervals, while the native differences telescope.
Consequently, neither summing these residuals nor treating them as foreign
CPU estimates is valid. Whole-command and interval screens also use
different native scopes (worker descendants versus native step aggregate).

The observed read overhangs and small outside-frontier counter differences
make bracket mismatch a concrete measurement hypothesis, not an established
causal explanation. Frontier counters have their own non-atomic windows;
they cannot simply be subtracted to repair the original screen. Counter
accounting delay and non-CPU interference remain unresolved.

Next test a separately versioned narrow-bracket collector under controlled
known workloads, retaining outer brackets for comparison and demonstrating
response to injected contention. Prespecify its measurement policy before
new overhead or scientific runs. Do not change this panel's thresholds,
discard failed tasks, or retrospectively admit timings. The complete
OrthoFinder overhead panel remains unavailable, and 0/27 scientific
scaling runs are admitted.

## Reproduction

```bash
python -m pytest -q tests/unit/test_summarize_pressure_panel_flags.py
python -m benchmark_tools.summarize_pressure_panel_flags \
  --audit benchmark_tools/results/dgx_pressure_overhead_audit_21889_20260919.json.gz \
  --audit-sha256 d656ae37dcb64d617a132c82391745e33bfd5f078f6c1df7186156f0022ea006 \
  --output /new/path/pressure_flags.json
```

Requires the unchanged local evidence archive referenced by the audit.
