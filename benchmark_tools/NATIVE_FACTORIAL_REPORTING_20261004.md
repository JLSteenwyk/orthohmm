# Native Factorial Reporting

`export_native_factorial_progress.py` creates a JSON/TSV/Markdown snapshot of
the thirteen remaining full-native identities. Supply explicit checksum-bound
terminal reviews and optional score/recovery receipts; do not discover files
by filename or borrow old cached scores for unfinished native output. Each
unreviewed identity stays `no_supplied_terminal_review`, with null metrics.
That means absent evidence in this snapshot, not a scheduler-state assertion.

The reporter checks direct plan/review/score bindings and terminal resource
scopes. It does not launch jobs, score predictions, replay raw accounting,
perform new admission, authorize the next identity or certify publication
readiness. Successful OrthoBench scores retain development exposure and
reference-group co-membership semantics. QfO resources can be reported without
scores, but OrthoBench F1 cannot be substituted for QfO native metrics.

The explicitly adopted failed-wrapper recovery remains
`failed_wrapper_science_recovered`. Its wall/CPU/peak observations are visible
failed-attempt values, never a clean-success timing. The recovery's score has
no reference-coverage numerator: retain null rather than deriving coverage
from the equality of the complete input/prediction universes.

```bash
/usr/bin/python3 -I -S -B benchmark_tools/export_native_factorial_progress.py \
  --plan benchmark_tools/results/native_factorial_receipt_amendment_20261004/plan.json \
  --plan-sha256 6c87babcbb5581830e0b9e7b9bf9aaba30a85bde1c4ab465e561017e67e9c89b \
  --attempt benchmark_tools/results/native_factorial_terminal_review_22427.json \
    c3e55f042a230e54b531f429e207a9f2928709a66695a308869dd7feff90ef13 \
    benchmark_tools/results/native_factorial_recovery_22427.json \
    0c919cfb97f497a98ca7ad6e614cb370731ddb58a63bdc4550fa0e640f4cbc15 \
  --attempt benchmark_tools/results/native_factorial_terminal_review_22428.json \
    a4d0d1c8893a70a1a08fef0528da1ad0a3e3def54e81d051388152af2dfdab19 \
    benchmark_tools/results/native_factorial_orthobench_score_22428.json \
    da7937a9bdb5a5bfd2771395d5cfa9ea7e92dc3b6f9b5f88d36b8de80ecfbe68 \
  --output benchmark_tools/results/native_factorial_progress_20261004_v2
```

Use `--attempt REVIEW SHA256 - -` for an unscored terminal attempt. Existing
destinations are refused; later snapshots retain previous versions. The
three previously associated configurations, original cached stage costs and
main scaling panel remain separate; no pooled medians, component-overhead
estimate or efficiency ranking is calculated here. Wall excludes preparation,
conversion and scoring; CPU includes the wrapper bracket; lifetime peak
includes the native-step launcher, not algorithm-only RSS.

Timing measurements were collected on a shared Threadripper while other
analyses were running. Competition for CPU, memory bandwidth and I/O may
have affected elapsed times, with an unknown and potentially tool-dependent
impact. These are observed shared-host timings, not estimates of isolated
performance. The disclosure is also generated into every Markdown snapshot.
