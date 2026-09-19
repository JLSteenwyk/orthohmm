# Full-Node Control Description

Completed the prespecified descriptive comparisons from the pinned nine-trial
audit. Retained all189 intervals and separately summarized the19 fully
enclosed common-work intervals per trial. No interval was relabeled or
removed. These are engineering descriptions, not independent-interval tests.

## Common-Work Medians

Ranges below span the three trial medians for each condition, not confidence
intervals. The complete per-trial min/median/max values and all interval
records are in the machine-readable report.

| Quantity | Steady | Churn | Contended |
| --- | ---: | ---: | ---: |
| Signed residual, average cores | 0.05019-0.05241 | 0.13085-0.13993 | 0.99949-1.00133 |
| Native average cores | 19.94775-19.94978 | 19.87831-19.88481 | 18.99912-19.00026 |
| Host process creations per window | 0 | 9,617-9,944 | 0 |
| Host context switches per window | 424-454 | 24,910-25,915 | 571-589 |
| Native CPU pressure some, microseconds | 2,710-2,857 | 43,542-46,877 | 49,986-50,235 |
| Observer/batch CPU seconds | 0.04905-0.05112 | 0.05011-0.05046 | 0.99811-1.00034 |

Within-block churn-minus-steady median residual differences are
0.085411,0.078448 and0.087623cores. Contended-minus-steady differences are
0.949300,0.947385 and0.949024cores. All signs are retained. Other
within-block metrics are recorded without hypothesis tests or pooling
intervals as independent replicates.

## Interpretation

This bounded process-creation workload increases context switching, pressure
and the residual, but its narrow residual remains below the unchanged
threshold. Creation rates at roughly9.6k-9.9k per host window do not reproduce
the historical native-tool flags. Conversely, the known competitor produces
a roughly0.95-core residual increase and matching reduction in native CPU.
These comparisons do not identify the cause of historical flags.

Native pressure increases in both churn and contended controls, despite
different competitor status. Therefore pressure alone is not a specific
foreign-workload detector. The named outside-target frontier excludes the
observer's batch step; its near-zero value during contention does not mean
the known batch competitor was absent. The hierarchy's batch-step measure
records that competitor. Do not conflate these scopes.

Host/frontier/pressure measurements use distinct windows, so descriptive
differences are not a causal accounting subtraction or runtime correction.
Whole-worker creation totals/rates are reported separately from common-window
host counters. None of these quantities grants scientific timing admission.

## Evidence And Next Step

Report: `full_node_control_description_21918_20260919.json`, SHA-256
`a0629de5fab89e1eb263e57e79bba394c0b6dfb71c393ab91f88435221c9de2e`.
Source audit SHA-256:
`bd8e11ead6ee2e979784c8cd97730d3cf82ee62e4343c9cdf18c0475d2221cd3`.
The script rechecks the complete frozen trial order, audited common-window
indices and original completion-witness hashes before and after reading.
Twenty-eight focused description tests pass, including nine new cases.

The full-node control experiment and its prespecified descriptive reporting
are now complete. Do not repeat controls merely to seek an explanation that
fits. The remaining timing work is to quantify collector overhead and assess
the native residual at whole-command and longer aggregation scales, retaining
the original interval flags. A prospective scientific inclusion policy must
distinguish observed interference from non-atomic accounting uncertainty;
neither silently relaxing thresholds nor requiring an unexplained diagnostic
to identify a cause is justified by these controls alone. No new inclusion
rule is adopted here, and the27 matched scaling runs remain outstanding.
