# Gated Canonical QfO Phylogeny Launcher

Added `admit_qfo_raw_tree_reuse.py` and `run_qfo_canonical_phylogeny.py` for
the second arm of the frozen downstream comparison. No full canonical plan
or job has yet been created.

Admission requires both native 22329 and readback 22330 to be uniquely
COMPLETED with exit 0:0 and their expected CPU allocations. It verifies the
frozen submission, readback plan, completion receipt, exact native admission,
all eight historical/fresh scientific reports, complete input universe and
nested source/artifact identities. Valid differences from historical output
do not prevent reuse admission; they must remain visible in the result.

The launcher pins that admission, the admitted canonical candidate partition,
constraints, original sequences, environment and external tools. It copies
only input-identical raw-tree artifacts through the tested cache helper,
infers the species tree afresh, reruns reconciliation and membership filtering,
and records actual reuse and output identities. Changed families are inferred
afresh; raw-tree reuse must equal the selected cache count. No requeue/retry,
parameter optimization or scoring is introduced.

Validation: 57 focused tests pass, covering bound completion receipts,
contradictory artifact identities, scope mutations, raw-tree cache selection,
and pair comparison. A [live-state probe](qfo_canonical_gate_probe_20260927.json)
confirmed pending readback 22330 is rejected before a canonical output
directory is created. This is a negative gate test, not completed full-data
admission or native execution. The earlier installed cache fixture validates
the native reuse path; this full launcher still awaits its real prerequisites.

Next recheck jobs 22329/22330, diagnose any failure without native retry, and
prepare/pin/commit the canonical execution plan only after successful admission.
Then submit once and independently read back the full canonical outputs and
all prespecified contrasts. Existing pinned sources were not modified.
Historical accuracy scores and the DGX deferral remain unchanged.
