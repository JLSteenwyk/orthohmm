# Independent Admission For Integrated Job 22337

Prepared while native job 22337 is still running, without inspecting a new
completed benchmark result. This implements the gates and comparison from the
[presubmission protocol](INTEGRATED_FULL_OB_PROTOCOL_20260927.md); it does not
change the method, expected endpoints or historical results.

`admit_integrated_full_ob.py` refuses any job other than 22337, any nonterminal
or unsuccessful accounting row, and a different CPU/node/memory allocation.
It revalidates the frozen run plan, launcher, 215 file records, presubmission
protocol, scheduler submission and final execution records. The observed
scheduler reports requested memory as `128G`; `128Gn` is accepted as the
equivalent per-node spelling, but `128Gc` is not.

The verifier reconstructs all eight expected stage commands/environments and
checks each start record, successful completion and exact log. It rejects
failure markers, omitted/reordered stages, changed controller/native command
scope and missing launcher resource logs. All private FASTA, reference and
low-certainty files must have exactly the original pinned contents and
inventory. These checks do not rerun inference or change a running job.

After successful terminal admission, it independently checks both installed
distribution inventories against pip reports, binds installed wheels to the
hash locks and audits package payload bytes. It then reruns all four frozen
scientific readers using the installed reader interpreter and pinned exported
source, writing separate reports. The original reports are not overwritten.

Finally, it parses both historical/current root partitions against the full
251,378-gene universe and recomputes all 70 family scores using the frozen
scorer. Each recomputed score must match its own recorded score. Differences
between the new and historical methods' outputs are reported, not rejected
or hidden. The plan and admission records are checked again afterward.

## Execution After Completion

```bash
python -m benchmark_tools.admit_integrated_full_ob \
  --directory benchmarks/work/publication_integrated_full_ob_20260927 \
  --job 22337 \
  --output benchmarks/work/publication_integrated_full_ob_20260927/independent_admission
```

The output directory must be fresh. A failed admission is evidence to
investigate, not authorization to restart native inference or loosen gates.
Thirty-three focused admission, launcher and controller tests pass. Read-only
inspection of the live job confirmed its six completed installation stages,
93 private input files and native command match the frozen plan. This is not
terminal admission and does not validate a final score.

The admission script has not yet been run on a completed full job. Package
audits retain exclusions for installer-generated metadata/bytecode and
relocated non-site payloads. Successful admission would establish same-host
workflow reproduction only, not independent biological accuracy, controlled
timing, cross-platform portability or publication readiness.
