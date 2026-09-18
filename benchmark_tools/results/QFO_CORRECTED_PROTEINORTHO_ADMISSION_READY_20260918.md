# Corrected Proteinortho Admission Workflow

Implemented `benchmark_tools/admit_qfo_corrected_proteinortho.py` to bind
the native graph to the frozen corrected-input plan and runtime manifests,
completed successful execution, exact command and working directory,
scheduler provenance, byte-identical input copies, final output inventory,
and retained logs/resource record.

The admission checks every recorded file hash before and after graph
validation, rejects extra/missing files in the final native output tree,
and invokes the complete-input graph ownership validator. It does not run
conversion, reference filtering or scoring, and it cannot overwrite an
existing admission record.

## Validation

- 37 focused tests passed across the admission, graph ownership and
  corrected Proteinortho runner modules.
- The real command was tested against the still-running corrected run.
  It correctly failed with `Native execution has not succeeded` before
  expensive output validation. No admission report was created.
- No inference output, frozen runner, previous converter or score changed.

After the inference finishes, run from the repository root:

```bash
/home/bizon/anaconda3/bin/python benchmark_tools/admit_qfo_corrected_proteinortho.py \
  --root . \
  --output benchmarks/work/qfo_corrected_proteinortho_admission_20260918.json
```

This is a graph provenance/representation gate, not independent proof of
algorithmic completeness or correctness. The group table is bound by its
recorded hash but is not admitted for group-based scoring. Resource records
describe shared-host inference, not matched dedicated timing. Missing or
failed output must be investigated rather than bypassing the gate.

## Progress Ledger

Previous turn: progress, with graph validator and historical check pushed
as `8d31017`. This turn: admission implementation, focused tests and a
verified rejection of premature admission.

At the scheduler poll, Proteinortho `21708`, SonicParanoid `21710`, OrthoHMM
`21706_0`, factorial assessment `21711`, reconciliation `21671_3`, and DGX
scaling `21656_14` were running. BLAST `21713` was pending resources;
OrthoFinder and dependent factorial tasks were queued. No jobs restarted.

Next: run this gate after terminal inference, inspect any failure, freeze
conversion/filtering/scoring commands and validate each stage. Continue
the factorial and dedicated scaling workflows independently. Publication
readiness remains unproven; no corrected Proteinortho score exists yet.
