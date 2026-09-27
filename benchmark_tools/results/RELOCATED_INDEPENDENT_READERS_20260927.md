# Relocated Independent Scientific Readers

The [fresh relocated native fixture](FRESH_RELOCATED_RUNTIME_20260927.md)
previously depended on validators executed from the original checkout.
An exported 31-module static analysis-source closure now runs the same four
independent readers outside the checkout. Source bytes come from immutable
revision `487d759d29ba2e88eb83346690d7fde047dddfda`; eight exporter tests pass.
Export contents were rechecked after execution.

The [execution/readback receipt](relocated_independent_readers_20260927.json)
binds both attempts, the export manifest, source identity, five report
comparisons and trace files. Structure, sequence content, event/pair semantics
and hierarchy checks all pass on the existing 16-gene relocated native output.
All scientific fields in the four reader reports and final summary exactly
match the earlier readback. Source/report identities necessarily change with
relocation and are recorded separately; 136 referenced-file occurrences were
rehashed successfully. No inference or accuracy scoring was rerun.

## Isolation Observation

The first successful run used isolated Python with `-I -B`, but its file trace
showed startup probes of unrelated editable-install paths on the original
mount. That result is preserved and does not satisfy a no-mount-access claim.
The second run additionally used `-S`, disabling site initialization, and
explicitly inserted only the exported module directory and the known installed
dependency directory. Its file-syscall trace contains no literal original-mount
prefix. No unrelated installation or startup file was changed.

Native input: `/tmp/orthohmm-fresh-runtime-20260927/native`.
Export: `/tmp/orthohmm-independent-readers-20260927`.
Successful stricter readback: `/tmp/orthohmm-independent-readback-nosite-20260927`.

To reproduce the export into a fresh directory:

```bash
python -m benchmark_tools.export_publication_readers --repo . \
  --output /absolute/new/readers
```

Run the exported `audit_publication_pipeline.audit` under `-I -S -B` with
the export root and a trusted dependency directory explicitly inserted into
`sys.path`; pass the existing native run and a fresh report destination.
The recorded dependency directory is
`/home/bizon/anaconda3/lib/python3.10/site-packages`; the tested interpreter
is `/home/bizon/anaconda3/bin/python`. Dependencies include Biopython 1.86
and DendroPy 5.0.8. These shared installations are not newly packaged here.

## Scope

This closes the checkout-code dependency for this fixture's scientific
readback. It is not a new dependency environment, OS/architecture validation,
full-dataset relocation, proof of arbitrary helper entrypoints, filesystem
sandbox or complete runtime archive. The fixture has an empty satellite-merge
trace. Static import closure and literal-prefix tracing have bounded coverage;
no universal isolation, biological accuracy or publication-readiness claim.
