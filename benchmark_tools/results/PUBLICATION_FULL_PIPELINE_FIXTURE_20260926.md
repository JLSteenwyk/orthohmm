# Explicit Canonical Full-Pipeline Fixture

The new experimental entrypoint runs the installed frozen method from raw
FASTA through built-in HMM search, high-sensitivity grouping/profile refinement,
satellite-v2 candidate expansion and inferred species/gene phylogeny. It
canonicalizes directed hits only at the candidate expansion boundary, without
changing hit values, earlier search/grouping stages or installed scientific
sources. It checks the recovery dependency versions, frozen orchestration/
refinement source hashes and ordering-policy hash, rejects existing outputs,
and records one policy invocation with original/canonical array identities.
It is not a production default or a generally portable release.

The [executed result](publication_full_pipeline_fixture_20260926.json) covers
16 genes in four species, 96 candidate-input hits, three root groups and 36
ortholog pairs. Two families build the inferred species tree; one family is
reconciled, with four duplications and three speciations. All checkpoint reuse
counters are zero. Four independent scientific readers pass; 86 focused tests
pass. A post-run complete package audit has the identical hash to the earlier
recovery-install audit. No benchmark accuracy score was calculated or changed.

This fixture has an empty satellite merge trace. It validates invocation of
canonical ordering and the complete execution path, but not full-data merge
behavior. Earlier separate order-mechanism and full OrthoBench recovery
experiments are not substitutes for validation of this new entrypoint on a
merge-heavy dataset. A pinned full OrthoBench reproduction remains next.

## Preserved Failures

The first launch stopped before inference because the frozen argument
processor requires an existing output directory and raises SystemExit(0)
otherwise. Its failure receipt, log and launcher source are retained. The
corrected launcher creates a new inference directory and converts early
SystemExit into RuntimeError, so it cannot report process success that way.
The corrected attempt used a new directory, not an automatic retry or resume.

The first readback rejected a path to an empty constraint trace because no
membership constraints had been applied. Its source and partial reports are
retained. The corrected reader pins the empty trace but passes no constraints
to the existing event oracle. It then passed without rerunning inference.

## Invocation

With the verified recovery venv, MAFFT/FastTree and input directory available:

```bash
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  PYTHONHASHSEED=0 MAFFT_BINARIES=/absolute/mafft/libexec/mafft \
  /absolute/recovery/venv/bin/python -I \
  /absolute/repo/benchmark_tools/run_publication_pipeline.py \
  --input /absolute/proteomes --output /absolute/new/run --cpu 32 \
  --aligner /absolute/mafft/bin/mafft --tree-builder /absolute/FastTree
python -m benchmark_tools.audit_publication_pipeline \
  --directory /absolute/new/run --output /absolute/new/readback
```

The executed fixture used one CPU and the same external MAFFT/AVX2 FastTree
as the earlier full recovery experiment. Native timing is shared-host fixture
timing, not controlled comparative performance. Source/argument/tool/input
records are retained in started.json, and completion remains pending scientific
readback until the second command passes. The wrapper is a benchmark harness
requiring the checkout and a separately verified installation; it is not a
standalone installation recipe or a full dependency attestation by itself.
