# Reproducible Zombi Infrastructure Smoke Test

The complete tree/genome/sequence workflow succeeds with pinned Zombi
`8db13ee4ba007f46c17f38586d31e5aa617c1647` and the explicit sequence-seed
adapter in `run_zombi_seeded.py`. This is an infrastructure result, not a
completed simulation accuracy or robustness benchmark. Native evolutionary
truth still requires independent checking and conversion.

Zombi is described in [Davin et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC7031779/)
and distributed in the [authors' repository](https://github.com/AADavin/Zombi).
The local checkout is unmodified. Its sequence engine is Pyvolve, not a
hand-written substitution simulator.

## Seed Integration

The pinned Zombi `SequenceSimulator.py` seeds Python and NumPy globally,
but calls `pyvolve.Evolver` without the explicit `seed` argument. Installed
Pyvolve 1.1.0 creates a new NumPy generator in that situation, so the native
parameter seed alone does not establish reproducible sequence simulation.
The exploratory 1.0.3 package also reseeds without an explicit seed; it was
not retained for the final smoke test.

The adapter subclasses the Pyvolve evolver and supplies an integer from
the first 128 bits of SHA-256 over `master_seed:sequence_file_basename`
when no explicit seed is already present. It preserves every other argument
and any explicit seed. Basenames make results independent of output-directory
paths. This is a documented integration change, not unmodified native Zombi
sequence execution. The tested modes are T, G, and S only; no claim is made
about advanced simulator modes.

Source revision, tracked Python cleanliness, Pyvolve version, parameter/master
seed equality, and fresh stage paths are checked before execution. Existing
stages are rejected because native Zombi can delete prior outputs. Future
datasets must use unambiguous, unique family filenames, and must continue to
record the adapter version as part of simulation provenance.

## Observed Results

The small panel uses four extant species, ten initial genes, WAG protein
evolution, sequence length 100, scaling 0.2, genome-level duplication rate
0.2 and loss rate 0.1, and zero transfer/origination/rearrangement rates.
Unchanged native defaults, including event extension parameters, are retained
in each generated parameter file and their hashes. These settings are a
feasibility test and are not a frozen scientific robustness panel.

| Run | Master seed | Recorded products | Duplication events | Loss events | Terminal F events |
| --- | ---: | ---: | ---: | ---: | ---: |
| repeat_a | 20260916 | 84 | 1 | 0 | 41 |
| repeat_b | 20260916 | 84 | 1 | 0 | 41 |
| independent | 20260917 | 84 | 2 | 1 | 44 |

All 84 stage files match byte-for-byte between repeat_a and repeat_b.
The independent seed produces different sequence files. Each run emits 20
FASTA files, including complete and pruned outputs; these are not 20
independent gene families. Complete output includes ancestral sequences and
must not be passed directly to orthology inference as an extant proteome.
Native event-code counts above are descriptive and do not replace a validated
truth adapter.

The committed [machine-readable report](zombi_smoke_20260916.json) records
all commands, generated/default parameter hashes, output hashes, driver
hashes, source pin, Python and package versions. Raw files and logs are in
`benchmarks/work/zombi_smoke_v2`. The first smoke run remains in
`zombi_smoke_v1`; it used 32-bit family seeds and is superseded by the
128-bit adapter test, not silently overwritten.

## Environment And Reproduction

The isolated overlay virtual environment uses Python 3.10 with Pyvolve 1.1.0,
ETE3 3.1.3, NumPy 2.2.6, SciPy 1.15.3, Biopython 1.86, and NetworkX 2.8.8.
It inherits system-site packages and is therefore not yet a portable fully
locked environment. The base benchmark Python environment was not changed.
No complete environment-variable dump is recorded; the workflow sets
PYTHONHASHSEED=0 and BLAS/OpenMP thread limits to one.

```bash
benchmarks/work/zombi_env_v1/bin/python benchmark_tools/smoke_zombi.py \
  --source benchmarks/work/Zombi_publication_candidate \
  --output benchmarks/work/zombi_smoke_v2
```

The existing output directory is preserved and the command refuses to
overwrite it; choose a new directory for reproduction. Five focused tests
cover stable family seeds, invalid inputs, preserved explicit seeds, and
unchanged simulation arguments. The real nine-stage integration run also
finished successfully and verified exact repeated output hashes.

## Remaining Simulation Gates

Validate extant sequence IDs against genomes, event histories, gene-tree
leaves and reconciled XML; derive orthology from event-labeled ancestry,
not family cliques. Include singleton/extinct-family behavior and known
small truth fixtures. Freeze multi-seed scientific conditions before scoring,
including duplication/loss, divergence, missing data, and uneven taxon
sampling. Assess gaps/fragments separately because this tested sequence
path does not introduce indels. Then run matched methods, error analyses,
tree/parameter robustness, and cost measurements. None of these downstream
scientific requirements is fulfilled by this smoke test alone.
