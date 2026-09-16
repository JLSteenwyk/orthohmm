# Evolutionary Simulator Preflight

Status: infrastructure smoke test failed before simulation; no synthetic
accuracy, biological validation, or robustness result was produced.

## Candidate Selection

ALF supports genome evolution with substitutions, indels, duplications,
losses and other events, and emits evolutionary truth. This makes it a
candidate for the required simulation controls rather than a substitute
based on arbitrary sequence mutation. The evaluated source is the authors'
[repository](https://github.com/DessimozLab/ALF), pinned to
`b674ab10018c3c0fcc0806434dddf7b36136e2c1` in the ignored work directory
`benchmarks/work/ALF_publication_candidate`.
The model description is [Dalquen et al. (2012)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3341827/).

The test targets four species, ten root proteins, WAG substitution,
fixed seed 20260916, and no indels, duplications, losses, transfer, or
rearrangement. It is an infrastructure test, not the final scientific panel.
Successful output would still require ID, family, orthology, and sequence
truth validation plus a repeated-seed reproducibility check. In particular,
a no-duplication smoke test cannot validate duplication/loss scenarios.

## Observed Compatibility Failure

The existing QfO Darwin container reports Darwin 4.0 (2019-02-06). Running
the unmodified ALF source through that engine exited 1 after loading several
library modules, with `The program terminated incorrectly: Bad VectorABQMode`.
The ALF parameter-read and simulation-completion messages were never reached.
No scientific parameter or result was selected on the basis of this failure.

This is evidence against this exact source/engine combination being usable
as configured, not evidence that ALF is generally invalid. The next step
is to investigate a compatible official ALF distribution or evaluate another
published simulator with explicit gene-history truth. Do not patch the
biological algorithms merely to make this smoke test pass. Zombi is a
candidate alternative with a [published description](https://pmc.ncbi.nlm.nih.gov/articles/PMC7031779/)
and [author-maintained source](https://github.com/AADavin/Zombi); it has not
yet been installed or tested here.

## Reproduction

From the OrthoHMM repository root, after checking out the pinned ALF source:

```bash
mkdir benchmarks/work/alf_smoke_v1
singularity exec --cleanenv \
  --bind "$PWD/benchmarks/work/ALF_publication_candidate:/alfroot/simulator:ro" \
  --env DARWIN_REPOS_PATH=/alfroot \
  --env "ALF_SMOKE_PARAMS=$PWD/benchmark_tools/alf_smoke_params.drw" \
  --env "ALF_SMOKE_OUTPUT=$PWD/benchmarks/work/alf_smoke_v1/" \
  qfo_benchmark/scoring/container_cache/qfobenchmark-darwin-2022.1.img \
  darwin -E < benchmark_tools/alf_smoke_driver.drw \
  > benchmarks/work/alf_smoke_v1/native.log 2>&1
```

The existing directory is preserved; use a new output directory for another
attempt and change the output environment variable and log destination
consistently. The source checkout was not edited. Inspect native error
messages as well as process status before accepting any future completion.

| Artifact | SHA-256 |
| --- | --- |
| `benchmark_tools/alf_smoke_params.drw` | `ffe637b6c4741efbb45d2e9b2adb588fa06b97635c82de04f3795862c8e183d3` |
| `benchmark_tools/alf_smoke_driver.drw` | `ba3470ae0be1d14925ca29814de750592b9f8cbd21b3f4d454cbf8b8e0883369` |
| [Native log](alf_smoke_native_20260916.log) | `c278b87f26d8f8eadd07cc5a642e1c2016868ba5fb5ce8409cd9e9afe275294f` |
| QfO Darwin container | `abed2b0ff8bb033a2372d0cd5282fd1d6c8ab5fa07496afc48c82f54ed65ea21` |

The simulation work package remains open: a compatible validated engine,
truth adapters, multiple independent seeds, duplication/loss and divergence
conditions, missingness and uneven sampling, and matched method evaluations
are all still required.
