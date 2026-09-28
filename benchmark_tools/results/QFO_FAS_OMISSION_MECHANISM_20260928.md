# Native FAS Omission Mechanism

The [requested-sample audit](QFO_FAS_SAMPLE_ATTRITION_20260928.md) established
missing saved scores, not their cause. This follow-up executes a controlled
mechanism probe inside the retained QfO FAS container, which contains
greedyFAS 1.18.7. The image SHA-256 is
`62b8f97bca23b67d451384f94108a5510d5fafc1bc398f95ab93532d324563fe`,
matching the frozen assessment environment. No historical inference or scoring
run is repeated with changed parameters.

## Confirmed Code Path

The historical command passes `--paths_limit 15`, which the native FAS option
conversion turns into a limit of `10**15` feature paths. In `fc_prep_query`,
path counts strictly above the limit return `None`. The caller `fc_main`
then returns a tuple of `NA` values for the pair instead of numeric scores.
The native bidirectional JSON writer serializes the result as `['NA', 'NA']`.
QfO's `load_precomputed_fas_scores` catches the numeric-conversion failure and
omits that pair. A mixed numeric/`NA` pair is also omitted rather than scored
using its one available direction. Two numeric directions remain present and
are averaged normally.

The path limit is distinct from `--max_cardinality 40`: exceeding the latter
can activate the priority strategy, whereas exceeding the path limit rejects
the protein before that strategy is selected. The additional `--pairLimit
30000` cannot by itself reduce a 9,000-line requested sample in the inspected
native `create_jobs` implementation, which samples only when the input count
exceeds that limit. This does not rule out other omission mechanisms.

## Executed Controls

Four injected-count boundary checks confirm admission just below and exactly
at `10**15`, rejection just above it, and admission with the native internal
disabled-limit value of 1. These checks mock graph preparation only; the cutoff
branch itself is native.

Two additional synthetic annotation fixtures do not mock graph construction.
Each has two mutually overlapping domain alternatives per nonoverlapping
region. Native construction counts `2**49 = 562949953421312` paths for 49
regions and admits that architecture. With 50 regions, native construction
counts `2**50 = 1125899906842624` paths and rejects the architecture. Native
pair processing then produces the expected `NA` result. Serialization and
loading controls verify both-direction `NA` omission, mixed-direction omission
and the retained numeric mean (0.2 + 0.8)/2 = 0.5.

The [probe result](native_fas_omission_probe_20260928.json) records source hashes,
boundary cases, native synthetic counts and loader results. A second execution
reproduces the entire JSON exactly. Its [execution receipt](native_fas_omission_execution_20260928.json)
records the command, terminal success, empty stderr and image/environment pins
checked around execution. Container-internal source paths refer to that image,
not host files. The probe modifies functions only in its own Python process
and writes temporary synthetic JSON; it does not modify the container or
benchmark inputs.

## Interpretation

This establishes one native mechanism by which a successful task can omit
requested scores. It does **not** establish that any specific historical
protein exceeded the path limit, that all 1,252 OrthoMCL omissions arose this
way, or that omitted scores are missing at random. It is not a remedy or a
recalculated benchmark result. Attribution still needs the original sampled
pair list/intermediate JSON or an independently justified reconstruction.
The previous conditional completion bounds and uncertainty limitations remain.

```sh
singularity exec --cleanenv \
  qfo_benchmark/scoring/container_cache/qfobenchmark-fas_benchmark-2022.1.img \
  python benchmark_tools/probe_native_fas_omission.py --output NEW_PROBE.json
```
