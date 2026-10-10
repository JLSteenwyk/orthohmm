# Native FAS Controlled Context Result

## Actual Outcome

The [corrected native execution](native_fas_context_probe_20261010_v2.json)
completed EXIT0 with 12/12 contexts, 40 pair evaluations, 34 numeric returns,
6 NA outcomes and no within-pair score or return-status differences.
Each focal pair has six observations, giving 20 comparisons against its first
observation. These are controlled annotation fixtures, not benchmark scores.

| Pair | Native Directional Strings | Native Loader Mean | Observations |
| --- | --- | --- | --- |
| A/D | 1.0 / 1.0 | 1 | 6 |
| B/E | 0.5 / 1.0 | 0.75 | 6 |
| C/F | 0.9985 / 0.9986 | 0.99855 | 6 |
| G/H | NA / NA | Omitted | 6 |
| I/J | 1.0 / 1.0 | 1 | 8 |
| K/L | 0.5 / 1.0 | 0.75 | 8 |

The unchanged native helper invoked `fas.runMultiTaxa` with the historical
bidirectional, no-config, merged-JSON, max-cardinality 40, paths-limit 15 and
pair-limit 30000 options. Contexts include four focal singleton runs, forward
and reverse full batches with one worker, rotated and reverse full batches
with two workers (hash seeds 0 and 1), and four focal-plus-companion subsets
with two workers. This changes worker settings, not proof that every pair
executed on each individual worker or every possible process assignment.

Native graph path counts were C/F=65536 and G=1125899906842624. C/F exercise
priority mode with nonidentical protein lengths; G/H has two NA directions
and is omitted by the native loader. Scores are the serialized native strings
and their loader mean, not an unrounded internal-kernel value.

## Preserved Failure And Correction

The [first attempt](native_fas_context_probe_20261010_v1.json) is retained
unchanged: EXIT1 in its first singleton, zero pair scores, and the native
error `pfam is missing in the seed annotation`. It used display-case tool
keys. Read-only inspection of native `fasInput.featuretypes` established
that it lowercases annotation-tool configuration labels. This was a fixture
schema defect, not evidence of context-sensitive scores or historical failure.

The separately named [correction](../probe_native_fas_context_lowercase.py)
records a display-case-to-lowercase mapping and uses the native configuration
for path preparation. Only tool keys changed: identities, feature values,
architectures and the context plan are preserved. It reuses the original
[fixture/capture/comparison/resource kernels](../probe_native_fas_context.py)
with their frozen SHA. No native function or child command was patched.
The first source was pushed as f8124a75 before execution; the correction was
pushed as dd10c4ce before its execution. The corrected focused/adjacent suite
passed 89 tests in 0.96s. Both sources and outcomes remain available.

## Readback And Scope

The [actual command and terminal record](native_fas_context_execution_20261010_v1.json)
contains both executions, tests and a separate reader invocation. That reader
independently parses all raw JSON texts, recomputes means/omissions, checks
all contexts and CLI worker settings, and compares repeated outcomes without
importing either producer. It also verifies the fixture-only transformation,
frozen original driver, corrected source and protocol bindings, native source
identities against the first result, and all false non-admission flags.
Readback EXIT0. Corrected JSON: 58118 bytes, SHA256
`845b4baea528f7633830e5a26ec58267da95165d8c4f79f2a0112c6063826f33`.

Execution used only the shared Threadripper, two CPUs for this driver/children,
4 GiB address-space limit per process, one numerical thread and a 120-second
child cap. This is not a whole-job peak-RAM observation or a tool timing rank.
No unrelated workload or DGX service was changed.

This finite successful-context evidence supports investigating the fixed-pair
assumption in the [prospective native sampling design](NATIVE_FAS_SAMPLING_DESIGN_PROTOCOL_20261010.md).
It does not establish invariance for all historical annotations, high-load
file races, batch crashes, missing historical pair identities or all possible
worker assignments. No historical random sample, annotation, benchmark score
or omission cause is reconstructed. No biological gene/family independence
is assumed or established. The expected random-denominator target remains
conditional on the stated design assumptions; its count/score-mean confidence
construction and joint projection have not yet been implemented or validated.
No native interval, sampling-law admission, independent biological validation
or publication readiness is claimed. All other endpoint uncertainty stays open.

A subsequent scoped join of the two already-retained population/stratum JSON
reports found all eight method identities consistent, positive k bounded by
P, c=min(M,9000), and 0<=r<=c. Full precomputed sum/P matches each retained
population mean; no large pair database or lookup was read. This supplies the
construction's input mapping, not confidence coverage. The retained audit
still reports historical parser hash identity and historical database hash
binding as false for all eight methods. Neither those provenance gaps nor
the fixed-pair assumption are repaired by matching summary arithmetic.
