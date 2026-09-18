# First Expanded-Candidate QfO Assessment Admitted

Cell `p0_c1_r0` completed scoring as job21681, COMPLETED0:0 on bizon with
8CPUs; scheduler elapsed39:15. Independent admission21683 completed0:0
in13seconds. These are workstation scoring/validation durations, not matched
inference timings or DGX scaling results.

The cell has profile refinement off, candidate expansion on and reconciliation
off. Initial HMM search remains present. QfO predictions are cross-species
clique pairs derived from complete expanded candidate groups, not native
phylogenetically inferred pairs. Of11,388,191 generated pairs,11,351,618
survive the frozen mapping;36,573 are excluded by that mapping.

## Validation

The existing frozen independent validator checked successful execution and
exact cell/job identity, pair-conversion provenance, frozen scorer/environment,
before/after input and output identities, all15native Nextflow tasks and
48native metric records (12aggregate metrics plus36SwissTrees family metrics).
All six aggregate endpoints agree with the native metric records.

Repeated the same validator from its pinned executor after admission:
`admission_2_recheck_20260917.json` is byte-identical to `admission_2.json`
(`cmp` exit0). No scoring or inference was repeated. The committed snapshot
`qfo_factorial_assessment_p0_c1_r0_20260917.json` has SHA-256
`01ebf5bc0043aac5269e1398610614ffde6bb32e9d84833e30aa4ee37505faaf`.
It retains scheduler accounting, pair counts, native metrics, task records,
execution/source/environment hashes and admission limitations.

## Retained Outcomes

Native/project endpoint values on the0-to1scale are GO0.471009680,
EC0.919343710, VGNC harmonicF1 0.647990059, SwissTrees harmonicF1
0.676019536, TreeFam-A harmonicF1 0.604915207 and FAS0.746946524.
The project-defined secondary six-metric mean is0.677704120; it is not an
official QfO F1. Full native axes and error bars remain in the JSON snapshot.
Native error bars are not paired configuration-difference intervals; FAS
retains its sampled-pair limitations.

This is an admitted cell, not the complete factorial result. Three of eight
cells now have admitted scores: two exactly reused unexpanded baselines and
this fresh expanded cell. The profile-on expanded R-off assessment21682 is
running, with validator21684 queued; four R-on cells still await reconciliation,
pair conversion and scoring. No inference settings, endpoints or statistical
protocol were changed after this outcome. The full42-endpoint SwissTrees
analysis must wait for all eight independently admitted cells and the exact
reference-count audit. Missing cells are neither imputed nor silently omitted.
