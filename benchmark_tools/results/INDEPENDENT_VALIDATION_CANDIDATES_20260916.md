# Independent Validation Candidate Audit

Status: data acquisition and overlap audit only. No candidate accuracy scores
have been calculated. This document is not a completed evaluation freeze.

## Development Exposure

- All 70 Open OrthoBench RefOGs are exposed, including both parts of the
  historical development/validation split. Parameter-selection reports used
  both parts; the labels cannot now support independent confirmation.
- QfO 2020's six retained metrics have been repeatedly inspected. The earlier
  satellite_v2 freeze preceded its first QfO result but does not make future
  tuning against these outcomes independent.
- Three Kingdoms BUSCO scores and the sample proteomes are exposed.
- Synthetic holdouts used to select reconciliation or refinement rules are
  development evidence, even when their seeds differed from earlier trials.
- This inventory is a lower bound on exposure, not proof that all other
  datasets or families are unseen. Broader history and resource audits remain.

## Candidate Decisions

| Candidate | Decision | Reason |
| --- | --- | --- |
| Open OrthoBench historical validation split | Not independent | Already inspected and used for selection. |
| Updated QfO alone | Insufficient for family-independent validation | Changing a proteome release does not remove shared species, sequences, reference families, or annotations. |
| paraBench | Not selected as the independent test | Its documented inputs derive from QfO 2019; overlap needs explicit treatment, not a new benchmark-name assumption. |
| AYbRAH | Secondary candidate only | Manual curation is valuable, but its construction used OrthoMCL and OrthoDB. Automatically assigned records cannot be treated as independent truth for those methods. |
| YGOB v7 | Preferred candidate for further auditing | Curated homology and syntenic evidence, matched sequences, and many yeast taxa not in the current primary benchmark inputs; not automatically family-disjoint. |

Sources: [paraBench repository](https://github.com/rderelle/paraBench),
[AYbRAH paper](https://doi.org/10.1093/database/baz022),
[YGOB paper](https://doi.org/10.1101/gr.3672305).

## YGOB Snapshot

Official archive: `http://ygob.ucd.ie/data/v7-Aug2012/`.
The [maintainer information page](http://ygob.ucd.ie/browser/info.html)
recommends v7 for publication while v8 homology curation remains unfinished.
The site works over HTTP from this machine; HTTPS retrieval failed. Recorded
SHA-256 values identify our downloaded bytes, not authenticated transport.
Raw files are excluded from Git; their checksums are in
`ygob_overlap_20260916.json` under `candidate_inputs`.

Reproduce acquisition from the repository root:

```bash
mkdir -p benchmarks/work/independent_ygob_v7
curl -fS -o benchmarks/work/independent_ygob_v7/AA.fsa http://ygob.ucd.ie/data/v7-Aug2012/AA.fsa
curl -fS -o benchmarks/work/independent_ygob_v7/Pillars.tab http://ygob.ucd.ie/data/v7-Aug2012/Pillars.tab
curl -fS -o benchmarks/work/independent_ygob_v7/README http://ygob.ucd.ie/data/v7-Aug2012/README
```

Verify the downloaded hashes against the recorded manifest before reuse;
do not silently replace this snapshot with changed upstream data. Redistribution
permissions remain to be checked before packaging raw files in an archive.
The audit JSON records the full command for `audit_ygob_overlap.py` and the
checksums of all development input files examined.

## Observed Structure And Overlap

- 14,101 pillar rows, 107,277 ON protein records in 20 species, and 7,389 OFF
  proteins. All ON proteins map to at least one pillar; three have internal
  stop characters and require a prespecified input policy.
- Two genes (`YOR011W`, `CAGL0F01419g`) occur in both rows 113 and 9896.
  Do not silently overwrite these memberships or treat the raw reference
  as an unambiguous partition.
- 5,813 ON proteins exactly match QfO 2020 input sequences, including 5,567
  S. cerevisiae proteins. Three Kingdoms matches 5,820, including 5,583 from
  S. cerevisiae. Counts include different candidate genes with the same sequence.
- There are no exact matches to the examined Open OrthoBench and test-sample
  inputs. This is not evidence of absence of homologous family overlap.
- Matching ignores case and one terminal translation stop only. Internal
  stops, sequence substitutions, and truncations are not normalized away.

## Required Before Evaluation

1. Complete the exposure inventory, including taxa synonyms and additional
   historical inputs. Exclude development-exposed taxa from the proposed
   confirmation set; S. cerevisiae is already known to be exposed.
2. Audit homologous family overlap with primary reference resources and
   dependencies used by competitors. Distinguish clade-transfer evidence
   from family-disjoint evidence; do not promise the latter without support.
3. Define the target biological relation. The v7 README says the two positions
   for post-WGD species are arbitrary, not assigned A/B tracks. Pillar clique
   expansion is therefore not a valid truth set for resolved pairwise orthology.
   A root-level homolog-group test may be appropriate, but must explicitly
   address the ancestral event and the method output level.
4. Prespecify ambiguous-pillar, OFF-feature, non-protein, and internal-stop
   handling from data properties alone. Report exclusions and coverage.
5. Pin the method configuration, comparator commands, species list, reference
   construction, primary statistic, uncertainty method, and failure criteria
   in a committed protocol before inspecting predictions or scores.
6. Reserve a new confirmation dataset if outcomes lead to method changes.

No outcome-dependent filtering or method selection has been performed here.
