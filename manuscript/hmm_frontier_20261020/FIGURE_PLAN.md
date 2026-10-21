# OrthoHMM2 Comparative Manuscript Figure Plan

Reporting analysis dated 20 October 2026. This is a post hoc presentation of
development-exposed results, not a preregistered hypothesis test or new accuracy
experiment. Preserve the completed campaign and its original stopping rules.

Editorial revision on 21 October 2026 uses OrthoHMM2 for the extended framework.
Only display names change; historical run identifiers, plotted values, method
inventory, endpoints and reporting decisions are unchanged.

## Evidence And Comparison Set

Use the final campaign archive delivered at
`02426dce9922e03e00773af20a23afc371fc1125`. Its 29 rows include repeated controls,
historical configurations and failed attempts, not 29 tools. The main figures
include all 12 distinct complete-vector configurations: four current OrthoHMM2
configurations (unweighted canonical, Henikoff, clustered Henikoff, projection
repair), retained sensitive and phylogenetic OrthoHMM2, full and sequence-only
OrthoFinder 3.1.5, SonicParanoid 2.0.9, ProteinOrtho 6.3.6, FastOMA 0.3.5 and
OrthoMCL 1.4. Collapse only documented repeated canonical controls. Retain all
29 original rows, including explicit NA, in Supplementary Table S1.

The point-estimate primary frontier maximizes OrthoBench, VGNC and SwissTrees
F1. Report pairwise projections separately from three-objective membership.
A sensitivity analysis maximizes all five F1 endpoints, adding supporting
TreeFam-A and supplementary Three Kingdoms. Never average endpoints, use ranks
as accuracy, select a different winner per dataset, or call functional means
orthology accuracy. Exact ties do not dominate one another. Use unrounded
numbers, not displayed rounding. A point-estimate frontier is descriptive:
membership does not establish statistically distinguishable superiority.

## Original Figures

1. HMM-centered workflow and four controlled experimental interventions.
2. OrthoBench F1 and precision-recall comparison of all 12 configurations.
   Do not invent missing fusion, fission, entropy or missing-group metrics.
3. Pairwise primary-benchmark frontiers and a magnified OrthoHMM2 development
   view. Export exact three- and five-objective membership and dominators.
4. Shared-host runtime, native CPU use and native-step peak memory over 4, 8
   and 12 proteomes. Reuse the final reviewed panel, not cached campaign costs:
   `threadripper_shared_panel_snapshot_20261004_v27/panel.json`, SHA256
   `1a586d5d85262267ddfbdd4fccc43a2d43208dae5b028368964d710d2b587187`.
   Show all 27 attempts, excluded/missing values explicitly, and median/range
   summaries only for cells with three eligible repeats. No isolated-efficiency
   claim, accuracy-runtime frontier, speedup extrapolation or pooling across
   different revisions/configurations. The panel has 24 eligible observations
   and six complete cells. These are descriptive shared-host measurements,
   not timings admitted under an uncontended-host protocol.
5. Six QfO panels: three reference-pair precision-recall plots and EC, GO and
   FAS means versus their recorded assessed/eligible pair populations. The FAS
   eligible pool is not its numerical sample size. Functional panels receive
   no accuracy-frontier interpretation; omit incomparable/unavailable intervals.
6. Supplementary Figure S1: complete eight-endpoint score matrix.
7. Supplementary Figure S2: original paired and phase-adjusted OrthoBench
   uncertainty for standalone candidates and registered combination contrasts.
   Reuse saved 100000-family-resample results; do not resample again or invent
   independent-pair intervals. Preserve all original multiplicity qualifications.

The figure types follow the workflow, orthogroup, resource and QfO comparisons
in Emms et al. (2026), DOI 10.1038/s41592-026-03126-6, and Emms and Kelly
(2019), DOI 10.1186/s13059-019-1832-y. All graphics are original; no published
image or published benchmark coordinates are reused. The current QfO plots do
not reproduce unrun species-tree discordance tests from those publications.

## Delivery And Review

Deliver a concise manuscript, PDF and editable DOCX, five main figures and two
supplementary figures in PNG/PDF/SVG, full source tables, exact frontier
membership, references, executable figure code, an input/output hash manifest
and a portable ZIP. Keep older manuscripts and unrelated user files unchanged.
Test dominance/ties/missing data and evidence-to-plot ownership. Inspect rendered
figures and PDF pages for overlap, clipping, missing assets and numerical drift.
Commit and push validated code and the final package, verifying the remote head.
Submission declarations requiring human approval must remain explicitly pending;
do not invent authorship, affiliations, funding or an archival DOI.
