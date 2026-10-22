# OrthoHMM2 Comparative Manuscript Package

Concise working manuscript drafted 20 October and revised 22 October 2026.
Six main figures cover the workflow, accuracy, resources, QfO and a new
original-versus-development runtime comparison; two supplementary figures
retain all endpoints and original uncertainty. The first five main figures
follow the panel types of the OrthoFinder papers.
Published figures and published benchmark coordinates are not reused.

## Read The Manuscript

- `manuscript.pdf`: formatted manuscript, references, figures and legends.
- `manuscript_generations_final.docx`: editable version, with embedded figure images.
- `manuscript_generations_final.html`: local browser version; keep it beside `paper.css` and assets.
- `manuscript.md`: editable source.
- `assets_v4/fig3_pareto.pdf`: vector Pareto figure.
- `assets_v4/`: the original seven figures in PDF, SVG and 220-dpi PNG, plus source data.
- `assets_speed/`: Figure 6 in all three formats, the complete nine-run timing
  snapshot, three small input proteomes, native records and pilot history.

The main text, methods and legends total approximately 3300 words, excluding
references. Jacob L. Steenwyk is the sole author; the byline contains no other
names, affiliations or contact details. Funding and competing-interest
disclosures still require confirmation. This is not an approved submission
or public archival deposition. No independent biological validation or new
default recommendation is claimed.

The current revision follows all three files in the supplied `writing_guides`
directory, including the explicit first-person and punctuation restrictions.
`validation/STYLE_REVIEW.json` records the guide fingerprints and checks.
The original ZIP remains an unchanged historical draft; the current download
is `OrthoHMM2_manuscript_20261022_speed.zip`.

OrthoHMM2 is the manuscript name for the extended development framework.
The original OrthoHMM preprint is cited with its original title and full author
list. Original run identifiers, archived scores, software version numbers and
source revisions are unchanged; this editorial naming is not a software release.

The bibliography has 17 cited entries with no unresolved keys, duplicate DOIs,
or unused entries. There are 20 reference mentions in 14 citation groups.
The initial citation revision added the 2024 QfO developments paper,
the original Schlicker similarity paper, and the FAS method paper. The original
OrthoHMM preprint retains its title and complete author metadata.
`validation/ORTHOHMM2_REVIEW.json` records naming, citation and data-preservation
checks. `validation/CITATION_AUDIT.json` records the subsequent review of all
17 DOI identities, bibliographic records and actual claim scopes. This review
corrects the 2016 QfO consortium ordering, restores four article locators lost
by citation rendering, adds finalized pagination for the OrthoFinder correction,
and clarifies the QfO adaptation of Schlicker similarity. Some original full
texts were unavailable; the report states the evidence and limits for each
reference. No biological results or benchmark scores changed.

## Scientific Interpretation

The three-objective frontier maximizes OrthoBench, VGNC and SwissTrees F1.
It contains five OrthoHMM2 configurations and full OrthoFinder among 12
complete-vector configurations. They are configurations, not six independent
tools. In the VGNC-SwissTrees F1 projection alone, full OrthoFinder is the sole
nondominated configuration. Point-estimate membership is not statistical
superiority. Adding supporting TreeFam and supplementary Three Kingdoms yields
nine nondominated configurations and illustrates objective dependence.

EC/GO/FAS conditional means are excluded from accuracy-frontier calculations.
Their populations differ, and the FAS eligible pool is not the numerical sample
size. Historical shared-host timing is descriptive; it does not supply timings
for the current weighted configurations or an accuracy-runtime frontier.

## Supplementary Tables

1. `assets_v4/scores_all_29.tsv`: exact original full-precision 29-row table,
   including repeated controls, failures and NA endpoints.
2. `assets_v4/pareto_membership.tsv`: primary, five-F1 and pairwise membership,
   with every dominator and explicit objective lists.
3. `assets_v4/timing_attempts.tsv`: all 27 resource attempts, eligibility,
   exclusions and missing pre-native measurements.
4. `assets_v4/plot_coordinates.tsv`: all 96 endpoint observations used in the
   main comparison, including precision, recall, populations and output semantics.
5. `assets_speed/runs.tsv`: all nine admitted fresh generation-timing attempts,
   with separate end-to-end inference and initial-search elapsed times.

`assets_v4/analysis.json` is the complete plotting snapshot. It preserves the
12 fixed configurations, original intervals and all resource attempts. The
figure plan was recorded before plotting but after the biological results were
known; it is a post hoc reporting plan, not a preregistration.

## Reproduce The Figures

From the extracted package, with NumPy and Matplotlib installed:

```sh
python code/build_hmm_frontier_manuscript.py \
  --snapshot assets_v4/analysis.json --output regenerated_assets
python code/plot_orthohmm_generations.py \
  --snapshot assets_speed/timing.json --output regenerated_generation_figure
```

The output must be new; existing attempts are not overwritten. This mode needs
no raw proteomes, network, GPU, external scientific executables or historical
repository runners. PDF/SVG files contain editable vector text and marks.
Backend metadata and rasterization can vary by software version; compare
`analysis.json`, exact tables and plotted coordinates, not PDF file bytes alone.

In the full repository, the `--repo` mode rechecks the admitted score/diagnostic
hashes and the complete final resource inventory before collecting the snapshot.
`REPRODUCIBILITY.json` records the render commands and software versions.
`PACKAGE_MANIFEST.json` fingerprints the delivered files. The ZIP contains only
this final package, not older figure attempts, large legacy datasets or browser
profiles. The three small generation-timing input proteomes are included.

## Provenance And Boundaries

The completed campaign was delivered at
`02426dce9922e03e00773af20a23afc371fc1125`. Its central archive is
`benchmark_tools/results/hmm_campaign_final_20261020/` in the OrthoHMM2 repository.
Current HMM experiment implementation is unchanged from scientific source
`86e39314ba462698b726f569a5d98a144c8d0086`; historical comparisons retain their
own source and conversion semantics. Some historical transitive provenance is
incomplete, as qualified in the manuscript. The plotting code does not perform
new inference, scoring, resampling, parameter selection or default promotion.

The historical resource-panel input is the final `threadripper_shared_panel_snapshot_20261004_v27`,
not an earlier incomplete snapshot. Its 27 reviewed attempts include 25 measured,
24 eligible, one excluded measured run and two pre-native aborts. Only six
method/size cells have three eligible repeats. Do not impute missing medians,
pool failed values, select fastest repeats, or infer uncontended speedups.

Figure 6 uses a separate fresh small bacterial experiment, not the historical
resource panel. Three unchanged complete proteome files contain 6432 proteins.
Original v0.1.0 uses phmmer/MCL; the frozen development checkout uses
high-sensitivity inference with and without inferred phylogeny. Nine successful
native attempts use a balanced sequential order, fresh outputs and the same
16-core affinity limit. Original nested search uses four workers with four
compute threads each; its later graph/MCL stages retain four workers. Modern
search uses 16 CPUs with four threads per worker. CPU affinity is not an
exclusive reservation. Shared-host median ratios do not assume equal outputs
or accuracy and do not establish general scalability or isolated speedups.
`assets_speed/PROTOCOL.md` and `validation/GENERATION_TIMING_REVIEW.json` record
the design, exact source bindings, successful-output checks and failed pilot.
No cached-search dependency smoke or failed/interrupted pilot contributes to
the final medians. The original seven figures and all accuracy data are unchanged.

Figures 1-5 are original alternatives to the presentation categories in
[Emms et al. 2026](https://doi.org/10.1038/s41592-026-03126-6) and
[Emms and Kelly 2019](https://doi.org/10.1186/s13059-019-1832-y).
Unrun species-tree discordance tests and unavailable seven-metric OrthoBench
diagnostics have not been invented. References were selected from the existing
reviewed project bibliography, with the original Henikoff weighting and HMMER
articles added from primary publication records.

The portable code is a copy of the tested repository builder; it is not a
standalone distribution of OrthoHMM2 or its full inference dependencies.
Manuscript publication rights and submission declarations require author approval.
