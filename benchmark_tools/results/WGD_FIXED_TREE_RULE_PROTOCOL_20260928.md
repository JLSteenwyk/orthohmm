# Fixed-Tree Root-Rule Diagnostic

Freeze this diagnostic before evaluating alternative-rule outcomes. This is a
post hoc, development-exposed mechanism analysis of the six previously selected
WGD examples, not independent validation, parameter optimization or default
promotion. Do not select or replace examples based on the diagnostic results.

## Fixed Inputs And Arms

- Use the seven families and all six cases in
  `biological_wgd_case_trace_20260917.json`, SHA-256
  `5b346d99d1ac063fdf169ecf76fde69390d876073e9c27d43968469936ac07cc`.
- Reproduce the retained trace and its source/admission/checkpoint checks first.
- Use the frozen `phylogeny.py` identified in that trace, SHA-256
  `216d97608e6dede8960f3b06fe7036bd23da8fae5677df88583ba9369ec3c1bf`.
- Hold candidates, rooted gene-tree bytes, species-tree bytes, leaf ownership,
  positive-paralogy pair rule and recorded satellite constraints fixed.
- Evaluate all four existing root rules in this order: `species_overlap`
  (baseline), `supported_children`, `confidence`, `mapped_event`.
- Do not reroot, realign, infer new trees, change branch support or alter
  unambiguous families that bypassed tree inference.

## Validation And Endpoints

Require exact preconstraint and final group-membership reproduction in the
baseline arm for all seven families. Recompute the original per-case separation,
homolog-support and coverage endpoints; group display IDs may differ, but
membership and all biological scoring fields must match the admitted cases.
The conflicting-reference-pillar case stays excluded from coverage/support.

For each alternative, retain preconstraint and final groups, all six scored
case rows and destinations of the five previously lost homologs. Verify that
ortholog/paralog pairs and their confidence annotations are unchanged across
arms: the intervention is only the root-group rule. Verify all input hashes
after analysis. Retain failures rather than reporting an incomplete panel as
complete. No confidence interval or significance test is planned for six cases.

## Interpretation

A changed group assignment on the identical trees demonstrates a conditional
software-rule effect. It does not prove that the recovered grouping is the
correct ancestral-copy assignment or that the topology is biologically correct.
Unchanged losses across these rules would not prove topology error: other
grouping rules, rooting choices and biological explanations remain possible.
Report gains and losses together, particularly anchor merging versus homolog
coverage. Do not recommend a global rule from these selected cases. Any later
method change requires separate evaluation under the full frozen protocols.
