# Four Native Cells: Conditional FAS Sampling Protocol

This is retrospective uncertainty reporting for already inspected, admitted
native QfO ablations, not an independent biological confirmation or a new
default comparison. Freeze source/tests and this plan before first numerical
application. No native scoring, RNG draws, historical retry or lookup recount.

## Actual Design And Target

Retain cells 6 P0/C0/R0, 7 P0/C0/R1, 8 P0/C1/R0, 10 P1/C0/R1. Keep
three failed/unavailable factorial cells unavailable. Initial HMM search stays
on; R-off group cliques and R-on inferred pairs retain different semantics.
Cell 7's failed timing stays null/ineligible regardless of accuracy reporting.

The pinned native scorer shuffles a precomputed finite population P and selects
k, then shuffles M annotation-eligible missing pairs and selects c=min(M,9000).
k uses its original floating-point ratio/round arithmetic, capped at P. All k
precomputed values are returned; r of the c missing selections return numeric
scores. Native raw rows are written precomputed-first, then missing. Recover
their exact sample means from that order, checking the six-decimal logged means,
counts and the full-precision aggregate; do not treat rounded logs as exact.

Condition on each retained pair set, lookup/annotations/options and fixed
numeric-return status and score in [0,1] per pair. Assume uniform sampling
without replacement within strata and independent stratum selections within
a cell. Biological scores may be arbitrarily dependent, omissions may depend
on score/features, and cells may be arbitrarily dependent. No gene-pair IID,
family independence or missing-at-random assumption is made. Uniform shuffle
and context-invariant return status remain assumptions, not historical RNG or
universal worker-behavior certification. The prior finite native context probe
supports only its tested cases.

G is the unknown number of numeric-return missing pairs. Conditional on R=r,
returned pairs are a uniform r-subset of G. Let muP and muG be the two unknown
population means, and R~Hypergeom(M,G,c). With k positive, the target is

    theta = E[Z] = w(G)*muP + (1-w(G))*muG
    w(G) = E[k/(k+R)]

Z is the actual native post-attrition raw mean, not full-eligible FAS or a
biological accuracy estimate. muP is NOT known for these new pair sets: do not
reuse the September population mean. Preserve observed Z separately; unknown
G/muP/muG preclude an exact theta point. G=0 implies w=1; empty returned
samples have no imputed mean and use the full [0,1] range.

## Simultaneous Interval Derivation

Allocate delta=.05/(4*3)=1/240 to three components per cell: inverted inclusive
hypergeometric tails for G, and two bounded finite-population sample means.
Use unchanged count_interval/mixture_weight kernels from the frozen previous
native sampling implementation. For n>0, either mean has two-sided error at
most delta within [mean-h,mean+h] intersect [0,1], with
h=sqrt(log(2/delta)/(2*n)). A known precomputed census has exact mean. For
the random returned stratum, this guarantee holds conditional on every r>0,
and [0,1] for r=0, hence also unconditionally. This follows from the ordinary
without-replacement Hoeffding bound in
[Bardenet and Maillard, Proposition 1.2](https://arxiv.org/pdf/1309.4029v2).
It does not require independent biological pair values.

w decreases with G; for fixed w the expression is linear in both means, and
for fixed means it is linear in w. Evaluate all eight corners of the two
G endpoints and the two endpoints for each mean. The resulting theta range
has failure probability at most 3*delta. Union-bound all four cell rectangles
at .05 and project ALL six pairwise differences as [L_left-U_right,
U_left-L_right]. No cross-cell independence or per-contrast selection. Preserve
favorable, unfavorable and zero-overlap results. This is simultaneous
conditional-design coverage, not family/generalization uncertainty, or an
interval for GO/EC/F1/the secondary aggregate.

## Inputs And Admission Scope

Snapshot native_qfo_scientific_scores_20261007_v3/report.json:
7aab1cbd31fb650df42e6e80a14ee0167b35a810f08202768c89c514755211a2.
Native scorer 1045c57f4d0f4787bec3d1f0690799df63c68dcc67ccfd5a925c338e3d33661d;
unchanged kernel 803ff0f455213e6a55b83e30c586ab1c57566830f8f877c916b78e84bf3658f9;
finite context result 845b4baea528f7633830e5a26ec58267da95165d8c4f79f2a0112c6063826f33.
The driver pins four original task logs, verifies current bytes, participant,
completed task identity, admitted raw/execution bindings and raw arithmetic.
Report whether logs were historically hash-bound; a newly measured log hash
must not impersonate an original admission binding. No reconstruction of the
full precomputed database or omitted pairs is attempted.

Before application, test invalid inputs, exact two-stratum subset ratios,
coverage over all selected subsets including varying precomputed samples,
score-dependent omissions, census/zero returns and corner projections. These
are implementation/design checks, not biological confirmation. Numerical
application uses SciPy 1.15.3 and preserves individual numerical failures plus
unavailable dependent contrasts at a fresh destination. Independently check
the four native count endpoints/weights and six differences after application;
no second producer run. Integrate the result and scope into existing reporting.
Other unresolved original requirements remain unresolved. No readiness claim.
