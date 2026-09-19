# First Corrected Reconciliation Admitted

Reconciliation array task 21760_0 (raw job 21853) completed with exit 0:0
in 1:14:36 on the shared local host. Independent admission 21761 completed
0:0 in 2:17. This is factorial p0_c0_r1: initial HMM search, no multi-sequence
profile expansion, no candidate expansion, with phylogenetic reconciliation.
It is not the publication satellite_v2 configuration (p1_c1_r1).

The native audit checks execution provenance, 157,585 recorded artifacts,
candidate ownership, partition integrity and native pair identities. It
reports 5,113,820 native ortholog pairs. All 984,137 candidate genes are
preserved across 397,041 root HOGs, from 394,328 candidate families; 692
source families are split and no cross-source merges are reported.
Group counts and pair volume are not biological accuracy or coverage of
reference ortholog relations.

The retained [admission receipt](qfo_corrected_factorial_native_admission_21761.json)
is an unchanged copy of
`benchmarks/work/qfo_corrected_factorial_native_admission_0_20260918.json`.
SHA-256: `97341e1b9ef6ac36e6b8c329aa3b5161690ccf00e617895f1d4ccb6ea127b16a`.

Pair conversion job 21766 follows this successful admission and uses native
phylogenetic pairs, not root-HOG cliques. Scoring job 21775 and independent
assessment 21776 remain downstream gates. No corrected accuracy value is
admitted by this native-output receipt. The other three reconciliation cells
remain separate serial tasks; no incomplete cell is assigned a zero score.

This validates output integrity, not independent tree correctness or
reconstruction of every orthology event. The shared-host elapsed time is
incremental reconciliation after candidate preparation, not a dedicated
end-to-end comparative runtime. No frozen prediction or default was changed.
