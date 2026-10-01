# Private QfO Deployment Admitted, Unscored

Control 22387 and corrected read-only admission 22389 are COMPLETED 0:0 on
bizon: respectively 32 CPUs/192 GiB/10:28 and 2 CPUs/64 GiB/4:04. These
shared-host elapsed times describe incremental execution/auditing, not controlled
comparative timing. Failed admission 22388 and its exact path diagnosis remain
preserved; neither native inference nor the control was retried.

The retained full admission is
`benchmarks/work/qfo_private_phylogeny_control_admission_v2_20261001.json`,
359,638 bytes, SHA256
`6865539281f9cb8c89c04b8ca16e2277056b4cc2b3293c330a12a3115879af97`.
Status `private_qfo_phylogeny_deployment_admitted_unscored`; recovered-CPM
inference authorization true, accuracy/scoring/controlled-timing/publication
flags false. Full native inventory, inferred-tree/parameter/input metadata,
complete partition and canonical native-pair validation pass. Retained 984,137
genes across 78 species, 366,068 root HOGs and 5,959,560 native ortholog pairs.

All six native exports are byte-identical to the admitted original p1_c1_r1
baseline: root HOGs, pairs, pair confidence, reconciliation events, hierarchical
groups and rooted species tree. Actual cache accounting: 351,739 candidates,
24,264 reconciled and 327,475 bypassed families; 24,264 raw-tree hits, zero
remapped hits; species-tree checkpoint reused with 26 selected families.
This is deployment parity under validated checkpoint reuse, not fresh all-tree
or search equivalence, an independent correctness proof or a new reference score.

The [independent readback](qfo_private_phylogeny_admission_readback_20261001.json),
SHA256 `713718ada456ad731df02b18f3c22219d8a76a1d7c2c45365d63f7df4ead18fc`,
uses isolated `/usr/bin/python3 -I -B` and only the standard library. It checks
both completed parents/steps, the exact report, submission/executor identity,
four source/protocol/script/test Git bindings, all six native comparisons and
1,198 unique bound records in two passes. It is file/source/completion readback,
not a third native group/pair/tree validation. Reproduce with:

```bash
/usr/bin/python3 -I -B benchmark_tools/readback_qfo_private_phylogeny_control.py --root "$PWD"
```

Initial compatibility test exposed Python 3.10 lacking hashlib.file_digest
(60 passed, one failed), while the isolated Python 3.12 readback passed. Replaced
only that standard-library call with streaming SHA256, then repeated readback
to bind the final source. All 61 focused readback/admission tests now pass.
Earlier 227 admission/original-gate/native-validator tests remain retained.

Next: separate explicit private recovered-high-CPM handoff with the already
admitted candidate/constraint bytes and unchanged frozen science. Its native
completion, output admission, lossless pair conversion, six reference metrics
and prespecified seven-arm/18-endpoint analysis remain unfinished. Do not infer
recovered-arm accuracy from baseline parity. Timing remains deferred; no DGX,
unrelated job/service action, shared package mutation or completed inference rerun.
