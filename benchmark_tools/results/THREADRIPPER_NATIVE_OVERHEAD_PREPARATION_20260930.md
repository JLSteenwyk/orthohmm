# Native Collector Control Preparation

The previous goal turn pushed user-facing documentation reconciliation at
`e718ace3`. This turn prepared an engineering control for one unresolved timing
gate, without rerunning successful diagnostic 22380 or any scientific analysis.
The current scheduler queue is empty, but this is not a quiet-host certificate.
The user's unavailable quiet window remains deferred; no further host-contention
probe, scheduler submission, DGX action or process interruption was performed.

## Actual Prepared Artifact

[Prospective protocol](THREADRIPPER_NATIVE_OVERHEAD_PROTOCOL_20260930.md)
and [generated plan](threadripper_native_overhead_plan_20260930.json):

| Item | Readback |
| --- | --- |
| Complete engineering pairs | 27 |
| Planned engineering tasks | 54 |
| Methods / sizes / repeats | 3 / 3 / 3 |
| Method-size-arm cells | 18, each with three repeats |
| Parent plan SHA256 | `a9358ac4f3c2f3eb9c1d7ce6a32528f8dd4bffef25d43cf95ac9765c5f2f3d05` |
| Protocol SHA256 | `a132652a641dc451ea7aac5ce2c4d62a2b9bd5e003e5f0dd1aa047d276440727` |
| Generated plan bytes | 447,512 |
| Generated plan SHA256 | `84affc2274bf3594f679661b95b49705fa36e2577c9b12594b4643e388dd3182` |
| Checked direct parent/protocol/source/helper pins | 9, all matching |
| Execution authorized | False |
| Inference directories created | None |

The 54 engineering tasks are distinct from the unchanged 27 production runs.
They compare completed inference under periodic native-point observation versus
a boundary-only control. Both retain common whole-host monitoring. The contrast
does not isolate common monitoring cost or validate total observer overhead.
First two arm orders are counterbalanced within each method/size; the third is
fixed prospectively by a documented hash seed. Preserve all attempts and stop
after failure; no selective retry, incomplete median or overhead subtraction.

## Executed Checks

54 focused tests passed across the new planner, existing local command builder
and input preparer. They check complete identities, changed resources/methods,
repeat input drift, storage isolation, symlinks/traversal, no-overwrite and
parent/protocol checksum rejection. Planner tests do not execute inference.

A separate stdlib readback, not a call to the builder, recursively normalized
each complete run object against its parent after replacing only fresh input
and run paths and restoring the parent index. All 54 objects matched. Verified
all pairs contain both arms, all 18 cells contain three repeats, all nine
direct pins match, and every future run/input directory is absent. This compares
planned metadata, not actual runtime/input file identities or execution outcomes.

```sh
python -B -m benchmark_tools.prepare_threadripper_overhead \
  --parent benchmark_tools/results/threadripper_private_commands_20260928.json \
  --parent-sha256 a9358ac4f3c2f3eb9c1d7ce6a32528f8dd4bffef25d43cf95ac9765c5f2f3d05 \
  --protocol benchmark_tools/results/THREADRIPPER_NATIVE_OVERHEAD_PROTOCOL_20260930.md \
  --protocol-sha256 a132652a641dc451ea7aac5ce2c4d62a2b9bd5e003e5f0dd1aa047d276440727 \
  --output-root /fresh/persistent/overhead \
  --input-root /dev/shm/fresh-overhead \
  --manifest /fresh/overhead-plan.json
```

Fresh roots are required; this command only writes a plan. Historical source
pins do not become new collector/runtime/source-readiness approval. Native
boundary implementation/replay, complete independent pair auditing,
environmental handoff, stable source/runtime recipes and a quiet window remain
required. No slowdown result, engineering budget pass, production admission
or publication readiness follows. Common monitoring cost also remains unisolated.

Additional public TreeFam filename/mirror searches yielded the previously
inspected Gerstein/InterMine/CyVerse leads and method-prediction papers, not a
new original download. No archive was redownloaded and no family labels,
reference, score or uncertainty result was substituted. Do not repeat those
same searches without a new lead. No person was contacted.
