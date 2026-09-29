# Refinement/GC Boundary Result

The single diagnostic in the
[prespecified protocol](QFO_CPM_REFINEMENT_GC_PROTOCOL_20260928.md) completed
with exit zero. Protocol, runner and 30 passing focused tests were committed
and pushed as `7fa613a6` before launch. No retry occurred.

The child used the retained interpreter and frozen refinement functions with
the debug allocator, one-thread numerical settings, CPU affinity `[0]`,
64-GiB address-space ceiling and 300-second CPU/360-second wall limits. GC
stayed enabled at thresholds 700/10/10. The approximately 74.91-second parent
elapsed field includes child execution and post-run identity checks on the
shared host; it is not controlled inference timing.

All 12 stage markers are present in order, without additional stderr text:

| Forced generation-2 collection | Unreachable objects reported | Completed |
|---|---:|---|
| Before refinement | 33 | Yes |
| After refinement returns | 0 | Yes |
| After writing the partition | 0 | Yes |

The final original-reader readback also completed. It reports 984,137 genes,
78 species, 390,845 groups and zero selected directed refinement hits. The
23,875,927-byte partition exactly matches the retained output SHA256
`f4c6f1973bc9636828baf1e6d9be3f416a18fc8ada5302fb502c082495fc1811`.
Zero selected directed hits does not mean no graph refinement: the full
RBH graph arrays were passed to the frozen refinement function.

The [result receipt](qfo_cpm_refinement_gc_20260928.json) retains all stage
records, input/source records, command, settings, raw-output identities and
the completion result. Independent post-run checking verified 273 referenced
records covering 259 distinct paths, including stdout, stderr, the partition
and five scientific modules. The child and parent are terminal.

## Interpretation

Unlike the earlier setup-only control, this observation executes refinement,
writing and readback with actual retained inputs. It shows that all these
boundaries can complete with the stated instrumentation. It neither reproduces
the original heap history nor rules out intermittent corruption. Logging and
three forced collections deliberately perturb allocation/collection timing;
their success does not prove a GC fix, attribute the crash to delayed GC, or
justify adding forced collections to production.

Original failed admission 22155 remains failed and high-CPM accuracy remains
missing. No scientific default, score, failure status or multiplicity treatment
changed. Further work needs a reproducing native trace or a specific new
causal hypothesis; repeating this successful diagnostic is not warranted.
