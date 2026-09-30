# High-CPM Private Runtime: Import Failure

## Outcome

Committed/pushed controller, 38 passing focused tests and the
[fixed protocol](QFO_CPM_PRIVATE_RUNTIME_PROTOCOL_20260930.md) as
`d9b31d242bed3fa5e3dfafc471b82ccc4f8513ef`, then executed one direct bounded
attempt on the shared Threadripper. The native CLI exited 1, without timeout,
at its frozen reader's `from Bio import SeqIO` import:

```text
ModuleNotFoundError: No module named 'Bio'
```

It did not reach checkpoint loading, refinement, writing or final readback.
There is no new partition or child completion JSON. This does not reproduce,
exclude or explain the historical GC crash. The existing private environment
validated for full OrthoBench inference does not contain Biopython required by
this separate historical benchmark helper. Inference-runtime validation alone
was therefore insufficient to establish this diagnostic's import closure.

The preliminary separate process did confirm Python 3.10.13, private NumPy
2.2.6, enabled GC and thresholds 700/10/10. Its success is not success of the
refinement CLI. The controller records one CLI attempt; the failure happened
before actual refinement. No retry, package installation, forced GC, debugger,
optimizer, search or scoring followed. Original environments/results remain
unchanged. No unrelated process/service or DGX action occurred.

## Retained Evidence

Local output root: `benchmarks/work/qfo_cpm_private_runtime_20260930`.
It retains the raw controller report, runtime probe stdout/stderr, original
CLI stdout/stderr and input symlinks. Controller report is 2,199,770 bytes,
SHA256 `77e5b69b34bfc54fae8f99a8daf1c255d5490e06748786b6aa88e2c1dab0c0e4`.
The CLI stdout is empty; stderr is 853 bytes, SHA256
`b5e6dd15848f3de037fc4bc9c0d81732e830da168266c5817ccdf5cc6feb15ec`.

The [independent receipt](qfo_cpm_private_runtime_readback_20260930.json)
is 6,241 bytes, SHA256
`5bd067a1e284e6798c07a2fb4872e6a9006d5d325bda7b0195afaee3e2dda659`.
An isolated system-Python stdlib readback, without controller/scientific imports,
rechecked all 8,332 bound input/source/base/package records and all four raw logs.
It independently checks the failed status, exit/timeout, exact traceback, probe
values, missing `Bio` package directory, retained inference package inventory,
absence of output/completion files and five helper/protocol/test Git bindings.
The receipt includes exact small log texts and binds the full local report and
canonical record inventory by digest; it is not a substitute for retaining that
full inventory. All evidence rechecks passed after the failed child.

## Interpretation And Next Gate

Preserve this one failed attempt, with no automatic rerun. Any further private
runtime arm would first require a separately prepared, pinned helper-complete
environment and prospective protocol, including Biopython and other actual
imports. Do not mutate the validated OrthoBench inference environment or add
the global site packages to bypass this failure. Even success under such a new
arm would not establish the corruption cause or repair original admission 22155.

No high-CPM seed/result or downstream dependency is admitted; all accuracy and
publication flags remain false. Shared-host elapsed values are descriptive only.
The deferred controlled timing window is not needed for this correctness work;
the matched timing panel, other uncertainty and final release requirements remain
open. No score, scientific default, endpoint or completed job changed.
