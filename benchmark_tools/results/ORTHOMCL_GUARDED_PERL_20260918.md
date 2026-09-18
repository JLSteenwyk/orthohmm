# Guarded Native Perl Launch

## Policy

`run_orthomcl_perl_script.pl` removes relative and executable-hook entries
from `@INC` in a `BEGIN` block, before compiling the selected native script.
It requires an absolute script path, preserves that path as `$0`, forwards
remaining arguments, and propagates native explicit exits and loading errors.
Launchers must still supply the clean environment and pinned absolute `-I`
module directory. This is module-lookup control, not a general sandbox.

The native OrthoMCL script/module, converter arithmetic, masking settings and
parameters are unchanged. The production launcher will use this wrapper to
avoid depending on the legacy build's implicit `.` module search path.

## Real Probes

The guarded module-load probe successfully matched 116 loaded module/mapped
files to the frozen runtime and explicitly bound probe script. No relative or
hook search paths remained. The runtime snapshot is unchanged:
`9cc29e84777f47064c23096ce36e32d5465616979b05febe510219d58b6b3239`.

Perl records the script evaluated through `do` in `%INC`; the first guarded
probe rejected it as outside the original runtime snapshot. The reviewed
second probe binds exactly that helper's pre-execution hash in addition to
the runtime, rather than allowing arbitrary external modules. Failed probe
artifacts are preserved separately.

The guarded native OrthoMCL/BioPerl fixture produced the same 225-byte BPO
as the Python converter: six directed records across three queries. Both
independent index checks passed, including seven offsets with EOF sentinel.
This supports unchanged behavior for the tested HSP/gap/cutoff fixtures, not
universal semantic equivalence or a completed production inference run.

Saved reports:

- `qfo_corrected_orthomcl_guarded_perl_probe_20260918.json`, with raw outputs
  in `benchmarks/work/orthomcl_guarded_perl_runtime_probe_v2_20260918`.
- `orthomcl_guarded_bpo_parity_20260918.json`, with raw outputs in
  `benchmarks/work/orthomcl_guarded_bpo_parity_probe_20260918`.

Older unguarded reports remain historical diagnostic evidence. Their helper
hashes refer to the earlier source revisions; the new reports record current
guarded helpers. The runtime inventory and native executables were not changed.

## Verification And Remaining Work

All 68 focused launcher, runtime, BPO-content and index tests passed with
installed native probes enabled. Tests verify argument/script identity,
native exit behavior, compile/load errors, refusal to load a module from the
working directory, explicit helper binding and native BPO/index parity.

No corrected conversion or inference job is authorized by these probes.
Configured native/pair-parallel sources, conversion/index execution provenance
and final-group validation still need integration behind search admission.
Actual production module maps and outputs must be validated separately.
The reports retain false accuracy/publication flags.
