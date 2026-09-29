# Shared Timing Runtime Drift

The planned native validation of collector v5 did **not** start. Fresh
runtime inventory regeneration stopped at its first manifest because the
entry set differed from the retained v4 binding. No v5 binding or new native
lookup receipt was written, and no fixture/production Slurm job was submitted.
The two collector source changes were the only allowed replacements; this
rejection was not bypassed by broadening that allowlist.

A subsequent read-only entry-name inventory found **549 additions and six
removals**, all under the shared `/home/bizon/anaconda3` prefix. The removed
entries are the `pyparsing-3.2.1.dist-info` directory and its files; a
`pyparsing-3.1.1.dist-info` directory is present instead. Current installed
metadata reports pyparsing 3.1.1, ProDy 2.6.1, propka 3.5.1, dm-tree 0.1.10
and ml-collections 1.1.0. This observation does not identify who changed them.

The [machine-readable diagnosis](threadripper_dependency_drift_20260928.json)
retains the complete added/removed path lists, metadata identities, prior
lookup references and per-file comparisons. All **ten pyparsing files**
listed in the retained OrthoHMM import trace now have different SHA-256
values. No pyparsing module occurs in the corresponding retained OrthoFinder
trace. Thus this is not merely a changed unused-package list. It does not,
however, establish that orthology predictions or inference speed changed.

The entry-name follow-up does not check all existing file contents. The
initial scan generated a full in-memory inventory but rejected it before
writing; those observed bytes were not retained and are not reconstructed
here. The rebinding helper now writes each observed inventory and a rejection
receipt before propagating validation failure. A regression test verifies
that failure leaves no accepted binding and preserves the changed inventory.
The expensive scan was not repeated to seek a passing result.

23 focused tests pass for rebinding, drift inspection, runtime checks and
fixture selection. The native core, scientific settings, source manifests,
shared packages and unrelated jobs were not modified. The proposed replacement
binding remains unissued; existing binding and lookup pins remain historical.

## Next Action

The [first private candidate](THREADRIPPER_PRIVATE_RUNTIME_20260928.md) now
installs successfully and restores the ten historical pyparsing import hashes.
It remains unapproved: packaging source/metadata disagree in the retained
baseline, and the PyYAML native extension differs. No shared package was reverted.

Prepare a private timing runtime with the intended dependency versions and
explicit import paths, using the existing reconstruction work where suitable.
Do not uninstall or downgrade unrelated users' packages in the shared prefix.
Any deployment amendment must preserve the frozen scientific configuration,
record all dependency differences, validate native lookup and output parity,
and pass fresh pre/post runtime checks before controlled timing. Quiet-host
eligibility, full-scale observer validation and terminal resource accounting
remain separate open requirements. This is a runtime-preflight failure, not
a failed scientific benchmark or one of the 27 production attempts.
