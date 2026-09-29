# Private Packaging Alignment And YAML Linkage

Resolved the private candidate's packaging discrepancy against the retained
runtime **file hashes**, rather than accepting the older distribution label.
The public packaging 26.1 wheel contains exactly the 20 historical package
payload files, all matching sizes and SHA-256 values. Only after that check,
the private candidate's packaging 26.0 was replaced using an offline,
hash-required, dependency-disabled installation. Shared packages were untouched.

The [alignment receipt](threadripper_packaging_alignment_20260928.json) retains
the wheel, historical inventory identity, all 20 payload comparisons, commands,
logs and installed-file rechecks. `pip check` passes. All 30 installed
distribution versions match the aligned selection, and the repeated native
import report contains no shared Conda Python search/module path. The
[v2 complete hash lock](threadripper_private_requirements_v2_20260928.txt)
replaces only packaging's version/hash. The v1 installation, lock and comparison
remain historical, not descriptions of the candidate's current state.

The [new import comparison](threadripper_private_import_comparison_v2_20260928.json)
rechecks the 40 frozen scientific-core files and now matches **709/710**
historical imported scientific/site-package files. The six omitted editable
startup finders remain explicit. The only remaining file difference is
`yaml._yaml`.

## YAML Investigation

The [linkage receipt](threadripper_yaml_linkage_20260928.json) records exact
extension hashes, interpreter hashes, commands, raw `readelf -d` outputs and
runtime version/mapping observations. Historical PyYAML 6.0.1 reports libyaml
0.2.5 and loads `/home/bizon/anaconda3/lib/libyaml-0.so.2.0.9`; its extension
declares `libyaml-0.so.2` and libc dependencies with a relative RPATH.
The private wheel also reports PyYAML 6.0.1/libyaml 0.2.5, but its extension
declares libpthread/libc dependencies and maps no separate libyaml file.
These observations explain a linkage difference, not semantic or performance
equivalence. No YAML binary was copied or replaced.

The first diagnostic prepended site-packages, accidentally shadowing the
standard-library pathlib with an obsolete installed backport. It failed on
`collections.Sequence` before producing a linkage result. The corrected
diagnostic appends the explicit site-packages directory after standard-library
paths, retains `-I -S -B`, and checks the exact loaded extension path. This
probe correction did not change a benchmark, native result or installed package.

20 focused tests pass for staging, import comparison and payload alignment,
including changed/missing/extra payload rejection and no-overwrite behavior.
The private candidate is not yet bound into the timing executor. Its YAML
build difference must be explicitly covered by the deployment amendment and
native validation, alongside full runtime/linkage checks. No inference,
production timing run, cross-host restoration or timing admission is claimed.
