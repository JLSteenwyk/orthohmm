# Native Checkout Test Setup

The inspected macOS jobs confirm the scoped Bash, affinity and RSS corrections,
but still fail native alignment and load-fixture tests. This milestone prepares
an explicit OpenMP compiler and freshly built checkout libraries for CI. It
does not change kernel sources, scientific settings or historical executors.

## Remote Evidence

Source-8b48 macOS Python 3.11 job 110680012860 reports 13,584 passed, 79 failed,
20 errors, 105 skipped and 30 warnings in 489.93 seconds. Homebrew Bash is
5.3.15 and system Bash is 3.2.57. All six direct shell checks, 238 historical
handoff cases and eight relocation-fixture cases pass. This confirms the scoped
prospective runtime behavior in this job, not the cause of earlier failures or
safe execution under arbitrary shells. The 203-case affinity/boot scope gives
182 passes and 21 explicit native capability skips, matching the earlier Linux
API-absence simulation. No native Linux capability is inferred on macOS.

The same job shows eight native load-fixture errors: `/usr/bin/gcc` invokes
Clang, which rejects `-fopenmp`. Nine profile/alignment cases fail with absent
checkout `pair_align.so`; the separate baseline compilation case skips. Package
installation builds into the wheel rather than supplying checkout libraries.
Tests import the checkout, so a compiler-only or wheel-only amendment is
insufficient. Both archived failures and their
[pinned logs](ci_native_checkout_setup_20261001.json) remain retained.

Source-3e126 macOS Python 3.11 full-test job 110682852675 confirms **53 passes
and one native RSS skip** across the 54-case metrics/resource-reader scope.
Both unavailable-memory admission rejections pass. Overall unit result remains
failed: 13,601 passed, 78 failed, 20 errors, 106 skipped and 30 warnings in
546.98 seconds. This is not full-suite or native-inference confirmation.

At 02:56:54 UTC on 2 October source-8b has five terminal test failures and
wheel/docs success. At 02:59:58 UTC source-3e has full-test failure, four live
fast tests and wheel/docs success. Only the identified logs are inspected;
no sibling cause or outcome is inferred and no job is restarted.

## Prospective Setup And Validation

Only the two macOS test-job definitions gain compiler and checkout-build steps.
Install the [Homebrew GCC formula](https://formulae.brew.sh/formula/gcc), require
exactly one executable numbered driver, and expose a `gcc` symlink in a fresh
job-local tools directory. Record its version. GCC's
[OpenMP flag](https://gcc.gnu.org/onlinedocs/libgomp/Enabling-OpenMP.html) links
the corresponding runtime; no change to application compiler selection is needed.
Use the existing baseline target for this job's installation and checkout build.
Homebrew's future resolved version is not a frozen scientific-runtime pin.

The checkout helper invokes unchanged `setup.py build_py` into a fresh output.
Reject inherited libraries, incomplete/extra output, load failure, nonbaseline
Viterbi and source drift before exposing the three libraries to checkout imports.
Never overwrite existing checkout libraries. The independent Linux installed-wheel
job and documentation job remain unchanged, as do native/scientific guards.

**70 focused Linux cases pass in 10.13 seconds**, no failures/errors/skips,
including 19 new helper checks. The staged committed-package test builds all
three libraries with existing GCC 13.3.0, loads them and performs native
center-star alignment in a fresh subprocess importing that staged checkout.
Other checks cover compiler/path failures, build/load failure, scalar ABI,
unchanged sources, profile expansion and native fixture execution. The earlier
68-case run is superseded, not additive. All builds use temporary copies; no
workspace library or shared package is replaced.

Structured YAML comparison confirms exactly two added steps and one baseline
install binding per macOS test job. Existing PyYAML 6.0.1 performs this read-only
check; no package installation is needed. Setup and all three C sources remain
byte-equal to the base commit. JUnit and source pins are in the receipt.

Actual new macOS compilation, loading and native test execution remain unverified
until the next automatic CI completes. Retained-reference paths, absent raw data,
GNU utility assumptions and other failures are separate remaining work. No
scientific score/default, controlled timing, DGX, host poll/question, unrelated
workload/service or local package change occurs. The publication goal remains
active with comparable resources, uncertainty, source/rights and release gaps.
