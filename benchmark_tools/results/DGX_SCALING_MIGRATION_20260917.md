# Dedicated DGX Scaling Migration

The user requested a different dedicated machine for timing and identified the
DGX reachable over Ethernet. Read-only connection checks succeeded using the
existing SSH alias `spark`, with strict known-host checking and batch-mode
authentication. No keys, network settings, scheduler settings or remote user
workloads were modified.

## Observed Destination

- Hostname: spark-7ff0, Ethernet10.10.10.2; local link10.10.10.1/24.
- Architecture: aarch64;20 CPU cores, comprising10 Cortex-X925 and10 Cortex-A725.
- Physical memory reported by `free`: approximately119GiB; approximately116GiB
  available at observation.
- NVIDIA GB10 GPU, no active GPU processes at observation.
- Slurm23.11.4, node spark-7ff0, partitions spark/all;20configured CPUs and
 106188MiB configured memory. Node was IDLE, CPUAlloc0, AllocMem0, load0.02.
- Existing miniforge installation under the remote account; tool-specific
  environments, compiler dependencies and destination working directory have
  not yet been established for this project.

These observations establish reachability and current capacity, not a reserved
quiet window or complete-run isolation. Slurm is shared with the existing
workstation; remote `squeue` also shows the workstation jobs. QfO21548_1 was
running on bizon when checked. Do not submit to an unconstrained default
partition and assume execution occurs on Spark.

## Prospective Timing Revision

The original32-CPU/128-GiB panel cannot run on this destination. Preserve its
unexecuted command/input manifests. Prepare a separate20-CPU/96-GiB CPU-only
panel on Spark, retaining the same27 run identities, input bytes, nested sizes,
method order, repeats and scientific settings. Record all heterogeneous CPU
affinity and actual limits. The GPU is not part of this comparison. Use fresh
native inference for allthree tools on the same destination; do not pool new
ARM timings with historical x86 measurements or label them32-CPU results.

This is a prospective infrastructure change before any scientific scaling
inference, not outcome-based tuning. New resource/command/runtime manifests
must be frozen and checked on the destination before execution. No scientific
scaling run has started there. The verified transfer and isolated installation
below supersede the original pre-transfer status.

## Architecture Gate

The frozen core revision remains7f3a9e40dd7e79f842cc2c11fb8b548f9a802806.
Its `hmm_viterbi.c` guards AVX2 kernels and contains scalar fallback paths.
However, the current native runtime builder unconditionally adds `-mavx2`,
and the frozen setup script skips that library when AVX2 is unavailable.
Neither is an ARM-ready validated build procedure. Existing `-march=native`
x86 shared libraries must not be copied as a usable runtime.

Develop a separately recorded ARM build recipe, compile the frozen sources
without unsupported x86 flags, verify allthree required libraries and symbols,
and run native profile/numerical/output checks before accepting that runtime.
Record scalar versus SIMD execution explicitly. ARM dependency availability
for OrthoFinder3.1.5, DIAMOND, MAFFT and FastTree must also be established.
Do not silently replace tools, versions or search algorithms to obtain a run.

`prepare_scaling_transfer.py` prepares byte-verified copies of the12unique
FASTA inputs and clean tracked workflow Python files, with portable relative
input paths and original nested membership. Its reference resource plan remains
the original unexecuted plan, not an authorization for destination execution.
It does not package machine-specific binaries, predictions or private keys;
nor does it claim a dependency lock or completed transfer. Native build,
destination-specific manifests, execution-wrapper integration and dedicated
workload admission remain outstanding.

## Verified Transfer

Created a new, previously absent project directory at
`/home/jlsteenwyk/projects/orthohmm-publication` on Spark. Approximately3.2TiB
was available on its filesystem. Copied the131MiB local
`benchmarks/work/scaling_transfer_v1` bundle using rsync with checksums,
strict known-host SSH verification, no deletion and no replacement of existing
files. No unrelated project or environment was modified.

All12 FASTA inputs and169 workflow files passed independent SHA256 checks both
locally and remotely. The remote bundle contains no symlinks. Its manifest SHA256
matches on both hosts:
`fcff891fc3d787df52585eb3788a3dbe6d400e8699c4f0f80ca2b66a4e6b0eb0`.
The workflow revision is2ad97c21e381a34c3df233efb8a023b27def4f0a.
Nested4/8/12-proteome inputs retain73266/165168/251378 proteins respectively.
The immutable manifest's original `prepared_not_sent` status describes its
creation; this transfer record documents the subsequent copy, not runtime
admission. Its32CPU128GiB reference plan remains unexecuted and superseded for
the prospective Spark panel.

Cloned the authorized public repository into the fresh `core` subdirectory
and checked out detached revision7f3a9e40dd7e79f842cc2c11fb8b548f9a802806.
Verified HEAD tree8138751a69846925d55f879fd5ea413b16907a3d, empty git porcelain
status and no `.so` files under the native kernel directory. No compiled x86
runtime or predictions were transferred. This establishes source/input
availability only; ARM dependencies, native builds, numerical equivalence,
target commands, resource accounting and isolation are still pending.

## ARM Build and Initial Tests

Created an isolated environment under `envs/orthohmm` using conda-forge
Python3.10.13 and pip26.0.1. Installed AArch64 wheels with exact core numerical
versions matching the original OrthoHMM environment: numpy2.2.6, numba0.65.0,
llvmlite0.47.0, igraph1.0.0, leidenalg0.11.0 and texttable1.7.0, plus
psutil7.2.2 and pytest9.0.2. `pip check` passed. Other original environments
were not changed. Transitive/system libraries are not claimed identical to
x86; a full destination lock remains required. The pip artifact report is
retained on both hosts (local benchmarks/work/dgx_orthohmm_pip_install_v1.json,
SHA b5537ed1861835eb91c524a943130c1ab369f6adccabb383800d9e3883a6d79a).

Separate builder `build_publication_arm_runtime.py` uses GCC13.3.0,
`-O3 -fopenmp -shared -fPIC -march=armv8-a` on unchanged frozen sources.
The baseline ARM ISA avoids assuming either heterogeneous core type. No CUDA
library is built. Its first attempt compiled HMM but failed an overly strict
symbol check: the multipair AVX2 entry point is correctly omitted by the C
source on ARM. The failed manifest remains
`publication_arm_runtime_failed_20260917.json` (SHA
9a65297bbb74ec50cee96b7f67334543d70a95bd799ac1e513e077bf9b8c3a2f).
The frozen engine already catches missing multipair and dispatches scalar C;
no inference source changes were needed.

Corrected builder e35157d records that symbol as optional, requires scalar HMM,
prefilter and pair-alignment entry points, and verifies `hmm_have_avx2()==0`.
A fresh frozen worktree `core_arm_v2` preserved the failed first build.
Allthree libraries compiled/loaded and the profile-construction probe passed.
Successful build manifest `publication_arm_runtime_20260917.json` has SHA
be945107f121e1a01367430ff94080bcd43922d8ef3ab9480cbc7851733778d8;
builder SHA45ad97fbf0115294935eaab6c9c3373cef1cf5785f25c5f3154ba1f0cee01b76
matches locally/remotely. Source diff remained empty after build.

On Spark, frozen test_profile.py/test_prefilter.py/test_profile_expansion.py
passed29tests in0.69s with OMP/OpenBLAS/MKL threads1. JUnit evidence retained
as benchmarks/work/dgx_arm_profile_tests_v1.xml, SHA
38d466e674cf7166ad0f98280034ebbbbd8b7ba30d00543a62e5bf597d3fa73e.
These tests and profile smoke are **not** cross-architecture score equivalence
or full pipeline validation. Numerical/decision/output comparisons, native
OrthoFinder dependencies, complete environment provenance and controlled
timing admission remain outstanding.0/27scientific scaling runs started.

## Native Score Portability Diagnostic

Probe97f9e33 ran on unchanged frozen x86/ARM runtimes with identical numerical
package versions and deterministic PCG64seed20260917. Nine queries of lengths
1/7/8/9/31/64/129/257/512 and45targets include identity, ambiguity, truncation,
insertion and random controls.405shuffled pairs were evaluated at band widths
0/1/8/64/128 using scalar C at1thread, selected native backend at4threads and
the frozen Numba reference at1thread. Native selection was multipairAVX2 on
x86 and scalarC on ARM. No benchmark labels or new method tuning were used.

| Band Width | x86 Native vs ARM Raw-Score Differences | Threshold Decisions Differ |
| --- | ---: | --- |
| 0 (unbanded) | 0/405 | No |
| 1 | 8/405 | One pair, at allthree tested thresholds |
| 8 | 3/405 | No |
| 64 (production default) | 0/405 | No |
| 128 | 0/405 | No |

Scalar C and Numba integer scores match exactly across both hosts at allfive
widths. At widths0/64/128, selected-native scores, normalized scores, E-values
and decisions also match exactly on this finite panel. Thresholds were
1e-3/1e-4/1e-5. The all-band admission gate correctly **failed**, retaining the
x86 discrepancy rather than relabeling the panel successful. This does not
establish arbitrary-sequence or end-to-end equivalence at64.

All differing pairs have both lengths<=50. Source inspection suggests a
short-pair banding inconsistency: scalar uses max(query_length,target_length)
to disable banding for short pairs, whereas multipair uses max(query_length,
maximum_target_length_in_batch). A short pair in a batch containing longer
targets can therefore remain banded. Pair96(query31,target42) scores37vs131
at width1 and changes allthree threshold decisions. This is a source-supported
mechanistic explanation, not yet an isolated patch/rescue experiment. No
frozen code or benchmark defaults were changed; preserve this as an open
correctness/release issue and test an isolated fix separately.

Machine-readable diagnostic: native_scoring_portability_20260917.json,
SHA f1cf901d14b82032955175e4b33821308f53597358f82af2c4cb677bcc9c1d18.
Raw work/native_scoring_x86_v1.json SHA
2ccab64a6f59cfba6e331a847deaabcecdd445d927e27aaa06dad1c9ef506e3d;
raw work/native_scoring_arm_v1.json SHA
ed3c20a9e8abc894f97944cc537aaf227205b2a047c8108d0161cc15ae08c62b
(both under benchmarks). Summary generated from raw records, not hand-entered
scores.22probe/summary unit tests pass. Broader default-band tests, complete
pipeline comparisons, comparator installation and timing gates remain due.

## Comparator Installation and Child-PATH Correction

Installed OrthoFinder3.1.5 in isolated `envs/orthofinder`, Python3.12.3/pip24.0,
from its official source-only release archive:
https://github.com/OrthoFinder/OrthoFinder/releases/download/v3.1.5/orthofinder-3.1.5.tar.gz
with verified upstream SHA
d292dc7ca650de6940996c980917f183366c873fe834cef81515f159e2a41bc5.
Version invocation and pip check pass. All original OrthoFinder environment
package versions match after pinning transitive dependencies, except the
irrelevant original installed orthohmm distribution is deliberately absent.
The minimal ARM environment additionally has its own packaging/build tools
and cloudpickle; inventories are retained, not claimed wholly identical.
ETE4.4.0 was built natively. A recursive comparison of installed OrthoFinder
package files excluding __pycache__ and bin found no differences between
source-only ARM installation and baseline Intel installation. This is package
source equivalence, not complete pipeline equivalence.

Built DIAMOND2.1.11 from unchanged upstream tag commit
add7f3bcb120704b8162b99d996fe9f5c110cec6 using systemGCC13.3.0 and CMakeRelease,
X86=OFF/AARCH64=ON,8parallel build workers. Build and version checks pass;
compiler emitted array-bounds warnings. Source checkout remains clean.
No search or performance claim is made from compilation alone.

**Important correction:** outer workflow PATH is not OrthoFinder's final
subprocess PATH. Newly added inspect_orthofinder_runtime.py observes the
imported, installed package's `parallel_task_manager.my_env`, executable hashes
and version commands. With the original manifest's environment overrides and
prepend paths applied, current x86 resolution is:

| Tool | Outer Workflow | OrthoFinder Child Environment |
| --- | --- | --- |
| DIAMOND | 2.1.11 | bundled2.0.13 |
| FastTree | 2.2.0 double precision | bundled2.1.11 SSE3 |
| MAFFT | 7.525 | 7.525 |
| MCL | 2.0 (old02-063installation) | bundled14-137 |

This **does not prove historical process execution**. Frozen manifests hash
both the outer entry points and bundled distribution, but their resolution
scope explicitly excludes execution tracing. Do not turn the outer-path
metadata into an unsupported claim that OrthoFinder historically used
DIAMOND2.1.11/FastTree2.2.0. Audit retained runtime evidence before assigning
historical child versions, and explicitly freeze actual child resolution for
future runs. Existing scores are not rescored or changed by this observation.

The source-only ARM installation currently resolves the newDIAMOND2.1.11 when
its build directory is supplied on PATH; MAFFT/FastTree/MCL are still missing.
That partial toolchain is **not admitted as baseline-matched**. Retain the
2.1.11build separately; obtain native2.0.13/FastTree2.1.11/MCL14-137 and verify
the full OrthoFinder dependency set before baseline-preserving timing. Shared
MAFFT7.525 and OrthoHMMFastTree2.2.0 builds also remain due. No silent upgrade
or algorithm substitution is authorized by these installations.

Evidence snapshots: orthofinder_child_resolution_x86_20260917.json SHA
f6e51c412a20d7775a1f532d91f9d430ae10971991ab19341512801efa2c5b11;
orthofinder_child_resolution_arm_20260917.json SHA
4d9e27134bb9cfb62cfb7f9421bca5c767b9e6ef499d1a039d332e55fa6cbf4e.
Local benchmarks/work retains pip-install-v1 SHA
72a80ca978ad207195394029fda6305add0c848a08cf64c0a234d60f2ae17fd6,
pip-pin-transitives-v1 SHA
1a190c830159c76e9fe3bc838df6da13596c4910f919d18fb10baa1e89afc22e
(both prefixed orthofinder), and dgx_diamond_2111_cmake_v1.txt SHA
42f3d395eaed8cffcb8dad0f097317da6a03798aab611ae28dee343800e867e9.
Three focused inspection tests pass; no benchmark inference launched.

## Companion Builds and Preserved Failures

DIAMOND2.0.13 source tag is commit2c66a0cc1d22cfd790141a94e22cb1359c40cdba.
GCC13.3/CMakeRelease/X86=OFF initially failed because upstream MemoryPool.tcc
uses uintptr_t without its standard header. Fresh build-v2 uses
`-DCMAKE_CXX_FLAGS=-include cstdint`, without modifying algorithm source.
Compilation and version checks pass. Built-in `diamond test --threads 4`
passes19/20 on ARM (twice) and on the bundled x86 baseline: XML format fails
on both; all other cases pass. Retain this common failure, not a20/20 claim.
This is not complete cross-architecture output equivalence.

Verified source downloads before building:

- MAFFT: https://mafft.cbrc.jp/alignment/software/mafft-7.525-with-extensions-src.tgz,
  SHA2876f4adc1a2de4ed206bc40896763bf208bf1a02bda52f8bfdd91cf52d73e4a.
- MCL: https://micans.org/mcl/src/mcl-14-137.tar.gz,
  SHAb5786897a8a8ca119eb355a5630806a4da72ea84243dba85b19a86f14757b497.
- FastTree2.1.11: upstream morgannprice/fasttree commit
  f76eebba0d594df96723696a1167b332e567b1cd, FastTree.c
  SHA04d14aa81962765b4d2e47a5a2ca6b97bed09ba0fac3f695c23df278616941e0.
- FastTree2.2.0: upstream tagv2.2.0 commit
  29c5e62fbcd93230ee325f9c6a17b81f00e3c72a, FastTree.c
  SHA975202a6b74c9996af871404ff043bb2152edcbda539035662514bc12d1f3431.

`build_dgx_external_tools.sh` records exact compilation/install commands and
rejects existing build/prefix directories. FastTree2.1.11 is single-precision
NoSSE3 on ARM;2.2.0 retains default double precision and is compiled with
OpenMP SIMD directives (not the multithreaded OPENMP variant). Compiler
warnings are retained. MAFFT core7.525 builds with its default multithread
support; unused RNA structural extensions are not installed. Initial script
completed both FastTree builds and MAFFT installation but failed at MCL's
2004system-detection scripts, which do not recognize aarch64.

Fresh MCLbuild-v2 with pinned system config.guess/config.sub passed configure
but failed at link time under GCC13's default no-common semantics. The final
`build_dgx_mcl.sh` builds in freshv3 directories using the same detection
scripts and `CFLAGS=-g -O2 -fcommon`, restoring the legacy common-symbol
behavior without algorithm source edits. Both failures remain on disk.
MCL14-137 installation/version check now pass. A6-node/7-edge graph smoke
atI1.2/4threads returns identical two clusters {a,b,c},{x,y,z} on x86 andARM;
this is only a small functional check.

Updated child-resolution snapshot
orthofinder_child_resolution_arm_v2_20260917.json (SHA
5f57fdcf6e81311f1c6c2ae5b8f1f4b435393c2c91a488442bf164cd5925161a)
confirms2.0.13/2.1.11/7.525/14-137 with explicit tool directories onPATH.
Inspector now includes FAMSA: **missing on ARM**. Frozen OrthoFinder defaults
to FAMSA, not MAFFT (process_args.py msa_program); bundled version reports
2.2.3-1669fc1. This exact companion and any remaining STAG/FastME dependency
checks are still required before end-to-end validation. Do not silently
replace its default alignment method with MAFFT.

Build/regression logs copied under benchmarks/work, with SHA256:

| File | SHA256 |
| --- | --- |
| diamond-2.0.13-build.log | de05768a62c02a0ecaaf456e05d4100b2820cccd2128c8e2ee1db55aa2abbc95 |
| diamond-2.0.13-build-v2.log | 3347fe3eb588d571899b13f866650d0580fad16e221b93b70d233129ab64494b |
| diamond-regression-v2.log | 1b19dfbe1c3bc2f8ee7b4cd9a051a96bcea7b51bc81b8c88d9b4f8bb8fc12e39 |
| external-tools-build-v1.log | f9e3fc7484cc81da8f913fc68f26b0b433c70146a04d35f46bd81077e2fccd6c |
| mcl-build-v2.log | cbb0459078ce7b43bacc8d2db8a613b96241c1be548c732bdb622bcdc50168e3 |
| mcl-build-v3.log | 546ab414c4d60764b8e9ac64dfd35f239709536692988abfcb56194f9a6502dc |

Both shell recipes pass bash syntax checks and runtime inspector tests3/3.
No benchmark inference or scientific timing was launched.

## FAMSA Native Build and FastME Source Gap

FAMSA now builds successfully from official refresh-bio/FAMSA commit
1669fc1444c8bc4000d71121ec2a7aa62d848b57, using GCC13.3 and
`make -j8 PLATFORM=arm8`. Its help reports 2.2.3-1669fc1 (2024-09-17),
matching the bundled x86 version. The build uses ARMv8 NEON and unchanged
algorithm sources. Generated objects/binary remain untracked in the source
checkout; this is not a claim of a clean post-build working directory.
Pinned submodules are atomic_wait 4d3bafcd562b4a47f0d5bac09fe642ec07e4eb3c,
libdeflate ced051e1c260c7008f312400a0471aecf37aa1e0, and
mimalloc 2765ec93026f445cad8f38e6b196dd226a1f6e61.

Remote binary: `famsa-source-v1/famsa`, SHA256
0a3988e9f7bae4dce39c730c0ceb3212d3ec9462ee585fe660f3861cbd38a7bb.
Build log copied to local `benchmarks/work/famsa-build-v1.log`, SHA256
108862a7a6ae41e8d9412399600506ccb561000366b003fe242409ea409540de.
Upstream compiler warnings remain in the log. Build/version success does
not establish alignment equivalence or complete-pipeline correctness.

OrthoFinder's bundled FastME reports 2.1.4. Its STAG code invokes FastME
with `-w O -s -n`; therefore FastME is now included in the runtime inspector
using `-V`. The old official download path
https://www.atgc-montpellier.fr/download/sources/fastme/fastme-2.1.4.tar.gz
returned HTTP404 on 2026-09-17. `git ls-remote --tags` of official
https://gite.lirmm.fr/atgc/FastME.git listed tags starting at v2.1.5, not
v2.1.4. This does not prove that older source is unavailable elsewhere or
in repository history. Exact-source recovery remains pending; no newer
FastME was substituted. Scientific scaling runs remain 0/27 launched.

## FAMSA Cross-Architecture Fixture Validation

`probe_famsa_portability.py` (SHA256
da62208d6bb6c3b534a228d32e867c953d3a6ab254f1b0f029020a23166c16a4)
was executed against bundled x86 FAMSA and the native ARM build. Three
hash-derived toy fixtures cover duplicated sequences, insertions/deletions
and truncation, and ambiguous residues. Each was run at1 and4 threads,
twice, yielding12 runs per host. All24 runs exited successfully, preserved
every input identifier and ungapped residue string, and produced equal
alignment lengths within each run. Canonical identifier-to-aligned-sequence
maps were identical across both hosts, both thread counts and both repeats
for each fixture. Both binary hashes remained unchanged; fixture contents
and probe script hashes matched. No accuracy or runtime endpoint is inferred
from these small fixtures; complete-pipeline portability remains unproven.

Full reports (including inputs, outputs, commands, versions and hashes):

- famsa_portability_x86_20260917.json SHA256
  4bcd6b87d94c21106e67828322e448e7c359317f60bb6119f78ef00c9647d19d.
- famsa_portability_arm_20260917.json SHA256
  dbeba47a5c432d5bba6f42519ad66f009c68ce019c076cd42f52af7086922d82.

Six focused fixture/parser tests pass, including duplicate/missing IDs,
unequal alignment lengths, changed residues and harmless record reordering.
Additional FastME source recovery checks found that the Bioconda initial
recipe and Galaxy depot begin at2.1.5, and examined GitHub mirrors retain
the same upstream2.1.5-onward history. Wayback endpoints returned503/429;
these failures do not prove that an archived2.1.4 source cannot be recovered.

## DGX Timing Collector Smoke

The prospective collector now retains absolute monotonic wrapper and
command launch/wait boundaries, plus start/finish times for each resource
observation. The recorded clock domain includes hostname and Linux boot ID;
these numbers must not be compared across hosts or boots. Host observations
retain the original observer PID and command boundaries, allowing later
interval replay without accidentally substituting the auditor's PID.
These changes do not alter old frozen measurements or their source hashes.

Slurm21625 and step21625.0 completed0:0 on spark-7ff0 with1CPU and1GiB.
The native command `/bin/sleep 5` had measured wall time5.040278321s;
six resource snapshots and four host snapshots were retained. Allocation
gates passed, all observations bracketed the command as expected, and all
three saved host intervals reproduced exactly using the recorded observer
PID. Maximum observed foreign average CPU use was0.009381cores; no large
persistent competitor was observed in this short smoke. This is not
certification of an exclusive host,20CPU/96GiB validation, a scientific
timing run, or complete resource-admission logic.

Retained report: dgx_timing_collector_smoke_20260917.json, SHA256
8e6c482aaf432169b9fd69092d85d3d8033a321abd60905ff8b82d92407222df.
Raw samples and logs are retained locally under
benchmarks/work/dgx_timing_collector_smoke_v1/timing_collector_smoke_v1
and remotely in timing_collector_smoke_v1. The report hashes those files
and collector/observer sources. Full end-to-end tool validation, prospective
DGX method/command manifests and independent timing admission remain open.
All79 focused collector, host-monitor, cgroup snapshot and resource-summary
tests pass, including native exit/failure, timeout cleanup and sampling error
retention. Existing measurements were not regenerated.

## Nested Inputs and Full Allocation Smoke

`materialize_scaling_transfer.py` checks the pinned transfer manifest, input
SHA256/size, exact inventory and nested membership before creating fresh
4/8/12-proteome directories. Every copy is rehashed. On the DGX it created
scaling_inputs_v1/{4,8,12};24 input copies preserve the original frozen
bytes. The destination manifest is retained as
dgx_scaling_inputs_20260917.json, SHA256
fa95fa60494c12d55f3bc277629c842ddbd2646a691c973882059708714e815e.
Five tests cover successful materialization, existing destination refusal,
changed hashes, extra files, symlink inputs and manifest mismatch.

Command configuration now accepts an explicit positive integer CPU count,
defaulting to32 for the original plan. A27-run comparison verifies that
20CPU changes only the allocation flags (-c or -t/-a), retaining all other
settings and run order. The original prepared manifests remain unchanged;
this does not yet establish a frozen, runnable DGX method manifest.

Slurm21626 and step21626.0 completed0:0 with an exclusive scheduler request
for20CPU/96GiB on spark-7ff0. Collector cpuset and inherited-memory-cap gates
passed. Five-second sleep wall time was5.027293531s, with six resource
snapshots and four host snapshots. All three host intervals replay exactly;
maximum observed foreign average CPU was0.033589cores. Scheduler exclusivity
does not exclude unscheduled processes, and this finite smoke does not
certify a quiet future benchmark window or measure full-load overhead.
Report: dgx_timing_allocation_smoke_20260917.json, SHA256
de5d60bfd8cb42133b6f6f1c119bcd782794f18b53145e84cd788a0aad3e7a10.

FastME2.1.4 archive was located at the alternate ATGC HTTP URL
http://atgc.lirmm.fr/download/sources/fastme/fastme-2.1.4.tar.gz.
The server reports1235934 bytes, last modified2015-10-05. HTTPS verification
failed (expired staging certificate for a different hostname); verification
was not disabled. HTTP retrieval remains unauthenticated, and no independent
historical checksum has yet been recovered. The initial30s request timed
out with117059 bytes; a resumed request is still incomplete. Partial-file
hashes and gzip failures are not final source evidence. The retrieval-only
script records a final hash only after expected length and gzip checks;
it does not execute the archive or certify authenticity.

The600s resumed request subsequently timed out (curl28), adding379089 bytes
to the initial117059. Only after this terminal result, Slurm21627 was
submitted to continue the same partial archive on bizon (1CPU/1GiB,
65-minute scheduler limit,60-minute curl limit). Retrieval script SHA256:
72f8addc2da1057abcfb29a9a1cb37e221853fc25343ae4824c1ae9e9db35bc3.
Log: benchmarks/work/fastme_source_214_21627.log. This is a download-only
job and must be independently checked before building. All1809 unit tests
pass after the input-materialization and CPU-configuration changes; the
download shell script passes bash syntax validation. No scientific scaling
run has been launched.

## OrthoHMM End-to-End Fixture

The first dataset in the frozen simulation execution order,
missing20_20261101 (8species/645proteins), was selected before running ARM
inference. Inputs were copied unchanged, with retained SHA256 lists checked
against local source files. Both methods retain4CPUs/4threads per worker,
the frozen core7f3a9e40, original scientific settings and native ARM libraries.
Phylogeny uses native MAFFT7.525/FastTree2.2.0 with inferred species tree.
This fixture is development-exposed and is not independent accuracy evidence.

Preserved setup attempts:

-21628 failed before inference: wrapper expected .fa instead of .fasta.
-21629 returned scheduler0:0, but native CLI only printed that output
  directories did not exist; neither method ran. This is a failed smoke,
  not successful inference. Fresh native directories and nonempty metrics/
  predictions are now explicit launch requirements.
-21630 completed high-sensitivity inference but phylogeny failed at the
  missing optional DendroPy dependency. High-sensitivity's group file was
  byte-identical to its x86 counterpart. Failed phylogeny metrics/log remain.
-Installed DendroPy5.0.8, matching the frozen x86 inventory, with --no-deps;
  pip check passed. Install report benchmarks/work/dgx_dendropy_install_v1.json
  SHA256 e366f34f7a1b4c5fba172c76710ffa3d2da87d15a2587119836ee500a57b0af6.
-21631 completed both fresh runs in orthohmm_pipeline_smoke_v4, scheduler
  0:0/23s. This elapsed time is not a matched performance measurement.
  Inputs, tracked core sources and all3 native libraries passed after-run
  checksum checks. Source code was not changed to repair the environment.

The retained runner's SHA256 is
b85cfe899870ddfc1ed2c18f6bbe5edded420564986a093123978ac95133e71f.
Earlier remote timing_recipe_v1/v2/v3 scripts and output directories remain
unchanged; the completed attempt used timing_recipe_v4. The final job log
is copied to benchmarks/work/dgx_orthohmm_pipeline_smoke_21631.log.

`audit_dgx_orthohmm_smoke.py` independently checks complete native metrics,
full unique input membership, declared counts, canonical unique sorted
cross-species pairs and transferred input hashes. Compared with retained
x86 predictions, all98 orthogroups agree for both methods; all98 satellite
root HOGs and1835 native ortholog pairs agree. Group labels/order do not
affect these canonical comparisons. No differences were discarded.
Report dgx_orthohmm_pipeline_smoke_20260917.json SHA256
c35b1f4228d6ff1af7d7459d2500a765604bf267f82fed54cdfe92d539eadf94.
Five focused audit tests pass. Gene/species tree and intermediate-score
equivalence are not established by this output comparison. Larger-data
portability, OrthoFinder end-to-end checks and all27 timing runs remain open.
