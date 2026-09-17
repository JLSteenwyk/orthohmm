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
run or package installation has started there. The verified transfer below
supersedes the original pre-transfer status.

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
