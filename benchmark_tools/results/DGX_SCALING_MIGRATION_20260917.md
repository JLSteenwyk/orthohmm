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
run, transfer or package installation has started there.

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
