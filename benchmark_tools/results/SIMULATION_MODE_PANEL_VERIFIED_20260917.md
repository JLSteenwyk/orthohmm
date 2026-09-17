# Unchanged-Tree Mode Controls

All 70 tasks in array 21334 terminated successfully. Independent auditor job
21367 completed in 1:24 using executor 63b1e908a5356c68f32f0484552b3973d398c848.
The [admission report](simulation_mode_panel_verified_20260917.json), SHA256
5f3e54e6df5c4e026d03a6bea3245da939e355a0ffdf9a50509cdabf6898e82d,
contains all 140 planned method slots: 132 equivalent and eight unavailable
because the original inferred-tree baseline failed. No new non-equivalence was
observed. Availability was determined by original completion, not accuracy.

OrthoHMM supplied-own-tree controls reproduced native pairs for all 67 available
baselines; OrthoFinder did so for all 65. The auditor also checked input,
command, source, scheduler, tree and retained upstream artifact provenance.
The pilot's documented representation exceptions remain explicit. These checks
do not establish accuracy, general tree invariance or equivalence for failed
baselines.

An additional [downstream partition audit](simulation_mode_partitions_verified_20260917.json)
compared 199 partitions: 67 OrthoHMM final orthogroup partitions, 67 OrthoHMM
root-HOG partitions, and 65 OrthoFinder final orthogroup partitions. All matched
independently of group names and ordering. Unrepresented input-gene sets were
compared explicitly, without singleton imputation. Source and execution output
inventories were checksum-checked; unknown or duplicated memberships reject.
Report SHA256:
bc091032e411ae878275a5965533f1edc45ce1087a54f97dddd867765f9a779d.

## Main Tree Experiments

Frozen executor f30eb87f107002e51f1a37b2935982a3a5facb04 ran cell 0 as
21405_0. Both methods passed native admission and retained the supplied topology;
the scheduler recorded completion in 36 seconds. This is a runtime check, not
independent scientific admission or a comparable timing measurement.

After review of all controls and partition checks, submitted the remaining cells
1-209 as array 21406, maximum two concurrent tasks, four CPUs and 16 GiB per task.
Both methods remain scheduled even where the original inferred baseline failed.
Together the two arrays cover 210 tree cells and 420 planned method runs.
No accuracy scores were inspected to authorize continuation. Complete-panel
independent admission and prespecified paired uncertainty analyses remain due.
