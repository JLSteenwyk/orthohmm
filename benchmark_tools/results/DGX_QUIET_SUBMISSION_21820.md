# Quiet-Control Submission Observation

The submission SSH command returned successfully with job21820 and exited
before the main-host UTC clock observation2026-09-18T21:45:58Z. The local
Slurm controller reported SubmitTime2026-09-18T17:45:58 and
EligibleTime2026-09-18T17:46:58 (local America/New_York time), with BeginTime
as the pending reason. The requested delay is60seconds, not an inferred
quiet period or proof that other users/services are idle.

Recipe27files were transferred and checked before submission; recipe manifest
SHA-25653a24cacbfefed66589d89ed9ae1bf4af701aa3e90453f219c2712611c8efb45.
Protocol/launchers were committed/pushed ascc3f63d and recipe verification
asbb4e8ae before this command. No DGX SSH/SCP/remote log reads will be issued
between submission completion and local-controller confirmation that all
three tasks are terminal. Main-host tests and local scheduler polls are
permitted; they do not execute a process on the DGX through SSH.

Requested partition=spark,node=spark-7ff0,exclusive20CPUs/96GiB,array0-2%1,
one-hour task limit,no requeue. Native timeout900seconds is unchanged.
Results and actual start/end scheduler observations remain separate evidence.
