# Predefined CPU Pressure Controls

The [protocol](NATIVE_PRESSURE_CONTROL_PROTOCOL_20260918.md) was committed
and pushed as `b357ead` before control outcomes. Protocol SHA-256:
`812bedfc6b9127053a1de426c17f36c864005164d6928da6f396db43c46109a7`.
Implementation `df9a66b` passed 52 focused control/native/host tests before
deployment. Seven committed files, including the protocol, were exported
to a fresh DGX directory; local and remote source archives matched SHA-256
`750b3be6ed2935381d997a368f21eaa9cf6ba1520f9327e3f773df3620a169d6`.

Job21867 and all nine native steps completed0:0 on spark-7ff0. The job used
two CPUs and256MiB and completed in16 scheduler seconds. All trials retained
native CPU0 and observer CPU1; the deliberate competitor ran on native CPU0
inside the same job's batch step. No external job or service was altered.
Native and competing burns passed the fixed0.75-0.90process-CPU-second
checks. All three injected-work overlaps exceeded0.25seconds.

## Results

All values below are native-step CPU `some` cumulative stalls, not wall-time
slowdowns. The 100ms response check was defined before these observations.

| Block | Quiet (us) | Native only (us) | Contended (us) | Contended minus native only (us) | Work overlap (s) | Response check |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| 0 | 423 | 613 | 709757 | 709144 | 1.417770 | pass |
| 1 | 359 | 613 | 727125 | 726512 | 1.446071 | pass |
| 2 | 466 | 502 | 716480 | 715978 | 1.427716 | pass |

All nine injection/observation records validated. A separate local process
replayed their raw pressure fields, identities, affinities, dose/overlap
checks and the complete fixed-order summary, and verified all seven source
hashes. No trial was repeated or discarded and no threshold was changed.

The [unchanged result](native_pressure_controls_21867.json) has SHA-256
`32575bbdde2c470c218bfe69dc6d955e7b399b6adf9bdb3b9e4b355a8001dfc1`.
The [detailed scheduler record](native_pressure_controls_21867_scheduler.txt)
has SHA-256 `4e244278fd11ab024a2379540353eae6463707dccffeb2af07d29a36d8ddb1f8`.
All per-trial snapshots, worker/injector records and logs remain under
`benchmarks/work/native_pressure_controls_21867/`. Collection followed
terminal accounting for the complete job and nine steps.

## Reproduction

Extract committed sources under the recorded remote recipe directory and
submit with a fresh output directory:

```bash
sbatch --parsable --job-name=dgx_native_pressure_controls \
  --partition=spark --nodelist=spark-7ff0 --cpus-per-task=2 --mem=256M \
  --time=00:05:00 --no-requeue \
  --chdir=/home/jlsteenwyk/projects/orthohmm-publication/native_pressure_controls_df9a66b \
  --output=/home/jlsteenwyk/projects/orthohmm-publication/native_pressure_controls_%j.log \
  --wrap='/usr/bin/python3 -B benchmark_tools/run_native_pressure_controls.py --output /home/jlsteenwyk/projects/orthohmm-publication/native_pressure_controls_v1'
```

The recorded output already exists and must not be overwritten. For local
replay, load the retained JSON, run `validate_trial(trial, job_id)` for every
trial and compare its result with `trial["validation"]`; run
`summarize(report["trials"])` and compare with the saved summary. Recheck
each recorded source hash against the frozen source revision.

## Scope

This supports the native PSI response to known, sustained, same-core CPU
competition in these controls. It does not calibrate brief bursts, memory
or I/O interference, observer overhead, thermal effects, all20CPU behavior,
or a threshold for accepting inference runs. Quiet and native-only are
different workloads; their durations are not interchangeable overhead arms.
Host/native PSI must not be subtracted to identify foreign work.

The pressure response is not an OrthoHMM accuracy or efficiency result.
Scientific timing admission remains false. Older panel failures and missing
evidence remain unchanged. A separately frozen complete-command monitoring
and overhead experiment, plus justified scientific inclusion rules, is still
required before controlled scaling claims.
