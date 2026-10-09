# Terminal Native Ablation Resource Observations

Timing measurements were collected on a shared Threadripper while other analyses were running. Competition for CPU, memory bandwidth and I/O may have affected elapsed times, with an unknown and potentially tool-dependent impact. These are observed shared-host timings, not estimates of isolated performance.

Native inference only; excludes preparation, conversion and scoring. CPU includes the task wrapper;
peak is native-step lifetime memory including its launcher, not pure algorithm RSS.

| Dataset | Cell | Job | Native Outcome | Accuracy Admitted | Resource Status | Wall (s) | CPU (s) | Peak (GiB) |
| --- | --- | --- | --- | --- | --- | ---: | ---: | ---: |
| orthobench | p0_c0_r0 | 22427 | failed_wrapper_science_recovered | True | failed_wrapper_observation_not_clean_success | 2289.7647 | 67271.6775 | 4.6140 |
| orthobench | p0_c0_r1 | 22428 | native_success | True | reviewed_native_command | 3254.9227 | 93796.6952 | 4.6299 |
| orthobench | p0_c1_r0 | 22429 | native_success | True | reviewed_native_command | 2306.6354 | 67315.2179 | 4.6278 |
| orthobench | p0_c1_r1 | 22430 | native_success | True | reviewed_native_command | 3522.0292 | 102817.7650 | 4.7245 |
| orthobench | p1_c0_r1 | 22431 | native_success | True | reviewed_native_command | 3460.6570 | 97852.1129 | 10.5875 |
| orthobench | p1_c1_r0 | 22433 | native_success | True | reviewed_native_command | 2578.5285 | 73225.2842 | 10.5704 |
| qfo_corrected | p0_c0_r0 | 22435 | native_success | True | reviewed_native_command | 38509.1572 | 1169184.7095 | 14.3655 |
| qfo_corrected | p0_c0_r1 | 22437 | native_scientific_outputs_recovered_measurement_failure_retained | True | measurement_failure_no_valid_resources | NA | NA | NA |
| qfo_corrected | p0_c1_r0 | 22444 | native_success | True | reviewed_native_command | 40676.3395 | 1226036.3493 | 18.0206 |
| qfo_corrected | p0_c1_r1 | 23891 | pre_native_cpu_binding_failure_reviewed_retained | False | pre_native_failure_no_resources | NA | NA | NA |
| qfo_corrected | p1_c0_r1 | 23902 | native_success | True | reviewed_native_command | 53995.3675 | 1585693.0692 | 18.0266 |
| qfo_corrected | p1_c1_r0 | 23985 | native_success | False | reviewed_native_command | 60568.5631 | 1607563.5080 | 19.0056 |
| qfo_corrected | p1_c1_r1 | 24036 | native_failure_retained | False | failed_native_command_not_successful_inference | 64526.1037 | 1739733.4571 | 19.0847 |

No new timing/accuracy admission, failed-attempt retry or isolated-speed ranking.
Index 7 has recovered accuracy but no valid resources; index 9 failed before native inference.
Index 11 has native observations but scoring OOM; index 12 has failed-command resources only.
Original cached factorial full costs and separate historical configuration associations remain unchanged.
