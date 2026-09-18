# FastOMA Task Validation

`audit_fastoma_tasks.py` provides a strict task-trace component for corrected
FastOMA native-output admission. It is not yet the complete admission driver.

The validator requires one successful occurrence of each fixed workflow stage:
input checking, root-HOG inference, root-HOG batching, sub-HOG collection,
native pair extraction and the native report. OMAmer task queries must cover
the exact staged proteome filenames, once each. Big/rest HOG task input names
must cover the corresponding batch output entries exactly.

Each trace hash resolves to one task directory. The task script, container
wrapper, exit file and task log are checksummed. Exit codes must be zero and
the expected native program must be present in the script. Pair extraction
must specify `--type ortholog`. Container wrappers must use the immutable
FastOMA image, network isolation and finite positive CPU/memory limits within
the configured 180-CPU/700-GiB task pool. Privileged containers are rejected.

Failed, cached, duplicated and retried trace entries are not silently accepted.
The pinned workflow itself permits retries; if any occur, they require separate
review of their attempts and retained outputs before native admission. This
strict component does not disable retry execution or discard failure evidence.

## Evidence

All 84 focused task, XML and launcher tests pass. The task tests include a
synthetic complete task chain and missing/foreign/duplicated query/batch cases.
Synthetic fixtures are not biological run evidence.

The task-directory/wrapper inspector was also applied to the actual fresh
Nextflow probe task at
`benchmarks/work/fastoma_clean_launch_probe_20260918/work/31/f6edec*`.
It confirmed exit zero, the pinned image and limits of one CPU and
268,435,456 bytes. Its checked artifact record is
`fastoma_task_inspector_probe_20260918.json`, SHA-256
`20ba46ea4d454e46e431b8f499b03c7c5408b8ca97105699c62793ce1af4a566`.
The probe does not satisfy or exercise the full biological trace requirements.

## Remaining Gates

The corrected biological run is still pending upstream input staging. Complete
admission must additionally bind successful inference accounting and execution
provenance, verify published output inventory and staged inputs, apply these
task checks, and validate native XML/table/pair contents. Script-program presence
does not authenticate every scientific argument or prove biological correctness.
Wrapper limits are not measured resource consumption or aggregate memory use.
Neither this component nor the probe admits accuracy scores.
