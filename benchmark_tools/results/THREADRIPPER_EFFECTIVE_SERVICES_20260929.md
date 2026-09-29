# Effective Local Service Launch Evidence

The new read-only collector supplements the earlier
[unit-file inventory](THREADRIPPER_SERVICE_INVENTORY_20260928.md). It records
selected effective D-Bus properties, declared launch references and file
fingerprints before ordinary-service policy review. It does not approve any
service, construct `ordinary_processes`, change configuration or launch timing.

The [retained receipt](threadripper_effective_service_review_20260929.json)
binds three separate captures and the final source revision `f9c5d85d`.
The final capture lasted 11.59 seconds on `bizon`; boot identity was unchanged.

| Scope | Service names | Timer names | Property parse gaps | Unit rows with review gaps |
| --- | ---: | ---: | ---: | ---: |
| System | 239 | 17 | 0 | 33 |
| User | 102 | 2 | 0 | 4 |

These are loaded names, not counts of distinct programs or simultaneously
running services. Both final before/after name/state inventories agree. Reads
are sequential, so this is not continuous configuration stability.

## Retained Fields

The collector uses structured `busctl --json=short` replies rather than parsing
human-readable command lines. The systemd v255
[interface documentation](https://raw.githubusercontent.com/systemd/systemd/v255/man/org.freedesktop.systemd1.xml)
defines `ListUnits`: the numeric job ID follows the unit object path.
The local installed busctl reports systemd 255.4-1ubuntu8.6.

For all names it retains unit ID/state, fragment/drop-in paths and pending
reload status. For services it also retains working/root directories, configured
user/group, main PID/cgroup, restart setting, EnvironmentFile declarations and
all six standard launch-command arrays. Literal arguments are excluded: only
their count, serialized byte length and SHA-256 are retained. Unselected
properties, environment values, full query responses and stderr are discarded.

File observations cover 330 fragment references, 11 drop-in references,
327 absolute executable references and 28 EnvironmentFile references.
All 531 distinct successfully fingerprinted files were independently reread
after collection and still match their retained identities. Paths, hashes,
sizes and resolved locations are retained; file contents are not serialized.
Reads use 1-MiB chunks with a 128-MiB per-file envelope. Bus calls have a
cooperative 180-second budget, not a guarantee against blocked filesystem I/O.

## Explicit Review Gaps

- 27 bare executable references remain unresolved. They are retained as valid
  launch declarations, not resolved using the observer's PATH. They have no
  guessed executable fingerprint.
- 14 declared-file references remain unverified, including absent optional
  environment files, missing `rc.local` and masked unit paths resolving to
  nonregular `/dev/null`. Their role and observed failure are preserved.
- One user unit reports a pending daemon reload. No reload or repair was
  attempted; this observation does not imply approval to change the manager.
- Loaded main executable images, shared libraries, interpreted scripts,
  transient scopes and timer scheduling are not validated by this capture.
  The launch digest alone does not establish a program's purpose or safety.
- No numerical background CPU/pressure policy, process classification, native
  handoff, whole-run environment pass, causal overhead bound or quiet window
  follows from matching these file identities.

The live source queries were read-only. No service start/stop/restart/mask,
signal, scheduler configuration change, DGX operation or remote connection
was made. Do not replace unknown or changing entries with automatic exceptions
to make a production preflight pass.

## Retained Parser Corrections And Tests

The first capture failed both scope parsers before reading unit properties:
the implementation placed the job ID in the wrong struct field. It remains
at `benchmarks/work/threadripper_effective_service_inventory_20260929.json`.

The second capture reached 358 loaded names but rejected valid bare executable
names, leaving 23 launch-property parses missing; it also observed a camera
service state change. Its input and exact source snapshot remain retained.
The final `_v3.json` capture follows corrections to both parser assumptions,
retains bare names as unresolved, and has no missing selected-property parses.
Earlier captures were not overwritten or relabeled complete.

All 41 focused new/local-inventory tests pass. Coverage includes actual file
hashing without contents, nonregular/missing/oversized files, valid integer
job IDs, malformed property/command structures, bare commands without guessed
resolution, observation errors, state changes, pending reloads, cooperative
deadline refusal, local-host restriction and secret-value canaries excluded
from both successful and failed reports. These are collector tests, not service
classification or native timing admission.

All 822 source pins for pending diagnostic 22378 still match. This collector
is a new module outside that fixed source set; scientific settings and existing
collector/executor files are unchanged.

Capture future evidence only to a fresh destination:

```bash
python -B -m benchmark_tools.observe_threadripper_service_launches \
  --output /absolute/fresh/effective-services.json
```

Keep future observations separate; do not discard the retained failures,
resolve bare commands by assumption, or interpret CLI exit alone as a reviewed
background policy or controlled timing eligibility.
