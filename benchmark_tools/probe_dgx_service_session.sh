#!/usr/bin/env bash
# Bounded maintenance probe only; never launches or admits benchmark timings.
set -euo pipefail

[[ ${ORTHOHMM_APPROVED_SAMWISE_STOP:-} == 20260920 ]]
[[ $(hostname) == spark-7ff0 && $(id -u) == 1000 ]]
unit=samwise-daemon-samwise.service
source_unit=/home/jlsteenwyk/.config/systemd/user/$unit
enabled_link=/home/jlsteenwyk/.config/systemd/user/default.target.wants/$unit
mask=/run/user/1000/systemd/user.control/$unit
expected=cd2e76207b56adc59dbdaba796995f59ea7dea1320abc97c8fd2ecf5015cb84a
[[ $(sha256sum "$source_unit" | cut -d ' ' -f 1) == "$expected" ]]
[[ $(readlink "$enabled_link") == "$source_unit" ]]
[[ -d /home/jlsteenwyk/Desktop/Samwise ]]
[[ ! -e "$mask" && ! -L "$mask" ]]
[[ ! -e /run/user/1000/systemd/user/$unit && ! -L /run/user/1000/systemd/user/$unit ]]
ctl() { timeout 15s systemctl --user "$@"; }
[[ $(ctl show "$unit" -p LoadState --value) == loaded ]]
prior=$(ctl show "$unit" -p ActiveState --value)
case "$prior" in
    active|activating|failed) restore_start=1 ;;
    inactive) restore_start=0 ;;
    *) printf 'Unsupported prior state: %s\n' "$prior" >&2; exit 1 ;;
esac
owned=0
stopped=0
cleanup() {
    local result=$? cleanup_failed=0
    trap - EXIT
    set +e
    if [[ $owned == 1 ]]; then
        if [[ -L "$mask" && $(readlink "$mask") == /dev/null ]]; then
            unlink "$mask" || cleanup_failed=1
        else
            printf 'Mask changed or disappeared; refusing to delete it.\n' >&2
            cleanup_failed=1
        fi
    fi
    ctl daemon-reload || cleanup_failed=1
    if [[ $(sha256sum "$source_unit" | cut -d ' ' -f 1) != "$expected" ||
          $(readlink "$enabled_link") != "$source_unit" ]]; then
        printf 'Persistent configuration changed; refusing automatic start.\n' >&2
        cleanup_failed=1
    elif [[ $stopped == 1 && $restore_start == 1 && $cleanup_failed == 0 ]]; then
        ctl start "$unit" || cleanup_failed=1
    fi
    printf '\nRESTORATION (start requested is not application health)\n'
    ctl show "$unit" -p LoadState -p ActiveState -p SubState -p MainPID -p FragmentPath || cleanup_failed=1
    sha256sum "$source_unit" || cleanup_failed=1
    [[ ! -e "$mask" && ! -L "$mask" ]] || cleanup_failed=1
    printf 'probe_exit=%s cleanup_failed=%s prior_active_state=%s\n' "$result" "$cleanup_failed" "$prior"
    if [[ $cleanup_failed != 0 ]]; then exit 1; fi
    exit "$result"
}
trap cleanup EXIT
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM

printf 'BEFORE\n'
date -u +%FT%TZ
ctl show "$unit" -p LoadState -p ActiveState -p SubState -p MainPID -p NRestarts
# Set before stop so cleanup also handles a stop observation timeout.
stopped=1
ctl stop "$unit"
mkdir -p "$(dirname "$mask")"
ln -s /dev/null "$mask"
owned=1
ctl daemon-reload
ctl reset-failed "$unit"
for sample in 0 1 2 3; do
    if [[ $sample != 0 ]]; then sleep 5; fi
    printf '\nMASKED_SAMPLE_%s\n' "$sample"
    date -u +%FT%TZ
    state=$(ctl show "$unit" -p LoadState -p ActiveState -p SubState -p MainPID -p FragmentPath)
    printf '%s\n' "$state"
    [[ $(readlink "$mask") == /dev/null ]]
    grep -qx 'LoadState=masked' <<< "$state"
    grep -qx 'ActiveState=inactive' <<< "$state"
    grep -qx 'SubState=dead' <<< "$state"
    grep -qx 'MainPID=0' <<< "$state"
    grep -Fxq "FragmentPath=$mask" <<< "$state"
done
