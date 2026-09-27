# Current Unit Regression Result

At revision `3eb446e80ec59bdef8c18d9be2cd72c33191fed6`, the complete
`tests/unit` suite passed: **10,874 passed, 10 skipped, 22 warnings**, exit
zero in 325.62 seconds. No fix or retry was required. The previous 10,819-pass
receipt remains historical; this run includes subsequent ELF inventory,
base-archive staging and GO/EC pair-panel additions.

The [machine-readable receipt](publication_unit_refresh_20260927_v3.json)
records the exact command, interpreter, pytest version, revision, JUnit hash
and all skip names/reasons. All ten skips are opt-in installed-native checks,
not successful integration tests. The 22 warnings are eleven each from
existing invalid escapes in frozen `parser.py:24` and `writer.py:57`; no
historical scientific source was edited to suppress them.

Scoped tracked source status was clean before/after for `tests/unit`,
`benchmark_tools/*.py` and `orthohmm`. Unrelated sample changes were preserved.
The completed job 22337 has separate independent admission evidence; these
tests neither rerun nor strengthen its biological scope. Test elapsed time
is not an inference benchmark. Scientific uncertainty, controlled timing,
rights review and complete release requirements remain open.
