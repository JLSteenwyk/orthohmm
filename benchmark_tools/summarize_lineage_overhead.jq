{
  status: "retained_lineage_overhead_summary",
  source_audit: {path: "dgx_lineage_overhead_audit_21999_20260919.json.gz", sha256: $sha256},
  tasks: (.runs | length),
  validated_tasks,
  failed_or_unvalidated_tasks,
  complete_numeric_panel: .paired.complete_numeric_panel,
  numerical_budget_met: .paired.numerical_budget_met,
  equivalent_work_and_duration_met: .paired.equivalent_work_and_duration_met,
  methods: [.paired.methods[] | {
    method, available_pairs, expected_pairs,
    median_percent: (if .median == null then null else .median * 100 end),
    numerical_budget_met
  }],
  pairs: [.paired.pairs[] | {
    method, pair, equivalent_work, duration_gate_met, numerical_pair_budget_met,
    percent: (if .wall_ratio_minus_one == null then null else .wall_ratio_minus_one * 100 end)
  }],
  runs: [.runs[] | {
    index, method, mode, status, native_wall_s, whole_command_screen_passed,
    interval_screening_available,
    original_flags: (if .flagged_intervals == null then null else (.flagged_intervals | length) end),
    narrow_flags: (if .narrow_flagged_intervals == null then null else (.narrow_flagged_intervals | length) end)
  }],
  temporal_issues,
  observed_screens_all_pass,
  scientific_timings_admitted,
  environmental_validity_established,
  publication_ready
}
