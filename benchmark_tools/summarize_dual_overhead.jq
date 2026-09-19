{
  status: "retained_dual_overhead_summary",
  source_audit: {path: "dual_overhead_audit_21920.json.gz", sha256: $sha256},
  tasks: (.runs | length),
  validated_tasks,
  failed_or_unvalidated_tasks,
  complete_numeric_panel: .paired.complete_numeric_panel,
  numerical_budget_met: .paired.numerical_budget_met,
  available_pairs: ([.paired.pairs[] | select(.wall_ratio_minus_one != null)] | length),
  planned_pairs: (.paired.pairs | length),
  methods: [.paired.methods[] | {
    method, available_pairs, expected_pairs,
    median_percent: (if .median == null then null else .median * 100 end),
    numerical_budget_met
  }],
  narrow_flags_in_validated_periodic_tasks: (
    [.runs[] | select(.status == "validated" and .mode == "periodic") |
      (.narrow_flagged_intervals | length)] | add
  ),
  failed_indices: [.runs[] | select(.status != "validated") | .index],
  temporal_issues,
  scientific_timings_admitted,
  environmental_validity_established,
  publication_ready
}
