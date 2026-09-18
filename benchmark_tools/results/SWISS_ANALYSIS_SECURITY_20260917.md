# SwissTrees Analysis Dependency Remediation

The post-push warning was verified through a read-only GitHub API snapshot:
13 open alerts (10 high, 3 medium), all for Pillow 12.2.0 in
`benchmark_tools/swiss_analysis_requirements.txt`. The alert snapshot retains
advisory IDs, reported version ranges, patched versions and source URLs, without
credentials or exploit payloads. Earlier zero-alert snapshots predate this new
analysis manifest; they are not current repository security claims.

All 13 advisories identify 12.3.0 as the first patched release. The official
[Pillow 12.3.0 release](https://github.com/python-pillow/Pillow/releases/tag/12.3.0)
was checked. Updated only the analysis pin from 12.2.0 to 12.3.0 and regenerated
the distribution-hash lock. Other ten package versions are unchanged. No frozen
inference environment, DGX recipe, scientific configuration or original result
was edited. The old reproduction report remains unchanged as historical evidence.

Hash-verified installation and `uv pip check` succeed. The clean committed-source
export at fd5f1ac was rerun with the patched analysis environment:

- Scientific JSON, including every interval and family contrast, matches exactly.
- Generated Markdown remains byte-identical.
- PDF, SVG and PNG generation succeeds; no bitwise image-equivalence claim.
- All 29 focused environment/reproduction/figure/bootstrap tests pass.

`audit_swiss_analysis_environment.py` verifies that the recorded successful
environment has exactly the expected pinned packages, then applies the existing
structured version-range evaluator to all 13 alerts. All are outside their
reported ranges. This is not comprehensive security certification or analysis
of exploitability, and it does not dismiss GitHub alerts.

Evidence:

- `dependency_alerts_analysis_env_20260917.json`: original open-alert snapshot.
- `swiss_analysis_security_audit_20260917.json`: tested environment/range audit.
- `swiss_relocated_reproduction_patched_20260917.json`: successful rerun,
  SHA256 7ea14f0bf13c257e0095f3caf4214a72e3d0325c2ccea4119335ef7cf800ca32.
- Updated lock SHA256:
  dc4b01ae88bd0d5b961d1a2027c2149e98c365a15b5e2b112f3184c986fe8863.

Use the current hash-pinned requirements when following the relocated workflow.
Do not install the historical 12.2.0 environment for routine reproduction.
Remote closure requires a separate API check after the fix is pushed.
