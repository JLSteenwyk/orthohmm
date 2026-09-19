# Root Context Control Deployment

The frozen 12-trial protocol is unchanged. This is an engineering control,
not a scientific accuracy or resource benchmark.

- Source commit: `02500e3`.
- Source archive: `benchmarks/work/root_context_controls_02500e3.tar`.
- Archive SHA-256: `f569a488f34c1aad1a55d2fc3f4852d9f8c03cadf4988fe34351e30004ba3fa5`.
- Recipe: `root_context_controls_recipe_20260919.json`.
- Recipe SHA-256: `4f98afeefda23e45b8e88d0a5a748643a8d6cbd73add44226f6dfd57ea4a93e8`.
- Verified all 508 deployed files against archive bytes and sizes; 511 total
  recipe records include directories.
- DGX root: `/home/jlsteenwyk/projects/orthohmm-publication`.
- Fresh source: `root_context_controls_recipe_v1`.
- Fresh output: `root_context_controls_v1`.
- User manager observed before submission:
  `/user.slice/user-1000.slice/user@1000.service`.

The Slurm launch script requests exclusive 20 CPUs, 96 GiB, 15 minutes,
and no requeue. All trial execution/validation failures stop further work;
unrun conditions remain explicit. No selective replacement is authorized
by this deployment. Low aggregate user-CPU response is retained as a failed
positive-control response without changing its threshold.

The scheduler queue must be checked again immediately before submission.
Retain queue and submission receipts. After submission, do not inspect or
transfer over SSH until the job is terminal; local Slurm controller queries
do not execute work on the DGX. Retrieve all raw evidence after termination,
including partial and failed outcomes. Run the new whole-panel audit against
the recipe hash and terminal scheduler record before reporting validation.

Tests before deployment: 461 focused tests passed. These are not evidence
that the deployed panel succeeds. Context overhead, actual native-tool
attribution, non-CPU isolation and scientific timing inclusion remain open.
