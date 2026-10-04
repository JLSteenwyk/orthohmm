# Resource Manuscript Section

The [section generator](render_shared_resource_section_20261004.py) renders
the reviewed panel into prose and three tables, preserving unavailable cells,
all attempt exclusions, missing endpoints, explicit units and contention limits.
It reuses the existing full snapshot validator, recomputing counts and cell
arithmetic before rendering. It never submits inference or performs raw replay.

Eight focused [tests](../../tests/unit/test_render_shared_resource_section.py)
pass in 17.81s. Actual v24 data exercise four complete cells and five missing
summaries in each metric table. Altered counts/median and incorrect table pins
are refused. A synthetic completed-panel fixture checks six complete cells and
retained failure gaps; it is not evidence of completed native runs. Real local
build readback checks emitted prose against the same actual snapshot.

Use an externally pinned table and fresh output directory:

```sh
python -B benchmark_tools/results/render_shared_resource_section_20261004.py \
  --table benchmark_tools/results/threadripper_shared_panel_snapshot_20261004_v24/panel.json \
  --sha256 4816200c3550d80507e63f191480f09924057348494cc48d21eaa5a5c34576b2 \
  --output /absolute/fresh/resource-section
```

Outputs are `resource_section.md` and a manifest binding the table, actual
sources and prose bytes. Direct retained evidence is checked; local ignored
review artifacts are required. This is not standalone archive verification or
native reproduction. The current main manuscript remains its dated 3 October
revision; integrate the generated section and final resource figure after the
full 27-attempt panel is reviewed, then render/inspect the whole manuscript.
Do not claim the partial section is a completed study. The same generator uses
completed-panel wording only when all planned attempts are actually reviewed;
incomplete eligible-repeat cells still have no median/range.
