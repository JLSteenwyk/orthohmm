# Corrected Comparator Pair Conversion Prepared

Implemented `prepare_qfo_corrected_comparator_pairs.py` for admitted corrected
Proteinortho and SonicParanoid outputs. It requires an explicitly supplied
admission SHA-256 and the correct method-specific admission status. The
inventory must describe all 78 corrected proteomes, 984,137 input accessions
and 3,003 species pairs. Original-release admission records are rejected.

Proteinortho uses its native post-clustering graph, not search edges or group
cliques. SonicParanoid uses native species-pair tables with the existing
converter's repeated-relation removal. Converted distinct pair counts and
SonicParanoid duplicate counts must match the admitted inventory. The
unchanged native converters remain separately pinned in the recorded source
inventory. Global uniqueness follows from unique per-species-pair output and
the admitted single-species ownership of each accession, not an assumption
that arbitrary concatenated pair files are unique.

The reference mapping is selected from the frozen scoring environment
manifest, SHA-256
`e86545fd04cb644ed642cf2a31b4993fb2225ca4c05743cd03e661db129e82bc`.
Corrected validated predictions must have zero mapping losses. Any loss is
an error to investigate, not a reason to alter reference denominators.
Conversion creates a fresh method-specific directory under
`benchmarks/results/qfo_corrected_comparator_pairs_v1/`. Partial files remain
separate until count, mapping and pre/post input/source hash checks succeed.
Failure records are preserved and existing output directories are never reused.

## Execution After Admission

```bash
python benchmark_tools/prepare_qfo_corrected_comparator_pairs.py \
  --root . --method proteinortho \
  --admission benchmarks/work/qfo_corrected_proteinortho_admission_20260918.json \
  --admission-sha256 REVIEWED_ADMISSION_SHA256
```

For SonicParanoid, use method `sonic` and its independently reviewed admission
report/hash. Do not fill the placeholder from an unfinished or unreviewed
file. Participant IDs are `qfo_corrected_proteinortho` and
`qfo_corrected_sonic`. Conversion does not run inference or QfO assessment;
successful converted pairs still require separately frozen scoring and native
assessment validation. No corrected pair conversion has been launched yet.

## Validation and Ledger

26 focused tests passed across this workflow and its filter/converters.
Fixtures cover original-release rejection, wrong inventory/counts, duplicate
accounting, exact native conversion, successful stage completion, refusal to
overwrite, and preserved partial/failed output when a mapping entry is absent.
Integration fixtures isolate manifest loading; they do not replace production
checksum/provenance validation or claim a real corrected result.

Previous turn: progress, SonicParanoid admission implementation and scheduled
validation pushed as `8db0c49` and `5a932bc`. This turn prepares the next
workflow stage. At the latest scheduler check, Proteinortho 21708 remained
running at 1:16:10 and QfO assessment 21711 at 52:34; dependent admissions
21715/21712 remained pending. No live job was restarted.

Next: review terminal admission reports, bind their exact hashes, execute
conversion, independently inspect counts/hashes, then freeze corrected
assessment commands. The full publication goal remains incomplete.
