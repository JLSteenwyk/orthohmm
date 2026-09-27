# Integrated Installation, Inference, Readback And Scoring

The stdlib-only `run_integrated_publication_workflow.py` now coordinates fresh
offline installation of separate inference/reader environments, frozen native
inference, four independent scientific readers and the frozen OrthoBench
score formula. It checks the native asset and reader-export manifests, both
environment locks, supplied data digest/dimensions and file identities. It
refuses an existing output directory, records each stage, stops without retry
on failure, and kills the stage process group on timeout.

Reference files are copied separately from FASTA inputs. Inference receives
only the FASTA directory, never reference paths or scoring thresholds.
Scoring runs after successful native completion and independent readback.
There is no parameter search or production default change.

## Executed Integration Test

The [execution receipt](integrated_workflow_fixture_20260927.json) records all
eight successful stages outside the checkout: venv creation, hash-required
binary-only offline installation and pip check for each environment; native
inference; and combined independent readback/scoring. The actual fresh root
was `/tmp/orthohmm-integrated-workflow-20260927-v2`, with two CPUs.

This execution used the 16-gene, four-species installation fixture, not full
OrthoBench. It produced three root groups, 36 ortholog pairs, four duplications,
three speciations, two marker families, one reconciled family and two bypassed
families. No native checkpoint was reused. The three synthetic reference
groups follow the fixture's family labels and exercise scoring integration;
they are [test inputs](../fixtures/integration_references/family0.txt), not
independent biological truth or publishable accuracy evidence. The fixture
has no satellite merges.

Post-run independent audits find the exact installation-report inventories
(11 inference packages and five reader packages) and 5,754 installed package
payload files matching their wheels. This is a post-run audit, not a claim
that the controller audits package bytes before and after every stage.
Generated metadata/bytecode and relocated non-site payload exclusions remain.
The controller/install/native/reader file trace contains neither the original
repository prefix nor the main installation's site-packages prefix. Base Python
and OS libraries remain shared; this is not sandbox or cross-host validation.
Twenty-two focused workflow, installer and entrypoint-export tests pass.

## MAFFT Relocation Correction

Initial preflight stopped before installation or native execution. The
historical relocated asset tree contains two absolute convenience links,
`mafft/bin/mafft-distance` and `mafft/bin/mafft-profile`, resolving to the
original build prefix. Earlier inference traces did not use these links;
their bounded execution results remain valid, but the entire old asset tree
was not independently portable. Historical files, receipts and archive remain
unchanged.

The runner's `relocate_assets(source, output)` helper makes a fresh copy and
replaces only those two links with relative paths to their already copied
helpers. The new tree is `/tmp/orthohmm-integrated-assets-20260927`. Payloads
match the original pins; the validator now rejects escaping asset links.
Tests verify both relative replacements and preservation of the original
links. The initial driver/trace and subsequent successful driver/trace are
retained separately; there was one actual inference attempt.

## Full OrthoBench Entry Point

The [full-data manifest](integrated_orthobench_data_20260927.json) was prepared
from the separately acquired upstream checkout after revalidating its frozen
Git revision, blobs and input hashes: 12 FASTA files, 70 reference files and
11 low-certainty files, covering 251,378 genes. It is pinned by SHA256
`fcae062525eec61a11de876ae798acffc0f2fe9614466c8ce339a6c214666061`.
Its paths describe this acquisition; on another machine regenerate the same
file identities with local paths, then pin that new manifest's digest.

With validated local assets and a trusted installer, the complete command is:

```bash
/base/python3.10 -I -S -B /export/run_integrated_publication_workflow.py \
  --assets /relocated/native-assets --readers /exported/readers \
  --reader-wheels /local/reader-wheels --reader-lock /local/reader-requirements.txt \
  --data /local/full-data-manifest.json --data-sha256 VERIFIED_MANIFEST_SHA256 \
  --base-python /base/python3.10 --installer-python /trusted/pip/venv/bin/python \
  --output /fresh/full-workflow --cpu 32
```

The new controller's full-OrthoBench mode has **not** been executed. The
previous separately admitted full native run and full scoring restoration
remain the dataset-scale evidence; their results are not transferred to a
new controller execution. Acquisition/bootstrap, full workflow validation,
archive/deposition, data rights, controlled timing and remaining scientific
uncertainty requirements are still open.
