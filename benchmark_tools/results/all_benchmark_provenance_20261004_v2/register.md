# All-Tool Benchmark Provenance

Descriptive retained evidence; no comparative speed ranking.

| Dataset | Method | Declared version | Resource scopes | Input evidence |
|---|---|---|---|---|
| OrthoBench | OrthoHMM high sensitivity | 0.5.0 | Cached-hit downstream replay only; initial search excluded | Cached-hit identity verified; historical FASTA bytes not proven |
| OrthoBench | OrthoHMM phylogeny satellite_v2 | 0.5.0 | Historical full inference, shared host; scoring excluded | All 12 historical input hashes match retained FASTAs |
| OrthoBench | OrthoFinder 3.1.5 full | 3.1.5 | GNU-time full native run; shared host | Staged bytes and processed protein sequences match |
| OrthoBench | OrthoFinder 3.1.5 sequence-only checkpoint | 3.1.5 | Conversion only; sequence-only inference time unknown | Staged bytes and processed protein sequences match |
| OrthoBench | SonicParanoid 2.0.9 | 2.0.9 | Tool-reported run duration | Staged bytes match |
| OrthoBench | ProteinOrtho 6.3.6 | 6.3.6 | No verified per-invocation duration | 177 proteins differ; 869 asterisks deleted |
| OrthoBench | FastOMA 0.3.5 final orthologous groups | 0.3.5 | Nextflow workflow duration including failed attempts | IDs and residues match; renamed/reordered files |
| OrthoBench | OrthoMCL 1.4 | 1.4 | July 2025 native-log interval; 8 BLAST threads | 177 proteins differ; 869 asterisks deleted |
| QfO | OrthoHMM high sensitivity | 0.5.0 | Full inference unestablished | Corrected-release inputs; retained pins/admissions, not fresh raw-input revalidation. |
| QfO | OrthoHMM phylogeny satellite_v2 | 0.5.0 | Full inference unestablished | Corrected-release inputs; retained pins/admissions, not fresh raw-input revalidation. |
| QfO | OrthoFinder 3.1.5 full | 3.1.5 | full native inference; pair conversion | Corrected-release inputs; retained pins/admissions, not fresh raw-input revalidation. |
| QfO | OrthoFinder 3.1.5 sequence-only checkpoint | 3.1.5 | pair conversion | Corrected-release inputs; retained pins/admissions, not fresh raw-input revalidation. |
| QfO | SonicParanoid 2.0.9 | 2.0.9 | full native inference; pair conversion | Corrected-release inputs; retained pins/admissions, not fresh raw-input revalidation. |
| QfO | ProteinOrtho 6.3.6 | 6.3.6 | full native inference; pair conversion | Corrected-release inputs; retained pins/admissions, not fresh raw-input revalidation. |
| QfO | FastOMA 0.3.5 final orthologous groups | 0.3.5 | full native inference; pair conversion | Corrected-release inputs; retained pins/admissions, not fresh raw-input revalidation. |
| QfO | OrthoMCL 1.4 | 1.4 | recovered downstream mode4 only; excludes BLAST/BPO; pair conversion | Corrected-release inputs; retained pins/admissions, not fresh raw-input revalidation. |
| ThreeKingdoms | OrthoHMM high sensitivity | 0.5.0 | measured | not_retained; historical consumption not proven |
| ThreeKingdoms | OrthoHMM phylogeny satellite_v2 | 0.5.0 | measured | all_12_staged_hashes_match; historical consumption not proven |
| ThreeKingdoms | OrthoFinder 3.1.5 full | 3.1.5 | measured | all_12_staged_hashes_match; historical consumption not proven |
| ThreeKingdoms | OrthoFinder 3.1.5 sequence-only checkpoint | 3.1.5 | derived from matching full run's MCL checkpoint | all_12_staged_hashes_match; historical consumption not proven |
| ThreeKingdoms | SonicParanoid 2.0.9 | 2.0.9 | full native inference | Contemporary matched12-proteome native admission; not historical Sonic timing. |
| ThreeKingdoms | ProteinOrtho 6.3.6 | 6.3.6 | measured | not_retained; historical consumption not proven |
| ThreeKingdoms | FastOMA 0.3.5 final orthologous groups | 0.3.5 | measured | all_12_staged_hashes_match; historical consumption not proven |
| ThreeKingdoms | OrthoMCL 1.4 | 1.4 | measured stage sum | all_12_staged_hashes_match; historical consumption not proven; later native copies/all.fa/species-map content audit available |

See register.json for exact scores, prediction/admission pins, commands, CPU/RSS scope and explicit gaps.

- 24 selected rows, not a new raw-input/prediction/native/scoring admission.
- Output/input hashes copied from retained admissions are explicitly inherited, not rechecked here.
- Unknown is not zero; mixed resource scopes and allocations do not establish speed/memory rankings.
- QfO GO/EC/FAS are similarities; the six-score mean is project-defined/secondary.
- The completed replacement Threadripper timing panel is separate and is not substituted for historical score-run resources.
