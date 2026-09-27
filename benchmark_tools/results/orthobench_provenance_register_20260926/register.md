# OrthoBench Score and Provenance Register

Descriptive retained evidence; not a controlled runtime ranking.

| Method | F1 (%) | Input sequence identity | Inference wall (s) | Timing scope |
|---|---:|---|---:|---|
| OrthoHMM high sensitivity | 70.3590 | Not consolidated here; retained replay evidence is separate | Unknown | No consolidated end-to-end time for this historical score row |
| OrthoHMM phylogeny satellite_v2 | 74.1061 | Not consolidated here; retained replay evidence is separate | Unknown | No consolidated end-to-end time for this historical score row |
| OrthoFinder 3.1.5 full | 72.7365 | Staged bytes and processed protein sequences match | 3,421.070 | GNU-time full native run; shared host |
| OrthoFinder 3.1.5 sequence-only checkpoint | 58.7060 | Staged bytes and processed protein sequences match | Unknown | Conversion only; sequence-only inference time unknown |
| SonicParanoid 2.0.9 | 46.7576 | Staged bytes match | 1,219.933 | Tool-reported run duration |
| ProteinOrtho 6.3.6 | 45.0573 | 177 proteins differ; 869 asterisks deleted | Unknown | No verified per-invocation duration |
| FastOMA 0.3.5 final orthologous groups | 30.9069 | IDs and residues match; renamed/reordered files | 3,023.000 | Nextflow workflow duration including failed attempts |
| OrthoMCL 1.4 | 55.0653 | 177 proteins differ; 869 asterisks deleted | 371,702.000 | July 2025 native-log interval; 8 BLAST threads |

OrthoFinder sequence-checkpoint conversion alone took 0.42 s; its inference time is unknown.
FastOMA used a supplied tree and retains 90 successful/10 exit-137 tasks. The successful final collection is scored.
OrthoMCL refers to July 2025 (8 BLAST threads), not April 2026 (32 threads).

Prediction hashes, native coverage, explicit singleton padding, available commands/settings and resource semantics are in register.json.
Proteinortho and July OrthoMCL inputs fail exact identity. Equal scores do not bound the effect of their preprocessing.
No score replacement, missing-time imputation or complete-publication-provenance claim is made.
