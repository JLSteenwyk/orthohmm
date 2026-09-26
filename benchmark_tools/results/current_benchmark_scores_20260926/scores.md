# Current Retained Benchmark Scores

| Method | OrthoBench | GO | EC | VGNC | SwissTrees | TreeFam-A | FAS | QfO_secondary_mean | ThreeKingdoms |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OrthoHMM high sensitivity | 0.703590 | 0.472271 | 0.932038 | 0.665933 | 0.685498 | 0.605008 | 0.776585 | 0.689555 | 0.826309 |
| OrthoHMM phylogeny satellite_v2 | 0.741061 | 0.490349 | 0.965650 | 0.901690 | 0.833513 | 0.614864 | 0.762993 | 0.761510 | 0.872133 |
| OrthoFinder 3.1.5 full | 0.727365 | 0.469548 | 0.936130 | 0.988546 | 0.848413 | 0.791918 | 0.691422 | 0.787663 | 0.988582 |
| OrthoFinder 3.1.5 sequence-only checkpoint | 0.587060 | 0.407043 | 0.753642 | 0.142755 | 0.690515 | 0.734407 | 0.561753 | 0.548352 | 0.989451 |
| SonicParanoid 2.0.9 | 0.467576 | 0.454385 | 0.872761 | 0.982794 | 0.798459 | 0.771956 | 0.736680 | 0.769506 | 0.991276 |
| ProteinOrtho 6.3.6 | 0.450573 | 0.486336 | 0.963168 | 0.954896 | 0.718111 | 0.643187 | 0.813595 | 0.763216 | 0.930463 |
| FastOMA 0.3.5 final orthologous groups | 0.309069 | 0.436894 | 0.897168 | 0.950096 | 0.780219 | 0.657902 | 0.654366 | 0.729441 | 0.861655 |
| OrthoMCL 1.4 | 0.550653 | 0.463755 | 0.921326 | 0.642027 | 0.766610 | 0.753685 | 0.733752 | 0.713526 | 0.988929 |

- All displayed values use 0-to-1 units; OrthoBench was converted from percent.
- OrthoBench is weighted orthogroup F1; Three Kingdoms is BUSCO-reference pair F1.
- QfO GO/EC/FAS are similarities, not F1. Its six-endpoint mean is project-defined and secondary.
- No cross-dataset mean or universal ranking is defined; prediction semantics differ.
- Three Kingdoms input-consumption gaps remain; only SonicParanoid uses the contemporary matched-input row.
- FastOMA uses supplied-tree configurations; the OrthoFinder sequence checkpoint is diagnostic.
- Five OrthoBench rows lack direct prediction-file hashes in this source report; full transitive provenance is not consolidated here.
- Development-exposed evidence; no inference, raw scoring, uncertainty or resource comparison rerun.
