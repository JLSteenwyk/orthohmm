# OrthoBench Search Decision Diagnostic

Before diagnostic execution, freeze the driver and use scientific core
7f3a9e40dd7e79f842cc2c11fb8b548f9a802806 with the existing checksum-bound
CPU-only publication runtime. Do not rebuild libraries or alter inference.

Scope: all 81,466 distinct directed pairs from the retained 70-family trace,
covering 1,944 unique reference genes. Preserve both directions, same-species
pairs and overlapping-family labels. Reference membership determines the
diagnostic query set, not a new accuracy endpoint or tuning set.

Search reference queries against complete target proteomes, retaining input
order. Fixed settings: BLOSUM62, k=4, minimum total hits=4, minimum diagonal
hits=1, diagonal bin width=10, candidate cap=100, band width=64, full alphabet,
E-value threshold strictly below 1e-4. These match retained run metadata;
that metadata does not authenticate the historical cache executable/runtime.

Persist every unfiltered candidate index, returned score and E-value in
per-direction NPZ files, with complete query/target ID arrays. Watched-pair
TSVs distinguish prefilter exclusion, scored nonsignificance and acceptance.
Retain historical hit-presence disagreements without retrying parameters.
The final report binds input, source, runtime and output checksums. Independent
raw-output audit remains required before interpreting summaries.

Query subsets are repacked contiguously for the profile builder. A native
fixture verifies subset/full-query candidate, score and E-value equality
under the fixed cap and unchanged full target database. This bounded test
does not prove arbitrary batch invariance or historical equivalence.

The diagnostic does not score prefilter-excluded pairs counterfactually,
reconstruct past execution, establish causal accuracy effects, or supply
comparative timing. Local shared-host allocation: 4 CPUs, 16 GiB, no requeue.
The DGX timing panel is not accessed. No benchmark predictions are replaced.
