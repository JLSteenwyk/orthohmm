# Base Runtime Notices Collected

The [export/relocation receipt](reconstructed_base_notices_20260927.json)
records 32 notice candidates, totaling 312,355 bytes, from the 19 exact
Conda archives used to reconstruct the base interpreter. This complements
the earlier wheel notice export; it does not replace it or clear either set
for redistribution.

Both metadata and package payload tar streams were inspected. The existing
notice-name predicate was reused, with archive SHA256 verification before
and after reading. Candidate members must be regular files with safe paths;
duplicate ZIP/tar member names are rejected. Original provider license
declarations and embedded `info/about.json` license fields remain verbatim.
Repeated notices in distinct streams are retained separately, not silently
deduplicated. Only notice candidates were exported, not executable payloads.

Local export: `benchmarks/work/publication_base_notices_20260927`.
Relocated copy: `/tmp/orthohmm-base-notices-20260927`.
Both complete inventories and all 32 file lengths/SHA256 values match.
Their identical 45,646-byte `NOTICE_INDEX.json` has SHA256
`354ae13e199dce29641d42b86b3b4c21e46ba66b891163e8dc641ed51bd7234d`.
The index retains package/archive identities, selected stream/member names,
original bytes/digests, declarations, helper-source identities and the
executed collection command.

`_libgcc_mutex` has no candidate notice and no recorded license declaration;
this gap is retained rather than assigned a guessed license. The other 18
archives contain one or more candidates. Provider declarations span several
license families and GCC runtime exceptions; their presence alone does not
resolve component attribution, corresponding-source obligations, compatibility
or the terms governing distribution of these particular package artifacts.

The candidate heuristic can miss material or include a file that is not a
license. This export does not cover operating-system libraries, external
MAFFT/FastTree notices, every statically linked component or datasets. No
public deposit, binary redistribution, runtime change or scientific score
change occurred. Full package review and the other publication gates remain.
