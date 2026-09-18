# Publication Data Rights Review

This is a partial, primary-source inventory for archive preparation, not a
legal opinion or blanket redistribution clearance. Scientific analysis and
redistribution are distinct questions. The repository MIT license does not
automatically cover acquired datasets, upstream scoring programs or binaries.

## Findings

| Material | Provider declaration / observed gap | Publication archive action |
| --- | --- | --- |
| QfO 2020.2 reference deposit | CC BY 4.0, declared for 25 deposit files | Match exact files; retain attribution, license link and modification notices. |
| UniProt QfO proteomes | Current database policy: CC BY 4.0 | Bind exact original/corrected releases and notices before bundling. |
| OrthoBench Revisited | No explicit grant located in reviewed sources | Hold new raw/scorer redistribution; clarify release-specific terms. |
| YGOB v7 | No explicit grant located in reviewed information page/README | Hold new raw sequence/pillar redistribution; seek maintainer clarification. |
| Kuzmin WGD deposited tables | CC0 declaration | Retain TableS7 identity, DOI and scientific attribution; does not cover YGOB. |
| BUSCO lineage datasets | CC BY-ND 4.0; software separately MIT | Review exact lineage and any modified reference material; prefer source acquisition instructions meanwhile. |
| Three Kingdoms input proteomes | Not cleared by this review | Trace provider/release-specific terms independently of BUSCO. |

The [QfO deposit metadata](https://zenodo.org/api/records/15087752) declares
CC BY 4.0; its scope must not be extended to missing original TreeFam inputs.
The [UniProt policy](https://www.uniprot.org/help/license), also retrieved via
the [official API](https://rest.uniprot.org/help/license), applies to
copyrightable database content and notes that other rights may exist.
The observed page revision is dated 2024-12-18; this is not itself an audit
of the notices shipped with every historical input archive.

The [OrthoBench repository](https://github.com/davidemms/Open_Orthobench)
README provides citations and benchmarking instructions. The reviewed scorer
has no explicit license declaration. Its complete recursive Git tree
`872d6f30592ab5ff837224db16a514b3f2bb916a` contains 2,133 entries and no paths
matching license/licence/copying/copyright. The downloaded README and scorer
match their tree-recorded Git blob hashes. This bounded search does not prove
that no permission exists elsewhere, nor does an article's license necessarily
settle the separate repository contents.

The [YGOB information page](http://ygob.ucd.ie/browser/info.html) and
[v7 README](http://ygob.ucd.ie/data/v7-Aug2012/README) were retrieved and
reviewed. Citation and copyright information are present, but no explicit
redistribution grant was found in these two sources. Retrieval used HTTP;
the snapshots' hashes document retained bytes, not authenticated transport.

The [WGD deposit](https://zenodo.org/records/3975054) declares CC0 through
its [metadata](https://zenodo.org/api/records/3975054). This applies to the
deposited tables, not the article PDF or independent YGOB material.

The [BUSCO license section](https://busco.ezlab.org/#license) distinguishes
software MIT from dataset CC BY-ND 4.0 and requests a release-appropriate
paper citation. The [CC BY-ND summary](https://creativecommons.org/licenses/by-nd/4.0/)
permits sharing under its conditions but restricts distributing modified
licensed material. Before packaging extracted reference content, assess its
specific scope or obtain clarification. This review does not declare every
numerical result to be an adaptation, or prohibit reporting benchmark scores.

## Retained Evidence

`publication_data_rights_20260918.json` records seven material categories,
nine source URLs, snapshot sizes/SHA-256 values, observed scope, archive
actions and unresolved questions. Snapshots total 833,426 bytes and remain
local under `benchmarks/work/publication_rights_evidence_20260918/`; no full
webpage or upstream program is newly committed. The register separates
provider declarations from project packaging decisions and unknowns.

The local figure-evidence archive has not been deposited publicly. Its
existence is not evidence that every embedded reference-derived field is
cleared for redistribution. A final file-level review must cover such fields
as well as any added raw inputs. Do not remove or silently rewrite already
retained scientific evidence to resolve a packaging question.

## Remaining Actions

1. Match intended archive entries to source releases, terms and attribution.
2. Resolve OrthoBench data/scorer and YGOB redistribution questions with
   maintainers or use approved source-acquisition-only packaging where needed.
3. Review exact Three Kingdoms proteome sources and BUSCO lineage notices.
4. Inventory software/binary/container licenses and required bundled notices
   separately, including compiled dependencies and scoring installations.
5. Review the final selected archive, rather than treating this register as
   clearance for arbitrary future contents. Retain explicit exclusions and
   reproducible acquisition commands for materials not distributed.

No email, issue, permission request or external archival submission was sent.
