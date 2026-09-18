# Corrected QfO Numeric Source Admission

Job21792 completed0:0 in01:42:57 with2CPUs/192GiB. The
[retained source-equivalence report](qfo_sequence_numeric_admission_20260918.json)
has425395bytes and SHA-256
`b6e72f054051c2e5fa1c846124119f02f7819913810c66e4e94928cc8e2cb730`.
[Scheduler evidence](qfo_sequence_numeric_admission_21792_scheduler.txt)
records completion and the frozen executor command.

| Prespecified checkpoint | Genes | Directed hits | Independent admission |
|---|---:|---:|---|
| All hits | 984137 | 593510904 | Matches reconstructed source hits |
| Top100 | 984137 | 321164891 | Matches reconstructed source hits |

The numeric checker reconstructs source tuples independently and checks both
checkpoints under the frozen conversion/cap semantics. This admits source
equivalence, not graph feasibility, orthology accuracy or matched biological
sensitivity. It does not turn the top100 diagnostic into a substitute for
an unsuccessful all-hit run. Original source files/checkpoints remain in place.

Completion released graph-payload review21798 and matched Three Kingdoms
SonicParanoid inference21795. Both were confirmed running. Payload review
must finish and actual memory needs be assessed before graph submission;
no sequence-graph inference or accuracy score is claimed here.
