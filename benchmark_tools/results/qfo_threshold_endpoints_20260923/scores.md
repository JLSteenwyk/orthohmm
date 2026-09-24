# Corrected QfO Threshold Neighborhood

| Arm | Status | GO similarity | EC similarity | VGNC F1 | SwissTrees F1 | TreeFam-A F1 | FAS | Secondary mean |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| control | admitted | 0.490349 | 0.965650 | 0.901690 | 0.833513 | 0.614864 | 0.762993 | 0.761510 |
| cpm_low | not_admitted | NA | NA | NA | NA | NA | NA | NA |
| cpm_high | not_admitted | NA | NA | NA | NA | NA | NA | NA |
| norm_low | admitted | 0.490215 | 0.965991 | 0.901665 | 0.833513 | 0.614227 | 0.763237 | 0.761475 |
| norm_high | admitted | 0.490528 | 0.965826 | 0.901489 | 0.833513 | 0.613955 | 0.763752 | 0.761511 |
| margin_low | admitted | 0.489988 | 0.965104 | 0.902777 | 0.828391 | 0.614703 | 0.761009 | 0.760329 |
| margin_high | admitted | 0.489926 | 0.966236 | 0.900569 | 0.835648 | 0.617854 | 0.766285 | 0.762753 |

GO/EC similarity and FAS are not F1; the six-score mean is a project-defined secondary summary.

Missing CPM arms are not zero. These development-exposed results do not promote defaults.

SwissTrees uses rounded native aggregate values here; raw-family uncertainty is retained separately.

No paired uncertainty for the other five endpoints or the secondary mean is established by this table.

FAS sampling is unseeded in the retained scorer; small differences cannot be assigned to parameter changes alone.
