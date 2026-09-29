# Remote CI Failed Before Test Setup

GitHub [run 36521885396](https://github.com/JLSteenwyk/orthohmm/actions/runs/36521885396)
at `81283f0b37d4edeb024f70c61e54e4b531ffaec4` is terminal with conclusion
failure. The [API observation](ci_remote_observation_20260929.json) preserves
all six job statuses and step outcomes. Docs succeeded; test-full and the
Python 3.10 fast job failed; the other three fast jobs were cancelled.

The downloaded test-full log (job 109256348324) identifies a DNS failure
resolving `files.pythonhosted.org` while pip fetched metadata for Cython 3.0.4
from the pre-existing root requirements. Pip exhausted its five displayed
attempts and exited one. This happened before `make install`, installation
of `tests/requirements.txt`, parser import or test execution. It therefore
does not establish failure or success of the new SQL parser dependency.
The other jobs' causes have not been inferred from this one log.

The raw log remains local, with size/SHA256 recorded in the observation.
Initial anonymous log download returned 403. The authenticated GitHub log
redirect required a separate request without the API Authorization header;
the resulting log was retrieved successfully. No credentials or signed log
URL are included in the retained report. The installed `gh` command is not
the GitHub Actions CLI, so inspection used GitHub's REST endpoints instead.

No workflow rerun, dependency substitution, TLS bypass or runner setting
change was made. The [local full unit run](PUBLICATION_TEST_REFRESH_20260929.md)
and separate 13-test parser follow-up remain their own evidence. They must
not be described as a passing remote CI run. Future pushes may produce new
observations, but this failed execution must remain in the provenance.

Remote test validation remains incomplete. The failure does not change
scientific scores, historical benchmark admission or the controlled-timing
requirements. No scientific inference job was launched by this inspection.
