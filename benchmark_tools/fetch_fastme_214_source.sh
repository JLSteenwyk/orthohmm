#!/usr/bin/env bash
# Retrieval only: a valid gzip is not authentication or build admission.
set -euo pipefail
destination=${1:?Usage: fetch_fastme_214_source.sh ARCHIVE_PATH}
test ! -L "$destination"
curl --fail --location --continue-at - --max-time 3600 --silent --show-error \
  http://atgc.lirmm.fr/download/sources/fastme/fastme-2.1.4.tar.gz \
  --output "$destination"
test "$(stat -c %s "$destination")" = 1235934
gzip --test "$destination"
sha256sum "$destination"
