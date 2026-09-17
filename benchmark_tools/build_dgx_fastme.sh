#!/usr/bin/env bash
# Match the bundled non-OpenMP FastME2.1.4; retrieval was unauthenticated HTTP.
set -euo pipefail
root=$(realpath "${1:?Usage: build_dgx_fastme.sh PROJECT_ROOT}")
test "$(uname -m)" = aarch64
test ! -e "$root/fastme-build-v1"
test ! -e "$root/fastme-prefix-v1"
cd "$root/external-sources-v1"
sha256sum --check <<'CHECKSUMS'
c47ae24b699c869db2d726adbe38e053179dc0c7042ccacc7a7314082ba3297c  fastme-2.1.4.tar.gz
CHECKSUMS
set -x
mkdir "$root/fastme-build-v1"
tar --extract --gzip --file fastme-2.1.4.tar.gz --directory "$root/fastme-build-v1" \
    --exclude='fastme-2.1.4/binaries' --no-same-owner
cd "$root/fastme-build-v1/fastme-2.1.4"
find src -type f -print0 | sort -z | xargs -0 sha256sum > ../source.before.sha256
gcc --version
./configure --disable-OpenMP --prefix="$root/fastme-prefix-v1"
make -j4
make install
sha256sum --check ../source.before.sha256
"$root/fastme-prefix-v1/bin/fastme" -V
sha256sum "$root/fastme-prefix-v1/bin/fastme"
