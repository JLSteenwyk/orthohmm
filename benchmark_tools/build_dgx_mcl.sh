#!/usr/bin/env bash
# Recover only MCL in fresh directories; retain the failed original build.
set -euo pipefail
root=$(realpath "${1:?Usage: build_dgx_mcl.sh PROJECT_ROOT}")
test "$(uname -m)" = aarch64
test ! -e "$root/mcl-build-v3"
test ! -e "$root/mcl-prefix-v3"
cd "$root/external-sources-v1"
sha256sum --check <<'CHECKSUMS'
b5786897a8a8ca119eb355a5630806a4da72ea84243dba85b19a86f14757b497  mcl-14-137.tar.gz
cf610daf8afdedbf2110abd79bdd4121d59080cab5ec46deaf67f97273bb6bda  /usr/share/misc/config.guess
deb02c26f43b2ea64276c9ede77ec0f53d08e6256710f3c0a12275712085c348  /usr/share/misc/config.sub
CHECKSUMS
set -x
mkdir "$root/mcl-build-v3"
tar --extract --gzip --file mcl-14-137.tar.gz --directory "$root/mcl-build-v3" --no-same-owner
cd "$root/mcl-build-v3/mcl-14-137"
cp autofoo/config.guess autofoo/config.guess.original
cp autofoo/config.sub autofoo/config.sub.original
cp /usr/share/misc/config.guess autofoo/config.guess
cp /usr/share/misc/config.sub autofoo/config.sub
# Restore the common-symbol behavior expected by this 2014 C release.
CFLAGS="-g -O2 -fcommon" ./configure --prefix="$root/mcl-prefix-v3"
make -j8
make install
"$root/mcl-prefix-v3/bin/mcl" --version
