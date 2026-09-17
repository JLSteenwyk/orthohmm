#!/usr/bin/env bash
# Sources are acquired separately; this recipe verifies bytes before building.
set -euo pipefail
root=$(realpath "${1:?Usage: build_dgx_external_tools.sh PROJECT_ROOT}")
test "$(uname -m)" = aarch64
src="$root/external-sources-v1"
build="$root/external-builds-v1"
prefix="$root/external-prefix-v1"
test ! -e "$build"
test ! -e "$prefix"
cd "$src"
sha256sum --check <<'CHECKSUMS'
04d14aa81962765b4d2e47a5a2ca6b97bed09ba0fac3f695c23df278616941e0  FastTree-2.1.11.c
975202a6b74c9996af871404ff043bb2152edcbda539035662514bc12d1f3431  FastTree-2.2.0.c
2876f4adc1a2de4ed206bc40896763bf208bf1a02bda52f8bfdd91cf52d73e4a  mafft-7.525.tgz
b5786897a8a8ca119eb355a5630806a4da72ea84243dba85b19a86f14757b497  mcl-14-137.tar.gz
CHECKSUMS
mkdir "$build" "$prefix"
mkdir "$prefix/fasttree-2.1.11" "$prefix/fasttree-2.2.0"
set -x
gcc -O3 -DNO_SSE -march=armv8-a -o "$prefix/fasttree-2.1.11/FastTree" "$src/FastTree-2.1.11.c" -lm
gcc -O3 -fopenmp-simd -funsafe-math-optimizations -march=armv8-a -o "$prefix/fasttree-2.2.0/FastTree" "$src/FastTree-2.2.0.c" -lm
tar --extract --gzip --file "$src/mafft-7.525.tgz" --directory "$build" --no-same-owner
make -C "$build/mafft-7.525-with-extensions/core" -j8 PREFIX="$prefix/mafft-7.525"
make -C "$build/mafft-7.525-with-extensions/core" PREFIX="$prefix/mafft-7.525" install
# RNA structural extensions are outside both protein benchmark command sets.
tar --extract --gzip --file "$src/mcl-14-137.tar.gz" --directory "$build" --no-same-owner
cd "$build/mcl-14-137"
./configure --prefix="$prefix/mcl-14-137"
make -j8
make install
"$prefix/mafft-7.525/bin/mafft" --version
"$prefix/mcl-14-137/bin/mcl" --version
"$prefix/fasttree-2.1.11/FastTree" -help
"$prefix/fasttree-2.2.0/FastTree" -help
