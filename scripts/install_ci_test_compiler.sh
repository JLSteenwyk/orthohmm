#!/usr/bin/env bash
set -euo pipefail

: "${RUNNER_TEMP:?Require the CI temporary directory}"
: "${GITHUB_PATH:?Require the CI PATH output file}"
brew install gcc
shopt -s nullglob
compilers=("$(brew --prefix gcc)"/bin/gcc-[0-9]*)
if [[ ${#compilers[@]} != 1 || ! -x ${compilers[0]} ]]; then
    printf '%s\n' 'Require exactly one executable Homebrew GCC driver' >&2
    exit 1
fi
tools="$RUNNER_TEMP/orthohmm-ci-tools"
mkdir "$tools"
ln -s "${compilers[0]}" "$tools/gcc"
"$tools/gcc" --version
printf '%s\n' "$tools" >> "$GITHUB_PATH"
