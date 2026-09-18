#!/bin/bash
#SBATCH --job-name=qfo_corrected_archive
#SBATCH --nodelist=bizon
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=08:00:00
#SBATCH --output=benchmarks/work/qfo_corrected_source_20260917/acquire_%j.log
set -euo pipefail

target=benchmarks/work/qfo_corrected_source_20260917/QfO_release_2020_04_with_updated_UP000008143.tar.gz
url=https://ftp.ebi.ac.uk/pub/databases/reference_proteomes/previous_releases/qfo_release-2020_04_with_updated_UP000008143/QfO_release_2020_04_with_updated_UP000008143.tar.gz
test -d "$(dirname "$target")"
printf 'Starting/resuming corrected archive acquisition at %s\n' "$(date -u +%FT%TZ)"
stat -c 'Existing bytes: %s' "$target"
curl --fail --location --show-error --silent --connect-timeout 30 \
  --max-time 25200 --continue-at - --header 'If-Match: "9ddf7056-5b478e2a3e88c"' \
  --output "$target" --write-out 'HTTP %{http_code}; transferred %{size_download}; seconds %{time_total}\n' "$url"
test "$(stat -c %s "$target")" -eq 2648666198
sha256sum "$target"
printf 'Download complete; archive/comparison validation still required: %s\n' "$(date -u +%FT%TZ)"
