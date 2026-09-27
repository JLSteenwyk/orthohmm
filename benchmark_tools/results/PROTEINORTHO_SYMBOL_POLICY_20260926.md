# Proteinortho Symbol-Policy Fixture

The previous [input audit](OB_MATRIX_PROVENANCE_20260926.md) found that 177
historical Proteinortho OrthoBench proteins differ by deletion of 869 `*`
characters. To distinguish observed differences from assumptions about tool
behavior, the pinned Proteinortho 6.3.6 container was inspected and exercised
on three tiny, isolated indexing-only fixtures.

## Source and Behavior

The container source matches the earlier corrected-QfO runtime snapshot:
SHA-256 `b948dcb6d94059580cd16da0426bc3c1b533ffb88d86560bdfe591ef9abb8ec9`.
Line 767 defines an amino-acid alphabet without `*`. Lines 5089-5095 reject
invalid sequence symbols and suggest a separate sanitization command. This
path does not silently rewrite such sequences.

Each case uses two synthetic 80-amino-acid proteins, with one extra `*` in
the two affected cases, and `-step=1 -cpus=1 -project=fixture`. The clean,
internal-asterisk and terminal-asterisk cases use separate fresh directories.

| Input case | Exit status | Input files changed | Outcome |
|---|---:|---|---|
| Clean amino acids | 0 | No | Indexing completed |
| Internal `*` | 1 | No | Invalid-symbol error and sanitization suggestion |
| Terminal `*` | 1 | No | Invalid-symbol error and sanitization suggestion |

The [receipt](proteinortho_symbol_policy_20260926.json) includes image/runtime
and native-source hashes, source excerpts, exact commands, environment,
input hashes, log hashes and all outcomes. All three expectations pass.
Only indexing was requested; no pair search, grouping or benchmark inference
was performed. All benchmark files were left untouched.

## Interpretation

This establishes a preprocessing constraint for this pinned tool: raw `*`
characters are rejected, not silently stripped by the tested code path.
It is consistent with a need for explicit sanitization before a successful
run, but it does not establish the historical transformation's author, time,
command, or exact binary identity. The two recorded historical invocations
must not be assigned specific outcomes from this new fixture.

The failed byte-identity gate remains failed. Native-to-normalized groups still
match, but neither matched-input promotion nor a claim that accuracy is
unaffected follows. Rerunning raw inputs alone would not supply a valid
counterfactual without a declared admissible-input/preprocessing policy.
No historical score or publication default was changed.

## Reproduction

```bash
python -m benchmark_tools.probe_proteinortho_symbols \
  --image /absolute/path/to/proteinortho_6.3.6--h2b77389_0.sif \
  --output /absolute/fresh/proteinortho-symbol-policy
```

The script requires the exact image and runtime hashes and refuses to reuse
an output directory. Eleven focused symbol/difference tests pass. An initial
test expectation counted the clean fixture as 76 residues; it was corrected
to its actual 80 residues before native execution. There were no native retries.
This fixture does not complete the broader provenance or publication goal.
