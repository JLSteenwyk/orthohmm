# Corrected FastOMA Assets

Pinned assets, not executable inference. Manifest SHA-256:
`bdb6878dbab9396c6cb5360d35625ebe44149cf422b1d18dbf9b0e693a2b7975`.
The corrected species tree is explicitly null and execution unauthorized.
FastOMA must receive the admitted corrected full OrthoFinder tree; the
original-release tree is not an acceptable substitute. This remains a
supplied-tree diagnostic, not independent end-to-end tree inference.

## Verified Identities

- FastOMA workflow revision `bf6dcbaa8cf516ab6f6e074dba37eceb59a9b80e`;
  FastOMA.nf, nextflow.config and conf/base.config tracked files unchanged.
- Container ID `sha256:dd9c319b85e65a38f45f033d85663211eedccd0bfcc9c25df8048924bb4862dd`;
  registry digest `dessimozlab/fastoma@sha256:67c802a71c2d150ba8825a8c5dbc5a8a68e1dd42554e6695038043d666ac358b`.
  A read-only, network-disabled container probe reports installed FastOMA 0.3.5.
- Nextflow 22.10.8 build 5860, verified offline with Java 17.0.18. Its
  launcher, Java executable and 140-record cached capsule inventory are recorded.
- OMAmer LUCA.h5: 9,375,445,304 bytes, SHA-256
  `53adf8307e244c3901e8e432e015f16f895efb777ea1f693e739169e0defb242`.
  The resolved symlink target is recorded and rehashed after the probes.

The container labels name source revision b54f6fce5335a9e7014fc35f1aa7de93e6212eb7.
Its workflow differs from the retained host workflow: the latter omits
`--min-sequence-length` from infer-roothogs. The historical QfO native
`.command.sh` confirms that omission. Preserve the executed workflow;
do not silently replace it with a workflow inferred solely from the image
label. The historical run log also records Nextflow 22.10.8, max_cpus 180,
max_memory 700.GB and a resumed workflow. A new corrected run must start
fresh, not resume those original-release tasks.

## Validation and Remaining Work

Eight focused identity tests passed. Real input/asset verification,
container package probe and offline Nextflow probe succeeded. No inference
or output-tree staging occurred, and unrelated FastOMA checkout files
were left untouched.

After the corrected OrthoFinder tree is admitted, freeze exact fresh-copy
commands, immutable image selection, effective Nextflow configuration,
runtime/resource enforcement and native-output validation. The old shell
pair converter silently ignores malformed non-two-column rows; corrected
conversion must validate rows rather than inherit that behavior. No
historical pair result is invalidated solely by identifying that parser
weakness, and no new score is available from this asset preparation.
