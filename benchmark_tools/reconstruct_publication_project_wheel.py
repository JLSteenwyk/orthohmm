"""Reconstruct the exact historical project wheel wrapper from verified payloads."""

import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys
import zipfile
import zlib

NAME = "orthohmm-0.5.0-cp310-cp310-linux_x86_64.whl"
EXPECTED = dict(bytes=144444, sha256="cfdfde5ed1be29e4080dd3571f5c0fc5fc5f45b57ebe5ee9b3c7559096c3b93d")
# Only ZIP metadata and payload identities from the hash-pinned historical wheel.
MEMBERS = [
  [
    "orthohmm/__init__.py",
    [
      2026,
      9,
      26,
      20,
      58,
      32
    ],
    2176057344,
    0,
    "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  ],
  [
    "orthohmm/__main__.py",
    [
      2026,
      9,
      26,
      20,
      58,
      32
    ],
    2176057344,
    127,
    "59690e02073880a5ec60087a93d77e3448a3370bb78fb17ee9ce9d5300d0d45a"
  ],
  [
    "orthohmm/accuracy.py",
    [
      2026,
      9,
      26,
      20,
      58,
      32
    ],
    2176057344,
    14046,
    "1a35944ab7fea859143f599b1262aa111272787f6f735b9acc37d299427f2ad6"
  ],
  [
    "orthohmm/args_processing.py",
    [
      2026,
      9,
      26,
      20,
      58,
      32
    ],
    2176057344,
    7908,
    "12648aeaa4e8d0b0af60f707e90ec2c95b6983c4498dacd9bccfade1839f4c54"
  ],
  [
    "orthohmm/externals.py",
    [
      2026,
      9,
      26,
      20,
      58,
      32
    ],
    2176057344,
    10166,
    "5304b4190b01a53c7c4e3235f9d1df9e7eff1f01d09633d07750103487268115"
  ],
  [
    "orthohmm/files.py",
    [
      2026,
      9,
      26,
      20,
      58,
      32
    ],
    2176057344,
    2977,
    "f1ec2acaea2566615b06b0a37e17b77f76bbd87fe6f2c76560a318343f6871a1"
  ],
  [
    "orthohmm/helpers.py",
    [
      2026,
      9,
      26,
      20,
      58,
      32
    ],
    2176057344,
    32293,
    "faf9e85c6af8ad5f5740c81157e8ce40bbd98c1731991d9f3fbe1f47a05f6149"
  ],
  [
    "orthohmm/leiden_worker.py",
    [
      2026,
      9,
      26,
      20,
      58,
      32
    ],
    2176057344,
    1411,
    "7d6d9fa5827bd99dcce28d2248688bbeb762b4b547bbc4e7491c16f0e2090c41"
  ],
  [
    "orthohmm/metrics.py",
    [
      2026,
      9,
      26,
      20,
      58,
      32
    ],
    2176057344,
    6812,
    "d19bc9c9777b2dd01943f6421930f9fab02b7153519ab76e6a40272328195a1a"
  ],
  [
    "orthohmm/orthohmm.py",
    [
      2026,
      9,
      26,
      20,
      58,
      32
    ],
    2176057344,
    30663,
    "2afb89b9dc683e64e58208188f720e07701d4760ff09d1ac3a53c7a8075b84bb"
  ],
  [
    "orthohmm/parser.py",
    [
      2026,
      9,
      26,
      20,
      58,
      32
    ],
    2176057344,
    15939,
    "9a959b9948d944fc5361a661ac3faad55575674993c746800daf6db7eccceba5"
  ],
  [
    "orthohmm/phylogeny.py",
    [
      2026,
      9,
      26,
      20,
      58,
      32
    ],
    2176057344,
    29732,
    "216d97608e6dede8960f3b06fe7036bd23da8fae5677df88583ba9369ec3c1bf"
  ],
  [
    "orthohmm/phylogeny_pipeline.py",
    [
      2026,
      9,
      26,
      20,
      58,
      32
    ],
    2176057344,
    41369,
    "44e00316546b5df78354badd2b1a6bb595b685e98b86f686f1a32543d7f15e4f"
  ],
  [
    "orthohmm/refinement.py",
    [
      2026,
      9,
      26,
      20,
      58,
      32
    ],
    2176057344,
    61183,
    "991f1eb6a5f73d0442529ed19095a34b7c6ba8bff8dfe43a1127e24ec73fb26d"
  ],
  [
    "orthohmm/version.py",
    [
      2026,
      9,
      26,
      20,
      58,
      32
    ],
    2176057344,
    22,
    "2c12b8ea17aeb6f9f72a6b02aca218bbc4508a46df9e365a8f6c45ad769e0b34"
  ],
  [
    "orthohmm/writer.py",
    [
      2026,
      9,
      26,
      20,
      58,
      32
    ],
    2176057344,
    2949,
    "d3dc7a318166272762395dfcb772d8779286bbf60023b2934a7d7bef0e0003fe"
  ],
  [
    "orthohmm/search/__init__.py",
    [
      2026,
      9,
      26,
      20,
      58,
      32
    ],
    2176057344,
    0,
    "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  ],
  [
    "orthohmm/search/engine.py",
    [
      2026,
      9,
      26,
      20,
      58,
      32
    ],
    2176057344,
    21146,
    "596e6e56b6dd3c2e72bbfa5bb67087285bd22833e4897443df0f9629a2f18b74"
  ],
  [
    "orthohmm/search/evalue.py",
    [
      2026,
      9,
      26,
      20,
      58,
      32
    ],
    2176057344,
    1865,
    "b395cd48874d75963512b5d1e249c5251796f869592a52983c8d55b98fe11841"
  ],
  [
    "orthohmm/search/matrices.py",
    [
      2026,
      9,
      26,
      20,
      58,
      32
    ],
    2176057344,
    27982,
    "36e59a6c649a72f9c3a5d6eca117b7fe84d5ff8e975a1d7561a82ee8d6c29fcb"
  ],
  [
    "orthohmm/search/msa_center_star.py",
    [
      2026,
      9,
      26,
      20,
      58,
      32
    ],
    2176057344,
    9180,
    "3790801fc294f2df9ae3c84897de410b1a1bb487b8ea4226d86c388c379f40ab"
  ],
  [
    "orthohmm/search/msa_profile.py",
    [
      2026,
      9,
      26,
      20,
      58,
      32
    ],
    2176057344,
    6728,
    "d4076b61905325f643dca1f7a6b4dbd6f2833398840fad7b79e3f255f2240468"
  ],
  [
    "orthohmm/search/prefilter.py",
    [
      2026,
      9,
      26,
      20,
      58,
      32
    ],
    2176057344,
    22855,
    "9842b1787d98a89bd819f9adbe340b6e2a170d63d16057645720a6bc9faf6d06"
  ],
  [
    "orthohmm/search/profile.py",
    [
      2026,
      9,
      26,
      20,
      58,
      32
    ],
    2176057344,
    6143,
    "20dca3dea06a85bdcf99e03e3abfab921fce38b9d7715097b0b4588a7cd7587e"
  ],
  [
    "orthohmm/search/profile_expansion.py",
    [
      2026,
      9,
      26,
      20,
      58,
      32
    ],
    2176057344,
    25077,
    "1025cfab4f80e3a85d6eaa3e768ae5e0e25dedb48fdba9deb08cb6e6aa7bb42a"
  ],
  [
    "orthohmm/search/sequences.py",
    [
      2026,
      9,
      26,
      20,
      58,
      32
    ],
    2176057344,
    4509,
    "e083679cb9ea780664943355d260f398a96d93a2db336bb615af1032766f5980"
  ],
  [
    "orthohmm/search/species_prune.py",
    [
      2026,
      9,
      26,
      20,
      58,
      32
    ],
    2176057344,
    5324,
    "d729b9aa8ddf311c52e1d46c51d89dc5a6ef5d64ccd979deca90b7d980905b9e"
  ],
  [
    "orthohmm/search/viterbi.py",
    [
      2026,
      9,
      26,
      20,
      58,
      32
    ],
    2176057344,
    17693,
    "ba8543538d9bcbd55c5eaf5c7c3d5725a62c071b00d925f5fabaa474001948ba"
  ],
  [
    "orthohmm/search/viterbi_cuda_ctypes.py",
    [
      2026,
      9,
      26,
      20,
      58,
      32
    ],
    2176057344,
    4420,
    "fba688e77941c92c5e4080ef5de3dc23e4dfee6c5c2dad5d7ccdf876f81f9d20"
  ],
  [
    "orthohmm/search/csrc/hmm_viterbi.c",
    [
      2026,
      9,
      26,
      20,
      58,
      32
    ],
    2175008768,
    46012,
    "b668f86de525776e939f1c3196b41f3c870dd22cc56e35a93dff6a957d01f67c"
  ],
  [
    "orthohmm/search/csrc/hmm_viterbi.so",
    [
      2026,
      9,
      26,
      20,
      58,
      46
    ],
    2180841472,
    21008,
    "eeb4985e6f35689497a9c6187db56f8094c2cb9d92af18b77c027af0dfbea1d2"
  ],
  [
    "orthohmm/search/csrc/hmm_viterbi_cuda.cu",
    [
      2026,
      9,
      26,
      20,
      58,
      32
    ],
    2175008768,
    13696,
    "aa0a7a6fa53da025f28607a055364c08c365e561033a48dcf2c4381d350301d2"
  ],
  [
    "orthohmm/search/csrc/kmer_prefilter.c",
    [
      2026,
      9,
      26,
      20,
      58,
      32
    ],
    2175008768,
    10657,
    "8124e2652cb382549f2a282c3122aaa9f5c393da2b304b5d4f6cb92fae243f01"
  ],
  [
    "orthohmm/search/csrc/kmer_prefilter.so",
    [
      2026,
      9,
      26,
      20,
      58,
      46
    ],
    2180841472,
    20632,
    "4bef1865bccee6c4a5df6237880cf22b272cd567a1dc2d918815a28c86c7ce5c"
  ],
  [
    "orthohmm/search/csrc/pair_align.c",
    [
      2026,
      9,
      26,
      20,
      58,
      32
    ],
    2175008768,
    9098,
    "c5a3c2bc869ec59ca8261c976d6e77e5938faf5bfd8571566315b81cf5b08e9a"
  ],
  [
    "orthohmm/search/csrc/pair_align.so",
    [
      2026,
      9,
      26,
      20,
      58,
      46
    ],
    2180841472,
    20784,
    "1a1036e00730d9548d49819cdd852dc6aca74ace53c27ca2c903720904bc1d02"
  ],
  [
    "orthohmm-0.5.0.dist-info/licenses/LICENSE.md",
    [
      2026,
      9,
      26,
      20,
      58,
      44
    ],
    2175008768,
    1078,
    "4b9b0c3ffc73fce4c47734d7ef12bd70b9d2cdc4a3bf76ac0f3a35edbe2c3d9e"
  ],
  [
    "orthohmm-0.5.0.dist-info/METADATA",
    [
      2026,
      9,
      26,
      20,
      58,
      44
    ],
    2175008768,
    12275,
    "25ae0849520f91a134bb74067e0b0115c4a9e0985a54e0ae5cfb4e1ac4717f37"
  ],
  [
    "orthohmm-0.5.0.dist-info/WHEEL",
    [
      2026,
      9,
      26,
      20,
      58,
      46
    ],
    2176057344,
    104,
    "2a82caad1e573e5d036e9be26a7b0fdd6e63bf3d0faf62a4005331c0752301bc"
  ],
  [
    "orthohmm-0.5.0.dist-info/entry_points.txt",
    [
      2026,
      9,
      26,
      20,
      58,
      44
    ],
    2176057344,
    52,
    "fc529f6e2c11e3a70a98adc08bb56089a9590905865f0b4af7ab67f3d3fe4f6c"
  ],
  [
    "orthohmm-0.5.0.dist-info/top_level.txt",
    [
      2026,
      9,
      26,
      20,
      58,
      44
    ],
    2176057344,
    9,
    "1f55bee29a4dfa2cb6a86ef39d8e491d30df6df5297e147f811dce7fef3b6556"
  ],
  [
    "orthohmm-0.5.0.dist-info/RECORD",
    [
      2026,
      9,
      26,
      20,
      58,
      46
    ],
    2176057344,
    3531,
    "96c253d79cf6581caaab051b38c7f54ac2a2c00f2fa5ce03894b98d2a4b884bf"
  ]
]


def record(path):
    path = Path(path).absolute()
    data = path.read_bytes()
    return dict(path=str(path), bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def save(path, value):
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")


def pack(path, payloads):
    with zipfile.ZipFile(path, "x", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
        for name, date, attributes, _, _ in MEMBERS:
            info = zipfile.ZipInfo(name, tuple(date))
            info.create_system = 3
            info.create_version = info.extract_version = 20
            info.external_attr = attributes
            info.compress_type = zipfile.ZIP_DEFLATED
            archive.writestr(info, payloads[name], compress_type=zipfile.ZIP_DEFLATED, compresslevel=6)


def run(candidate, output):
    candidate, output = candidate.absolute(), output.absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    if (output.resolve() != output or candidate.is_symlink() or not candidate.is_file()
            or candidate.name != NAME or candidate.stat().st_size > 10 * 1024 ** 2):
        raise ValueError("Require a regular bounded candidate and fresh canonical destination")
    watched = [record(candidate), record(__file__)]
    payloads = {}
    with zipfile.ZipFile(candidate) as archive:
        names = archive.namelist()
        if len(names) != len(set(names)) or set(names) != {r[0] for r in MEMBERS}:
            raise ValueError("Candidate member inventory differs from the historical wheel")
        for name, _, _, size, sha in MEMBERS:
            if archive.getinfo(name).file_size != size:
                raise ValueError("Candidate payload size differs: " + name)
            data = archive.read(name)
            if hashlib.sha256(data).hexdigest() != sha:
                raise ValueError("Candidate payload differs: " + name)
            payloads[name] = data
    output.mkdir(parents=True)
    save(output / "started.json", dict(inputs=watched, expected=EXPECTED, members=len(MEMBERS),
        metadata_only_recipe=True, historical_admission=False, publication_ready=False))
    partial = output / (NAME + ".partial")
    try:
        pack(partial, payloads)
        ref = record(partial)
        if {k: ref[k] for k in ("bytes", "sha256")} != EXPECTED:
            raise ValueError("Reconstructed ZIP bytes differ; preserve attempt without retry")
        for item in watched:
            if record(item["path"]) != item:
                raise ValueError("Supplied candidate or reconstruction source changed")
        target = output / NAME
        if target.exists() or target.is_symlink():
            raise FileExistsError(target)
        partial.rename(target)
        target.chmod(0o644)
        result = dict(status="historical_project_wheel_byte_reconstructed", inputs=watched, wheel=record(target),
            matched_payload_members=len(MEMBERS), metadata_only_recipe=True, original_wheel_required=False,
            environment=dict(python=platform.python_version(), executable=sys.executable,
                zlib_build=zlib.ZLIB_VERSION, zlib_runtime=zlib.ZLIB_RUNTIME_VERSION),
            attempts=1, retry=False, historical_wheel_reproduced=True, historical_admission=False,
            native_code_executed=False, installation_performed=False, scientific_inference_executed=False,
            controlled_timing=False, publication_ready=False, security_clearance=False, redistribution_clearance=False,
            limitations=["Exact artifact reconstruction, not independent new numerical or inference validation.",
                "Only historical payload identities pass; another compiler's different binaries are not substituted.",
                "Compression differences fail closed; no universal Python/zlib/toolchain portability claim.",
                "No compiled artifact publication, runtime/OS/source-rights/security closure or complete study release."])
        save(output / "complete.json", result)
        return result
    except BaseException as error:
        save(output / "failed.json", dict(status="historical_project_wheel_reconstruction_failed",
            type=type(error).__name__, error=str(error), attempts=1, retry=False, publication_ready=False))
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.candidate, args.output), indent=2, sort_keys=True))
