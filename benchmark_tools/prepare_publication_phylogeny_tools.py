"""Build/stage pinned native phylogeny tools offline without old installations.

Use an externally verified source component and public artifacts. This runs
the fixed MAFFT core make target and version probes, not scientific inference.
"""

import argparse
import json
import os
from pathlib import Path
import platform
import shutil
import string
import sys

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmark_tools import bundle_publication_source as source
from benchmark_tools import acquire_publication_phylogeny_tools as acquisition
from benchmark_tools.build_publication_fasttree import validate_help
from benchmark_tools.run_integrated_publication_workflow import record, save, stage

# Metadata only, derived from the admitted 2026-09-26 MAFFT source-build receipt.
HELPERS = {
    "addsingle": dict(bytes=715016, sha256="f1dc4d582a6dd5bc0ae2ee4d864bfb3feb1236db19dd174f2ea1997a5745e84f"),
    "contrafoldwrap": dict(bytes=370952, sha256="9c5c4b0ac2a9b7a9621d33fcd89738bd79b02df8636daea39ba46b751e1cc02f"),
    "countlen": dict(bytes=108920, sha256="38827b4193765e75ee254930b4ba443bcd050200df973891a1771979ad30e41e"),
    "disttbfast": dict(bytes=731400, sha256="91ad64cfbd13cefc90f6f2aa47a0f25979a4a35775be4fcefe4b586b66ada52c"),
    "dndblast": dict(bytes=387336, sha256="f774cd87f0e5f492c1296761515de43cd63f7fa2b5172578ac7e1fa3a3723a12"),
    "dndfast7": dict(bytes=383240, sha256="63786abcec720442f8845d761c70ae8689f7c113f51f65070bd38c5678628071"),
    "dndpre": dict(bytes=370952, sha256="4e8a89c1e725fa9d463c94f636993118fb79290af2e5dc5509a2bb962b413121"),
    "dvtditr": dict(bytes=731400, sha256="5ea9386614ae82f2915e0f4be3f98ff4c3e3739c10dc938cc5c9b7a8c265cc2f"),
    "f2cl": dict(bytes=182536, sha256="1cf39f5063fbd478e34bf470aa4b6db491ec4bf607137c43e9ee5383b1fa51c5"),
    "filter": dict(bytes=297336, sha256="56f72b8ff9a3878fa00e283425218b5a76716cc42e07689dd611dab66e82b80b"),
    "getlag": dict(bytes=649480, sha256="a788bc8e0eb6a35f0693e98da690533ec57cfc353f2b7f809c30c9fada4f62b8"),
    "hex2maffttext": dict(bytes=14472, sha256="069e138d9dad49ac6be2ca98dcc5aad22aed5de54542f217c6c08d10c5f8468f"),
    "mafft-distance": dict(bytes=383240, sha256="54d1dceea0454a5cf082fcb00d427d0ea6f8abe86c82f047a1c2fa3f201216c2"),
    "mafft-profile": dict(bytes=645384, sha256="d3de6cde0d4dce0887dd248244c26f1b152ad2de7efe19b4346b221841bdf081"),
    "mafftash_premafft.pl": dict(bytes=10475, sha256="0594ba064e06ad0b647a9927bc73541edd85fd18c09a65c7a5c849b318b93ef9"),
    "maffttext2hex": dict(bytes=14472, sha256="d50a29efa51a002700d4013229b34cc70b6659aef8a24e1b6b39eab121450067"),
    "makedirectionlist": dict(bytes=682264, sha256="b8a9458b09087f668f48fbdba20494d896bc98c33a215fa2590f03b826d5d440"),
    "mccaskillwrap": dict(bytes=370952, sha256="a4095c4b9da4d82818dea66e3018953e24e8862844dc9640a83efff3dd72e1f2"),
    "multi2hat3s": dict(bytes=395528, sha256="d77137be47c0100006fb9b72484189e8b1da3c1ad2a7c7338a8608e784811a28"),
    "nodepair": dict(bytes=739592, sha256="987e5d65c7668d98052adcf94262dc8d4c8c8a5fc7d5ff46acf9addde7c0e935"),
    "pairash": dict(bytes=669960, sha256="73184ab0a93deca54524f0e600e4f4cddfb085f34bad4c7045563b9ab78cafea"),
    "pairlocalalign": dict(bytes=706824, sha256="268ad1d72740d73d2298731f3769c37ed12fc76924ca0974e5de9a28f2ee852d"),
    "regtable2seq": dict(bytes=297336, sha256="7e5800124365538d9b0bab80df4a086f169a59850a45fba9a1fa65b146a703f9"),
    "replaceu": dict(bytes=297336, sha256="5638e6fd7278e3a80563659eb53f90e986eb5b8ea468c28676cd8074d2631383"),
    "restoreu": dict(bytes=297336, sha256="80c587171f2cc6c93cae3018f0ae1fc97c01032d4ecaddbf75d962d2886704e1"),
    "score": dict(bytes=366856, sha256="92e392b2e0684e579d9534b11cdaa8b8eb3a2d69fb6b3afa9ad9ed53a8e5db9b"),
    "seekquencer_premafft.pl": dict(bytes=15541, sha256="bfe2f4cd579c3734ac9a4979c9a1b577bef8993f5242c342019df2c490d82fd3"),
    "seq2regtable": dict(bytes=108920, sha256="3f268e55b08c30ed77d7f7ea2584c5f774058b2fe45234ef02c9317728bc52d4"),
    "setcore": dict(bytes=645384, sha256="2e749598ca68caa44eb5d5c600fe1da2f5bfe2b91b995401a140c7d8d137a046"),
    "setdirection": dict(bytes=297344, sha256="a86ffca59d331c198407f63adfc2815e8584f7d3064d552b8832b723b108dd04"),
    "sextet5": dict(bytes=383240, sha256="34dde941f162ad2abc0489997b86120d83c0bc15eb09c5d01004616e9393a819"),
    "splittbfast": dict(bytes=674056, sha256="c36661cfab0336ce26b10adbbc21e27ecf32373a2747ca1f3dd2cd7a711d3737"),
    "tbfast": dict(bytes=784648, sha256="f5f77a3ba787e06b663f7b63d363131c59ae1f1cea752a0b1af311bbcbb443b5"),
    "version": dict(bytes=14472, sha256="d454117bae91fa575f6092cc381753953784c190e1bbb11de40f1241070e4953"),
}
SOURCE_FILES = 173


def regular(path):
    path = Path(path)
    if path.is_symlink() or not path.is_file():
        raise ValueError("Require a nonsymlink regular file: " + str(path))
    return record(path)


def cpu_compatible():
    if platform.system() != "Linux" or platform.machine() != "x86_64":
        raise ValueError("Historical tool preparation requires Linux x86-64")
    flags = [set(value.split()) for line in Path("/proc/cpuinfo").read_text().splitlines()
             for key, separator, value in [line.partition(":")] if separator and key.strip() == "flags"]
    if not flags or not all("avx2" in row for row in flags):
        raise ValueError("Require AVX2 support for the retained upstream FastTree binary")
    return dict(system="Linux", machine="x86_64", required_feature="avx2", observed=True)


def preflight(args):
    output = args.output.absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    if (output.resolve() != output or type(args.timeout) is not int or args.timeout < 1
            or args.acknowledge_historical_runtime is not True):
        raise ValueError("Require canonical fresh output, positive timeout and historical acknowledgement")
    # Upstream make recipes interpolate PREFIX without quoting every shell use.
    if any(char not in string.ascii_letters + string.digits + "/._-" for char in str(output)):
        raise ValueError("Upstream make requires a shell-safe output path without spaces")
    cpu = cpu_compatible()
    component = args.component.resolve(strict=True)
    artifacts = args.artifacts.absolute()
    if artifacts.resolve() != artifacts or not artifacts.is_dir():
        raise ValueError("Require a canonical artifact directory")
    if output.is_relative_to(component) or output.is_relative_to(artifacts):
        raise ValueError("Output must be outside the source component and acquired artifacts")
    verified = source.verify(component, args.manifest_sha256)
    rows = acquisition.artifacts()
    files = []
    for row in rows:
        path = artifacts / row["relative"]
        if path.resolve() != path:
            raise ValueError("Artifact path alias/link is not admitted")
        item = regular(path)
        if any(item[key] != row[key] for key in ("bytes", "sha256")):
            raise ValueError("Changed frozen native artifact: " + row["relative"])
        files.append(item)
    tools = {}
    for name in ("gcc", "make", "ld"):
        entry = Path("/usr/bin") / name
        if not entry.is_file() or not os.access(entry, os.X_OK):
            raise ValueError("Require separately supplied system toolchain: " + name)
        tools[name] = entry.resolve(strict=True)
    if record(tools["gcc"])["sha256"] != args.compiler_sha256:
        raise ValueError("Compiler identity differs from external anchor")
    watched = [*files, *acquisition.source_records(), regular(__file__),
               regular(source.__file__), regular(sys.modules[stage.__module__].__file__),
               regular(component / "SOURCE_INDEX.json"), *[regular(path) for path in tools.values()]]
    return output, component, artifacts, verified, cpu, tools, watched


def inspect_helpers(prefix):
    directory = prefix / "libexec/mafft"
    if directory.is_symlink() or not directory.is_dir():
        raise ValueError("Missing regular MAFFT helper directory")
    if {path.name for path in directory.iterdir()} != set(HELPERS):
        raise ValueError("Changed MAFFT helper inventory")
    rows = []
    for name, expected in sorted(HELPERS.items()):
        path = directory / name
        item = regular(path)
        if ({key: item[key] for key in ("bytes", "sha256")} != expected
                or path.stat().st_mode & 0o777 != 0o755):
            raise ValueError("Built MAFFT helper differs from historical bytes/mode: " + name)
        if not name.endswith(".pl") and path.read_bytes()[:4] != b"\x7fELF":
            raise ValueError("Required MAFFT native helper is not ELF: " + name)
        rows.append(dict(**item, mode=0o755))
    return rows


def relative_links(prefix):
    rows = []
    for name in ("mafft-distance", "mafft-profile"):
        path = prefix / "bin" / name
        expected = prefix / "libexec/mafft" / name
        if not path.is_symlink() or os.readlink(path) != str(expected) or path.resolve() != expected:
            raise ValueError("Unexpected generated MAFFT convenience link")
        previous = os.readlink(path)
        path.unlink()
        target = "../libexec/mafft/" + name
        path.symlink_to(target)
        rows.append(dict(relative=path.relative_to(prefix).as_posix(), previous=previous, target=target,
                         payload=regular(expected)))
    for path in prefix.rglob("*"):
        if path.is_symlink() and (not path.resolve().is_relative_to(prefix) or not path.resolve().is_file()):
            raise ValueError("Installed tool link escapes private prefix or is broken")
    return rows


def inventory(directory):
    result = []
    for path in sorted(directory.rglob("*")):
        if path.is_symlink():
            result.append(dict(relative=path.relative_to(directory).as_posix(),
                               kind="symlink", target=os.readlink(path)))
        elif path.is_file():
            result.append(dict(relative=path.relative_to(directory).as_posix(), kind="file",
                               mode=path.stat().st_mode & 0o777, **regular(path)))
        elif not path.is_dir():
            raise ValueError("Nonregular installed tool entry")
    return result


def run(args):
    output, component, artifacts, verified, cpu, tools, watched = preflight(args)
    output.mkdir(parents=True)
    for name in ("home", "tmp"):
        (output / name).mkdir()
    environment = dict(PATH="/usr/bin:/bin", HOME=str(output / "home"), TMPDIR=str(output / "tmp"),
                       LANG="C.UTF-8", LC_ALL="C", MAKEFLAGS="", MFLAGS="", GNUMAKEFLAGS="",
                       OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    scope = dict(attempts=1, retry=False, historical_runtime=True, old_installation_required=False,
                 private_tools_only=True, shared_environment_modified=False,
                 scientific_inference_executed=False, historical_admission=False, controlled_timing=False,
                 security_clearance=False, redistribution_clearance=False, publication_ready=False)
    save(output / "started.json", dict(inputs=watched, source=verified, cpu=cpu,
                                      environment=environment, native_execution_planned=True, **scope))
    outcomes = []
    try:
        archive = artifacts / acquisition.artifacts()[0]["relative"]
        sources = acquisition.mafft.unpack(archive, output / "source")
        if len(sources) != SOURCE_FILES:
            raise ValueError("Wrong frozen MAFFT source inventory")
        core = output / "source" / acquisition.mafft.ROOT / "core"
        source_root = core.parent
        prefix = output / "tools/mafft"
        outcomes.append(stage(output, "compiler_version", [str(tools["gcc"]), "--version"],
                              environment, min(args.timeout, 30)))
        command = [str(tools["make"]), "-C", str(core), "-j2", "CC=" + str(tools["gcc"]),
                   "CFLAGS=-O3", "PREFIX=" + str(prefix), "install"]
        outcomes.append(stage(output, "mafft_core_build", command, environment, args.timeout))
        helpers = inspect_helpers(prefix)
        links = relative_links(prefix)
        launcher = regular(prefix / "bin/mafft")
        notices = output / "tools/notices/mafft"
        notices.mkdir(parents=True)
        notice_records = []
        for name in ("license", "license.extensions", "README.md"):
            original = regular(source_root / name)
            target = notices / name
            shutil.copyfile(original["path"], target)
            target.chmod(0o644)
            item = regular(target)
            if any(item[key] != original[key] for key in ("bytes", "sha256")):
                raise ValueError("Changed copied MAFFT source/notice")
            notice_records.append(item)
        fasttree = output / "tools/fasttree"
        fasttree.mkdir()
        for row in acquisition.artifacts()[1:]:
            original = artifacts / row["relative"]
            target = fasttree / original.name
            shutil.copyfile(original, target)
            target.chmod(0o755 if target.name == "FastTree" else 0o644)
            item = regular(target)
            if any(item[key] != row[key] for key in ("bytes", "sha256")):
                raise ValueError("Changed copied FastTree artifact")
        probe_environment = dict(environment, MAFFT_BINARIES=str(prefix / "libexec/mafft"))
        outcomes.append(stage(output, "mafft_version", [launcher["path"], "--version"],
                              probe_environment, min(args.timeout, 30)))
        if (output / "mafft_version.log").read_text().strip() != "v7.525 (2024/Mar/13)":
            raise ValueError("Unexpected MAFFT launcher version")
        outcomes.append(stage(output, "mafft_helper_version", [str(prefix / "libexec/mafft/version")],
                              probe_environment, min(args.timeout, 30)))
        if (output / "mafft_helper_version.log").read_text().strip() != "7.525":
            raise ValueError("Unexpected compiled MAFFT version")
        outcomes.append(stage(output, "fasttree_help", [str(fasttree / "FastTree"), "-help"],
                              probe_environment, min(args.timeout, 30)))
        validate_help(0, (output / "fasttree_help.log").read_text(), "")
        installed = inventory(output / "tools")
        for item in [*watched, *sources, *(row for row in installed if row["kind"] == "file")]:
            actual = regular(item["path"])
            if any(actual[key] != item[key] for key in ("path", "bytes", "sha256")):
                raise ValueError("Input/source/installed tool changed during preparation")
        if inventory(output / "tools") != installed or inspect_helpers(prefix) != helpers:
            raise ValueError("Prepared tools changed")
        if any(Path("/usr/bin", name).resolve() != path for name, path in tools.items()):
            raise ValueError("System toolchain alias changed")
        source.verify(component, args.manifest_sha256)
        result = dict(status="private_frozen_phylogeny_tools_prepared", inputs=watched, source=verified,
                      source_files=sources, helper_files=helpers, helpers_historical_byte_equal=True,
                      launcher=launcher, link_changes=links, notices=notice_records, inventory=installed,
                      stages=outcomes, cpu=cpu, environment=environment, probe_environment=probe_environment,
                      mafft=str(prefix / "bin/mafft"), fasttree=str(fasttree / "FastTree"),
                      source_unchanged=True, native_code_executed=True, **scope, limitations=[
                          "Same recorded compiler/host reconstruction, not a hermetic or cross-platform build.",
                          "The generated launcher embeds its fresh private prefix and is not the historical launcher bytes.",
                          "After relocation, MAFFT_BINARIES must explicitly identify the relocated helper directory.",
                          "FastTree is the exact upstream binary, not a newly source-reproduced baseline binary.",
                          "Version probes are not numerical alignment/tree or full benchmark-equivalence validation.",
                          "Optional MAFFT extension engines are not built; source and notices remain preserved.",
                          "OS/toolchain/security/rights closure, runtime assembly and scientific/timing admission remain open.",
                          "No network acquisition, old installation, scientific fixture rerun or public binary upload.",
                      ])
        save(output / "complete.json", result)
        return result
    except BaseException as error:
        save(output / "failed.json", dict(status="private_phylogeny_preparation_failed", inputs=watched,
                  stages=outcomes, error_type=type(error).__name__, error=str(error), **scope))
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--component", type=Path, required=True)
    parser.add_argument("--manifest-sha256", required=True)
    parser.add_argument("--artifacts", type=Path, required=True)
    parser.add_argument("--compiler-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--timeout", type=int, default=900)
    parser.add_argument("--acknowledge-historical-runtime", action="store_true")
    args = parser.parse_args()
    result = run(args)
    print(json.dumps(dict(status=result["status"], complete=record(args.output / "complete.json")), sort_keys=True))


if __name__ == "__main__":
    main()
