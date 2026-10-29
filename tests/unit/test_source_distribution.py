"""Check the source archive users install, including native build inputs."""

import os
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import sys
import tarfile


ROOT = Path(__file__).resolve().parents[2]


def test_source_distribution_excludes_non_application_payload(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    shutil.copytree(ROOT / "orthohmm", source / "orthohmm",
                    ignore=shutil.ignore_patterns("*.so", "__pycache__", "*.pyc"))
    for name in ("setup.py", "MANIFEST.in", "README.md", "LICENSE.md", "CITATION.cff", "requirements.txt"):
        shutil.copyfile(ROOT / name, source / name)
    for directory in ("tests", "docs", "scripts", "manuscript", "benchmark_tools"):
        path = source / directory
        path.mkdir()
        (path / "unrelated.py").write_text("must_not_be_distributed = True\n")
    (source / "orthohmm/search/csrc/inherited.so").write_bytes(b"stale binary")
    environment = {key: value for key, value in os.environ.items()
                   if not key.startswith(("COVERAGE_", "COV_CORE_"))}
    result = subprocess.run([sys.executable, "setup.py", "sdist", "--dist-dir", str(tmp_path / "dist")],
                            cwd=source, env=environment, capture_output=True, text=True, timeout=90)
    assert result.returncode == 0, result.stdout + result.stderr
    archives = list((tmp_path / "dist").glob("*.tar.gz"))
    assert len(archives) == 1
    with tarfile.open(archives[0]) as archive:
        names = {PurePosixPath(*PurePosixPath(member.name).parts[1:]).as_posix()
                 for member in archive.getmembers() if member.isfile()}
    forbidden = {"tests", "docs", "scripts", "manuscript", "benchmark_tools", "benchmarks"}
    assert not {name for name in names if name.split("/")[0] in forbidden}
    assert not {name for name in names if name.endswith((".so", ".pyc", ".dylib", ".dll"))}
    assert {"README.md", "LICENSE.md", "CITATION.cff", "orthohmm/__main__.py"} <= names
    assert {"orthohmm/search/csrc/" + name for name in
            ("hmm_viterbi.c", "kmer_prefilter.c", "pair_align.c", "hmm_viterbi_cuda.cu")} <= names
