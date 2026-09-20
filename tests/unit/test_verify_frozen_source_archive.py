import io
from pathlib import Path
import subprocess
import tarfile

import pytest

from benchmark_tools import verify_frozen_source_archive as module

ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture(scope="module")
def sources():
    return module.git_files(ROOT)


def make_archive(tmp_path, sources, fault=None):
    path = tmp_path / "source.tar.gz"
    with tarfile.open(path, "w:gz") as archive:
        for name, row in sources.items():
            if fault == "missing" and name == "setup.py":
                continue
            content = row["content"]
            member = tarfile.TarInfo(module.PREFIX + name)
            member.mode = row["mode"]
            if name == "setup.py":
                if fault == "content":
                    content += b"\n# altered\n"
                elif fault == "mode":
                    member.mode = 0o777
            member.size = len(content)
            archive.addfile(member, io.BytesIO(content))
        if fault in ("extra", "duplicate", "link", "traversal"):
            name = {"extra": "binary.so", "duplicate": "setup.py", "link": "link", "traversal": "../escape"}[fault]
            member = tarfile.TarInfo(module.PREFIX + name)
            if fault == "link":
                member.type, member.linkname = tarfile.SYMTYPE, "setup.py"
            archive.addfile(member)
    return path


def test_complete_frozen_archive(tmp_path, sources):
    result = module.verify(ROOT, make_archive(tmp_path, sources))
    assert len(result["files"]) == len(sources)
    assert result["python_files_syntax_checked"] == sum(n.endswith(".py") for n in sources)
    assert result["executable_benchmark_reproduced"] is False
    assert result["publication_ready"] is False


@pytest.mark.parametrize("fault", ["missing", "content", "mode", "extra", "duplicate", "link", "traversal"])
def test_corrupt_or_expanded_archive_rejected(tmp_path, sources, fault):
    with pytest.raises(ValueError):
        module.verify(ROOT, make_archive(tmp_path, sources, fault))


def test_cli_does_not_overwrite(tmp_path, sources):
    archive = make_archive(tmp_path, sources)
    output = tmp_path / "report.json"
    argv = ["--repo", str(ROOT), "--archive", str(archive), "--output", str(output)]
    assert module.main(argv) == 0
    original = output.read_bytes()
    with pytest.raises(FileExistsError):
        module.main(argv)
    assert output.read_bytes() == original


def test_real_git_archive_with_explicit_permissions(tmp_path):
    archives = [tmp_path / name for name in ("one.tar.gz", "two.tar.gz")]
    for path in archives:
        subprocess.run(["git", "-C", str(ROOT), "-c", "tar.umask=0022", "archive",
            "--format=tar.gz", "--prefix=" + module.PREFIX, "--output=" + str(path),
            module.REVISION, *module.PATHS], check=True)
    assert archives[0].read_bytes() == archives[1].read_bytes()
    assert module.verify(ROOT, archives[0])["status"] == "frozen_source_archive_matches_git"
