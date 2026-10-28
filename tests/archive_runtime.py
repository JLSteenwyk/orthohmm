"""Extract only regular test-package files on Python versions before PEP 706."""

from pathlib import Path
import shutil


def extract_regular_files(archive, destination):
    destination = Path(destination).resolve()
    members = archive.getmembers()
    for member in members:
        target = destination / member.name
        if (not (member.isfile() or member.isdir()) or Path(member.name).is_absolute()
                or not target.resolve().is_relative_to(destination)):
            raise ValueError("Unsafe test archive member")
    for member in members:
        target = destination / member.name
        if member.isdir():
            target.mkdir(parents=True, exist_ok=True)
        else:
            target.parent.mkdir(parents=True, exist_ok=True)
            with archive.extractfile(member) as source, target.open("xb") as output:
                shutil.copyfileobj(source, output)
