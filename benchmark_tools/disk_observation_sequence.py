"""Lazy views of numbered raw observations; no retained decoded history."""

from collections.abc import Sequence
import json
from pathlib import Path

from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.prepare_ob_candidate_neighborhood import record


class DiskObservations(Sequence):
    def __init__(self, directory, count=0, *, indices=None):
        if type(count) is not int or count < 0:
            raise ValueError("Require nonnegative observation count")
        self.directory = Path(directory)
        self.indices = range(count) if indices is None else indices
        self.writable = indices is None

    def __len__(self):
        return len(self.indices)

    def path(self, index):
        return self.directory / f"point_{self.indices[index]:06d}.json"

    def __getitem__(self, index):
        if isinstance(index, slice):
            return DiskObservations(self.directory, indices=self.indices[index])
        path = self.path(index)
        if path.is_symlink():
            raise ValueError("Observation must not be a symlink")
        return json.loads(path.read_text())

    def append(self, point):
        if not self.writable:
            raise ValueError("Cannot append to an observation view")
        path = self.directory / f"point_{len(self):06d}.json"
        save(path, point)
        self.indices = range(len(self) + 1)

    def records(self):
        return [record(self.path(i)) for i in range(len(self))]
