import csv
import io

import pytest

from benchmark_tools.review_fastoma_retries import partition_attempts, review
from benchmark_tools.audit_fastoma_tasks import audit_tasks, task_directory
from tests.unit.test_audit_fastoma_tasks import task_set, trace, wrapper


def retry_trace():
    rows = list(csv.DictReader(io.StringIO(trace()), delimiter="\t"))
    failed = rows[-1]
    successful = {**failed, "task_id": "10", "hash": "02/000010"}
    failed.update(status="FAILED", exit="1")
    rows.append(successful)
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue()


PAIRS = {"01/000009": "02/000010"}


def test_explicit_partition_retains_failed_attempt():
    rows, failed, successful, counts = partition_attempts(retry_trace(), ["a.fa", "b.fa"], PAIRS)
    assert len(rows) == 10 and len(failed) == 1 and len(successful) == 9
    assert counts["hog_rest"] == 1


@pytest.mark.parametrize("old,new", [("FAILED", "CACHED"), ("FAILED", "RUNNING"),
    ("FAILED\t1", "FAILED\t137"), ("10\t02/000010", "9\t02/000010"),
    ("02/000010", "01/000009"), ("10\t02/000010\thog_rest (1)", "10\t02/000010\thog_rest (2)")])
def test_unreviewed_or_ambiguous_attempts_rejected(old, new):
    with pytest.raises(ValueError):
        partition_attempts(retry_trace().replace(old, new), ["a.fa", "b.fa"], PAIRS)


@pytest.mark.parametrize("pairs", [{}, {"01/000009": "02/999999"},
    {"01/999999": "02/000010"}, {"01/000009": "01/000009"}])
def test_review_pair_must_match_exact_observed_attempts(pairs):
    with pytest.raises(ValueError):
        partition_attempts(retry_trace(), ["a.fa", "b.fa"], pairs)


def fixture(tmp_path):
    path, work, paths = task_set(tmp_path)
    path.write_text(retry_trace())
    before = paths["hog_rest (1)"]
    after = work / "02/000010abcdef"
    after.mkdir(parents=True)
    tree = tmp_path / "checked.nwk"
    tree.write_text("(a,b);")
    for directory, code, memory in ((before, "1", "12g"), (after, "0", "24g")):
        (directory / ".command.sh").write_text("fastoma-infer-subhogs --species-tree tree.nwk --input-rhog-folder batch1\n")
        (directory / ".command.run").write_text(wrapper(memory=memory))
        (directory / ".exitcode").write_text(code)
        (directory / ".command.err").write_text("failure" if code == "1" else "done")
        (directory / ".command.log").write_text("retained")
        (directory / "tree.nwk").symlink_to(tree)
    return path, work, before, after


def test_file_review_retains_both_attempts_without_admitting(tmp_path):
    path, work, _, _ = fixture(tmp_path)
    result = review(path, work, ["a.fa", "b.fa"], PAIRS)
    assert len(result["attempts"]) == 10
    assert result["retry_pairs"][0]["failed"]["trace"]["exit"] == "1"
    assert result["native_outputs_admitted"] is False
    assert result["accuracy_evaluated"] is False


def test_integrated_explicit_audit_keeps_default_closed(tmp_path):
    path, work, before, after = fixture(tmp_path)
    for directory in (before, after):
        script = directory / ".command.sh"
        script.write_text(script.read_text().rstrip() + " --output-pickles pickle_hogs\n")
    (after / "pickle_hogs").mkdir()
    rows = list(csv.DictReader(io.StringIO(trace()), delimiter="\t"))
    row = next(r for r in rows if r["name"] == "collect_subhogs")
    collector = task_directory(work, row["hash"])
    (collector / "pickle_folders").mkdir()
    (collector / "pickle_folders/1").symlink_to(after / "pickle_hogs", target_is_directory=True)
    with pytest.raises(ValueError, match="explicit review"):
        audit_tasks(path, work, ["a.fa", "b.fa"])
    report = audit_tasks(path, work, ["a.fa", "b.fa"], retry_pairs=PAIRS)
    assert report["status"] == "explicitly_reviewed_fastoma_task_trace_verified"
    assert len(report["tasks"]) == 9
    assert len(report["retry_review"]["attempts"]) == 10
    assert len(report["collection"]) == 1


@pytest.mark.parametrize("change", ["command", "tree", "exit", "cpu", "memory"])
def test_file_drift_rejected(tmp_path, change):
    path, work, _, after = fixture(tmp_path)
    if change == "command":
        (after / ".command.sh").write_text("fastoma-infer-subhogs --species-tree tree.nwk --different\n")
    elif change == "tree":
        (after / "tree.nwk").unlink()
        (after / "tree.nwk").write_text("(b,a);")
    elif change == "exit":
        (after / ".exitcode").write_text("1")
    else:
        (after / ".command.run").write_text(wrapper(cpu="2" if change == "cpu" else "1", memory="36g" if change == "memory" else "24g"))
    with pytest.raises(ValueError):
        review(path, work, ["a.fa", "b.fa"], PAIRS)
