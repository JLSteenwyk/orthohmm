from pathlib import Path

import pytest

from benchmark_tools.audit_fastoma_tasks import FIXED, IMAGE_DIGEST, PROGRAMS, audit_tasks, batch_name, collection_links, inspect_task, task_directory, validate_trace, wrapper_limits


def wrapper(cpu="1.0", memory="24g"):
    return f"docker run --cpus {cpu} --memory {memory} --network none {IMAGE_DIGEST} /bin/bash .command.sh\n"


def trace():
    names = sorted(FIXED) + ["omamer_run (a.fa)", "omamer_run (b.fa)", "hog_rest (1)"]
    return "task_id\thash\tname\tstatus\texit\n" + "".join(
        f"{i}\t01/{i:06x}\t{name}\tCOMPLETED\t0\n" for i, name in enumerate(names, 1))


def test_valid_trace():
    rows, counts = validate_trace(trace(), ["a.fa", "b.fa"])
    assert len(rows) == 9 and counts["omamer_run"] == 2


@pytest.mark.parametrize("before,after", [("COMPLETED", "CACHED"), ("COMPLETED", "FAILED"),
                                           ("\t0\n", "\t1\n"), ("hog_rest", "unknown"),
                                           ("01/000002", "01/000001"), ("2\t01/000002", "1\t01/000002")])
def test_bad_trace_rejected(before, after):
    with pytest.raises(ValueError):
        validate_trace(trace().replace(before, after, 1), ["a.fa", "b.fa"])


def test_missing_task_rejected():
    with pytest.raises(ValueError):
        validate_trace("\n".join(trace().splitlines()[:-1]), ["a.fa", "b.fa"])


def test_wrapper_limits():
    assert wrapper_limits(wrapper()) == {"cpus": 1.0, "memory_bytes": 24 * 1024**3}
    assert wrapper_limits(wrapper("16.0", "280g"))["cpus"] == 16


@pytest.mark.parametrize("cpu,memory", [("nan", "24g"), ("0", "24g"), ("181", "24g"),
                                        ("1", "0"), ("1", "701g"), ("1", "unlimited")])
def test_bad_limits(cpu, memory):
    with pytest.raises(ValueError):
        wrapper_limits(wrapper(cpu, memory))


@pytest.mark.parametrize("old,new", [(IMAGE_DIGEST, "fastoma:latest"), ("--network none", "--network host"),
                                     ("--cpus 1.0", "--cpus 1.0 --cpus 2.0"),
                                     ("docker run", "docker run --privileged")])
def test_wrapper_drift(old, new):
    with pytest.raises(ValueError):
        wrapper_limits(wrapper().replace(old, new))


def test_realistic_task_files(tmp_path):
    (tmp_path / ".command.sh").write_text("#!/bin/bash\nfastoma-check-input --proteomes proteome\n")
    (tmp_path / ".command.run").write_text(wrapper())
    (tmp_path / ".exitcode").write_text("0\n")
    (tmp_path / ".command.log").write_text("done\n")
    report, argv = inspect_task(tmp_path)
    assert len(report["files"]) == 4 and argv[0] == "fastoma-check-input"
    (tmp_path / ".exitcode").write_text("1")
    with pytest.raises(ValueError):
        inspect_task(tmp_path)


@pytest.mark.parametrize("value", ["../outside", "aa/../outside", "aa/123", "/tmp/abcdef"])
def test_bad_task_path(tmp_path, value):
    with pytest.raises(ValueError):
        task_directory(tmp_path, value)


def test_missing_and_ambiguous_task(tmp_path):
    with pytest.raises(ValueError):
        task_directory(tmp_path, "aa/123456")
    first = tmp_path / "aa/123456789"
    first.mkdir(parents=True)
    assert task_directory(tmp_path, "aa/123456") == first
    (tmp_path / "aa/123456abc").mkdir()
    with pytest.raises(ValueError):
        task_directory(tmp_path, "aa/123456")


def task_set(tmp_path):
    work = tmp_path / "work"
    path = tmp_path / "trace.txt"
    path.write_text(trace())
    rows, _ = validate_trace(trace(), ["a.fa", "b.fa"])
    paths = {}
    for row in rows:
        directory = work / (row["hash"] + "abcdef")
        directory.mkdir(parents=True)
        name = row["name"].split(" (", 1)[0]
        script = PROGRAMS[name]
        if name == "omamer_run":
            script += " search --query " + row["name"].split(" (")[1][:-1]
        elif name == "hog_rest":
            script += " --input-rhog-folder batch1"
        elif name == "extract_pairwise_ortholog_relations":
            script += " pw-rel --type ortholog"
        elif name == "batch_roothogs":
            (directory / "rhogs_rest/batch1").mkdir(parents=True)
        (directory / ".command.sh").write_text(script + "\n")
        (directory / ".command.run").write_text(wrapper())
        (directory / ".exitcode").write_text("0")
        (directory / ".command.log").write_text("done")
        paths[row["name"]] = directory
    return path, work, paths


def test_full_synthetic_task_chain(tmp_path):
    path, work, _ = task_set(tmp_path)
    report = audit_tasks(path, work, ["a.fa", "b.fa"])
    assert report["status"] == "fresh_fastoma_task_trace_verified"
    assert len(report["tasks"]) == 9


def test_absolute_batch_is_bound_to_exact_batching_directory(tmp_path):
    path, work, paths = task_set(tmp_path)
    batch = paths["batch_roothogs"] / "rhogs_rest/batch1"
    script = paths["hog_rest (1)"] / ".command.sh"
    script.write_text(f"fastoma-infer-subhogs --input-rhog-folder {batch}\n")
    assert audit_tasks(path, work, ["a.fa", "b.fa"])["process_counts"]["hog_rest"] == 1
    other = tmp_path / "other/batch1"
    other.mkdir(parents=True)
    script.write_text(f"fastoma-infer-subhogs --input-rhog-folder {other}\n")
    with pytest.raises(ValueError, match="batching output"):
        audit_tasks(path, work, ["a.fa", "b.fa"])


def test_batch_symlink_and_relative_escape_rejected(tmp_path):
    folder = tmp_path / "batches"
    folder.mkdir()
    other = tmp_path / "other"
    other.mkdir()
    (folder / "batch").symlink_to(other, target_is_directory=True)
    for value in (str(folder / "batch"), "../other", "missing"):
        with pytest.raises(ValueError):
            batch_name(value, folder)


@pytest.mark.parametrize("problem", [None, "missing", "duplicate", "foreign", "not_symlink", "broken"])
def test_exact_successful_collection(tmp_path, problem):
    output = tmp_path / "successful/pickle_hogs"
    output.mkdir(parents=True)
    collector = tmp_path / "collector"
    folder = collector / "pickle_folders"
    folder.mkdir(parents=True)
    link = folder / "1"
    link.symlink_to(output, target_is_directory=True)
    if problem == "missing":
        link.unlink()
    elif problem == "duplicate":
        (folder / "2").symlink_to(output, target_is_directory=True)
    elif problem in {"foreign", "broken"}:
        link.unlink()
        other = tmp_path / "failed/pickle_hogs"
        if problem == "foreign":
            other.mkdir(parents=True)
        link.symlink_to(other, target_is_directory=True)
    elif problem == "not_symlink":
        link.unlink()
        link.mkdir()
    if problem:
        with pytest.raises(ValueError):
            collection_links(collector, {str(output)})
    else:
        assert collection_links(collector, {str(output)}) == [{"path": str(link), "target": str(output)}]


@pytest.mark.parametrize("change", ["query", "batch", "program", "pair_type", "missing_batch_task"])
def test_full_chain_mismatch_rejected(tmp_path, change):
    path, work, paths = task_set(tmp_path)
    if change == "query":
        script = paths["omamer_run (b.fa)"] / ".command.sh"
        script.write_text("omamer search --query a.fa")
    elif change == "batch":
        (paths["hog_rest (1)"] / ".command.sh").write_text("fastoma-infer-subhogs --input-rhog-folder other")
    elif change == "program":
        (paths["check_input"] / ".command.sh").write_text("echo skipped")
    elif change == "pair_type":
        (paths["extract_pairwise_ortholog_relations"] / ".command.sh").write_text("fastoma-helper pw-rel --type paralog")
    else:
        (paths["batch_roothogs"] / "rhogs_rest/batch2").mkdir()
    with pytest.raises(ValueError):
        audit_tasks(path, work, ["a.fa", "b.fa"])
