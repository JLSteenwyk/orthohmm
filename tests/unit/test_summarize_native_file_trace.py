import pytest

from benchmark_tools.summarize_native_file_trace import summarize
from benchmark_tools.probe_threadripper_file_access import traced_command


def test_opens_failures_and_mutations():
    result = summarize([
        '1 openat(AT_FDCWD</tmp>, "x", O_WRONLY|O_CREAT, 0600) = 3</tmp/x>\n',
        '2 openat(AT_FDCWD</tmp>, "x", O_RDONLY) = 3</tmp/x>\n',
        '1 openat(AT_FDCWD, "absent", O_RDONLY) = -1 ENOENT (No such file)\n',
        '1 unlink("/tmp/x") = 0\n',
        '1 execve("/bin/tool", ["tool"], 0x10) = 0\n',
        '1 +++ exited with 0 +++\n'])
    assert result["opened_paths"] == [dict(path="/tmp/x", successful_opens=2,
        write_capable_opens=1, first_line=1)]
    assert len(result["successful_mutations"]) == len(result["successful_executions"]) == 1
    assert result["unparsed"] == [] and result["unfinished"] == {}


def test_interleaved_resume_and_incomplete():
    result = summarize([
        '1 openat(AT_FDCWD, "/tmp/a", O_RDWR <unfinished ...>\n',
        '2 openat(AT_FDCWD, "/tmp/b", O_RDONLY) = 3</tmp/b>\n',
        '1 <... openat resumed>) = 4</tmp/a>\n',
        '2 wait4(-1,  <unfinished ...>\n'])
    assert [p["path"] for p in result["opened_paths"]] == ["/tmp/a", "/tmp/b"]
    assert result["unfinished"]
    assert result["opened_paths"][0]["write_capable_opens"] == 1


def test_unknown_open_is_not_silent():
    result = summarize(['1 open("a", O_RDWR) = 3\n', 'not decoded\n'])
    assert len(result["unparsed"]) == 2
    with pytest.raises(ValueError):
        summarize(['1 <... openat resumed>) = 3</tmp/a>\n'])


def test_device_annotation():
    result = summarize(['1 openat(AT_FDCWD, "/dev/null", O_RDWR) = 3</dev/null<char 1:3>>\n'])
    assert result["opened_paths"][0]["path"] == "/dev/null"
    assert result["unparsed"] == []


def test_pipe_descriptor_is_not_a_filesystem_path():
    result = summarize(['1 openat(AT_FDCWD, "/dev/stderr", O_WRONLY) = 3<pipe:[123]>\n'])
    assert result["opened_paths"] == [] and result["unparsed"] == []
    assert result["nonfilesystem_opens"][0]["descriptor"] == "pipe:[123]"


def test_tracer_does_not_rewrite_native_arguments():
    native = ["/python", "-m", "orthohmm", "/input with spaces", "-c", "32"]
    wrapped = traced_command(native, "/trace with spaces")
    assert wrapped[wrapped.index("--") + 1:] == native
    assert wrapped[0] == "/usr/bin/strace"
