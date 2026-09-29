from benchmark_tools.bind_private_timing_runtime import root_partition


def test_retirement_is_prefix_scoped_and_deduplicated():
    roots = ["/home/bizon/anaconda3/lib", "/home/bizon/anaconda3-other/lib", "/usr/lib",
             "/home/bizon/.local/lib/python3.10/site-packages", "/usr/lib",
             "/home/bizon/.cache/matplotlib/fontlist.json", "/private/venv"]
    retained, retired = root_partition(roots)
    assert retained == ["/home/bizon/anaconda3-other/lib", "/private/venv", "/usr/lib"]
    assert len(retired) == 3


def test_non_python_user_paths_not_silently_removed():
    assert root_partition(["/home/bizon/other"])[0] == ["/home/bizon/other"]
