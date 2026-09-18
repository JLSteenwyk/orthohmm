import json
import os
from pathlib import Path
import shutil
import subprocess

import pytest

WRAPPER = Path(__file__).resolve().parents[2] / "benchmark_tools/run_orthomcl_perl_script.pl"
PERL = shutil.which("perl")
pytestmark = pytest.mark.skipif(PERL is None, reason="Perl required")


def test_native_script_identity_arguments_and_module_path(tmp_path):
    script = tmp_path / "script.pl"
    script.write_text('use JSON::PP; print encode_json({name=>$0,args=>\\@ARGV,inc=>\\@INC});\n')
    done = subprocess.run([PERL, "-I.", str(WRAPPER), str(script), "one", "two words"],
                          cwd=tmp_path, capture_output=True, text=True, check=True)
    result = json.loads(done.stdout)
    assert result["name"] == str(script)
    assert result["args"] == ["one", "two words"]
    assert result["inc"] and all(Path(p).is_absolute() for p in result["inc"])


def test_foreign_current_directory_module_not_loaded(tmp_path):
    (tmp_path / "ForeignProbe.pm").write_text('die "FOREIGN MODULE LOADED";\n')
    script = tmp_path / "script.pl"
    script.write_text("use ForeignProbe;\n")
    done = subprocess.run([PERL, "-I.", str(WRAPPER), str(script)], cwd=tmp_path,
                          capture_output=True, text=True)
    assert done.returncode != 0
    assert "Can't locate ForeignProbe.pm" in done.stderr
    assert "FOREIGN MODULE LOADED" not in done.stderr


@pytest.mark.parametrize("text,expected", [("exit 7;", 7), ("exit 0;", 0), ("undef;", 0),
                                          ('die "native failure";', None), ("this is not valid perl !!!", None)])
def test_native_exit_and_errors(tmp_path, text, expected):
    script = tmp_path / "script.pl"
    script.write_text(text)
    done = subprocess.run([PERL, str(WRAPPER), str(script)], cwd=tmp_path, capture_output=True)
    if expected is None:
        assert done.returncode != 0 and done.stderr
    else:
        assert done.returncode == expected


def test_relative_script_rejected(tmp_path):
    (tmp_path / "script.pl").write_text("1;")
    done = subprocess.run([PERL, str(WRAPPER), "script.pl"], cwd=tmp_path, capture_output=True)
    assert done.returncode != 0 and b"absolute" in done.stderr


@pytest.mark.skipif(os.environ.get("ORTHOHMM_LEGACY_BLAST_SMOKE") != "1", reason="Opt-in native Perl probe")
def test_installed_perl_probe_has_no_relative_lookup(tmp_path):
    from benchmark_tools.probe_orthomcl_bpo_parity import TOOL
    from benchmark_tools.run_qfo_corrected_blast import environment
    driver = WRAPPER.with_name("inspect_orthomcl_perl_runtime.pl")
    argv = [str(TOOL / "venv_orthomcl/bin/perl"), "-I" + str(TOOL), str(WRAPPER), str(driver)]
    done = subprocess.run(argv, cwd=tmp_path, env=environment(), capture_output=True, text=True, check=True)
    result = json.loads(done.stdout)
    assert "." not in result["search_path"]
    assert result["versions"]["orthomcl"] == "1.4"
