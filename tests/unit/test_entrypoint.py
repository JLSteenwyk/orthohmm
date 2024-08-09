import os
import subprocess

import pytest


class TestEntrypoint(object):
    @pytest.mark.slow
    def test_help(self):
        cmd = "orthohmm --help"
        exit_status = os.system(cmd)
        assert exit_status == 0

    @pytest.mark.slow
    def test_run(self):
        cmd = "orthohmm tests/samples/"
        exit_status = os.system(cmd)
        assert exit_status == 0

    def test_input_error(self):
        cmd = "orthohmm /file/doesnt/exist"
        response = subprocess.check_output([cmd], stderr=subprocess.STDOUT, shell=True)
        assert response == b"Input directory does not exist\n"

    def test_run_no_args(self):
        cmd = "orthohmm"
        exit_status = os.system(cmd)
        assert exit_status == 0
