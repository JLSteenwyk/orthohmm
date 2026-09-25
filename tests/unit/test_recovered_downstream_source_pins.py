import hashlib
import importlib
from pathlib import Path

import pytest


@pytest.mark.parametrize("consumer,constant,target", [
    ("admit_blast_recovery_bpo", "PREPARER_SHA", "prepare_blast_recovery_bpo"),
    ("recovered_orthomcl_inputs", "ADMITTER_SHA", "admit_blast_recovery_bpo"),
    ("admit_recovered_orthomcl", "RUNNER_SHA", "run_recovered_orthomcl"),
    ("prepare_recovered_orthomcl_pairs", "ADMITTER_SHA", "admit_recovered_orthomcl"),
    ("run_qfo_recovered_orthomcl_assessment", "CONVERTER_SHA", "prepare_recovered_orthomcl_pairs"),
    ("admit_qfo_recovered_orthomcl_assessment", "RUNNER_SHA", "run_qfo_recovered_orthomcl_assessment"),
])
def test_current_downstream_contract_matches_source(consumer, constant, target):
    module = importlib.import_module("benchmark_tools." + consumer)
    path = Path(module.__file__).with_name(target + ".py")
    assert hashlib.sha256(path.read_bytes()).hexdigest() == getattr(module, constant)
