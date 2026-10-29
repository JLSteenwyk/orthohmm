"""Keep the application checkout independent of the archived research suite."""

import ast
from pathlib import Path
import subprocess

from packaging.requirements import Requirement


ROOT = Path(__file__).resolve().parents[2]
RESEARCH = {"benchmark_tools", "manuscript", "benchmarks", "bacterial_scaling",
            "qfo_benchmark", "three_kingdoms"}


def test_tracked_application_tree_excludes_research_payload():
    paths = subprocess.check_output(["git", "ls-files", "-z"], cwd=ROOT, text=True).split("\0")
    assert not [path for path in paths if path and path.split("/", 1)[0] in RESEARCH]


def test_python_sources_do_not_import_research_drivers():
    for directory in ("orthohmm", "scripts", "tests"):
        for path in (ROOT / directory).rglob("*.py"):
            tree = ast.parse(path.read_text(), filename=str(path))
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    imports = [name.name for name in node.names]
                elif isinstance(node, ast.ImportFrom):
                    imports = [node.module or ""]
                else:
                    continue
                assert not [name for name in imports if name.split(".")[0] in RESEARCH], path


def test_runtime_manifest_matches_package_metadata():
    tree = ast.parse((ROOT / "setup.py").read_text())
    requires = next(ast.literal_eval(node.value) for node in tree.body
                    if isinstance(node, ast.Assign)
                    and any(isinstance(target, ast.Name) and target.id == "REQUIRES"
                            for target in node.targets))
    declared = [line.strip() for line in (ROOT / "requirements.txt").read_text().splitlines()
                if line.strip() and not line.lstrip().startswith("#")]
    assert {str(Requirement(line)) for line in declared} == {str(Requirement(line)) for line in requires}
    setup = next(node for node in tree.body if isinstance(node, ast.Expr)
                 and isinstance(node.value, ast.Call) and isinstance(node.value.func, ast.Name)
                 and node.value.func.id == "setup").value
    extras = ast.literal_eval(next(kw.value for kw in setup.keywords if kw.arg == "extras_require"))
    assert extras["phylogeny"] == ["dendropy>=5.0"]
    assert extras["hmmer"] == ["biopython>=1.85"]
