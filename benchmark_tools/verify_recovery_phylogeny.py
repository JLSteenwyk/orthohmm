"""Fresh recovery-install phylogeny fixture, retaining the historical template."""

import argparse
import copy
import json
from pathlib import Path
import subprocess

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.readback_canonical_ob_phylogeny import save

TEMPLATE_SHA = "03e0e517fe7158e85df5be2b60a7324d08c60d16911f546fd8e18abde7477213"


def fixture_plan(template, directory, python, runtime):
    if (template["fixture"] is not True or template["cpu"] != 1
            or template["attempts"] != 1 or template["checkpoint_reuse"] is not False
            or template["scoring"] is not False):
        raise ValueError("Require single-attempt, fresh, unscored fixture template")
    if runtime["executable"] != str(python):
        raise ValueError("Wrong recovery interpreter")
    site = python.parent.parent / "lib/python3.10/site-packages"
    for item in runtime["scientific_sources"] + runtime["dendropy_sources"]:
        if not Path(item["path"]).resolve().is_relative_to(site.resolve()):
            raise ValueError("Import escaped recovery installation")
    if not runtime["scientific_sources"] or not runtime["dendropy_sources"]:
        raise ValueError("Empty source inventory")
    plan = copy.deepcopy(template)
    plan.update(output=str(directory), inference=str(directory / "inference"),
                python=str(python), runtime=runtime)
    plan["checked_records"].extend(runtime["scientific_sources"] + runtime["dendropy_sources"])
    return plan


def verify(repo, directory):
    from benchmark_tools.audit_recovery_install import audit as install_audit
    from benchmark_tools.admit_canonical_ob_phylogeny import verify_artifacts
    from benchmark_tools import audit_phylogeny_structure as structure
    from benchmark_tools import audit_phylogeny_sequences as sequences
    from benchmark_tools import audit_phylogeny_events as events
    from benchmark_tools import audit_phylogeny_hierarchy as hierarchy

    if directory.exists():
        raise FileExistsError(directory)
    template_path = repo / "benchmarks/work/canonical_ob_phylogeny_smoke_20260926/plan.json"
    if record(template_path)["sha256"] != TEMPLATE_SHA:
        raise ValueError("Changed historical fixture template")
    template = json.loads(template_path.read_text())
    for item in [template["source"], *template["checked_records"]]:
        check(item)
    installation = repo / "benchmarks/work/publication_recovery_install_20260926"
    before = install_audit(repo, installation)
    python = installation / "venv/bin/python"
    code = ("import orthohmm,sys,json;sys.path.insert(0," + repr(str(repo)) + ");"
            "from benchmark_tools.run_canonical_ob_phylogeny import runtime_inventory;"
            "print(json.dumps(runtime_inventory()))")
    runtime = json.loads(subprocess.check_output([str(python), "-I", "-c", code],
                        env=template["environment"], text=True, timeout=60))
    plan = fixture_plan(template, directory, python, runtime)
    plan["checked_records"].extend([record(template_path), record(__file__)])
    directory.mkdir(parents=True)
    save(directory / "install_before.json", before)
    save(directory / "plan.json", plan)
    plan_record = record(directory / "plan.json")
    command = [str(python), "-I", plan["source"]["path"], "--run", plan_record["path"],
               "--plan-sha256", plan_record["sha256"]]
    with (directory / "native.log").open("x") as log:
        subprocess.run(command, env=plan["environment"], stdout=log, stderr=subprocess.STDOUT,
                       check=True, timeout=plan["timeout_seconds"])
    admission = verify_artifacts(directory, plan_record["sha256"], 4, 16)
    phylo = directory / "inference/orthohmm_phylogeny"
    reports = {"admission": admission}
    reports["structure"] = structure.audit(phylo, Path(plan["inputs"]))
    save(directory / "structure.json", reports["structure"])
    reports["sequences"] = sequences.audit(phylo, directory / "structure.json")
    reports["events"] = events.audit(phylo, directory / "structure.json", Path(plan["constraints"]["path"]))
    save(directory / "events.json", reports["events"])
    reports["hierarchy"] = hierarchy.audit(phylo, directory / "events.json")
    reports["install_after"] = install_audit(repo, installation)
    for name, report in reports.items():
        if name not in ("structure", "events"):
            save(directory / (name + ".json"), report)
    result = dict(status="fresh_recovery_phylogeny_fixture_verified", plan=plan_record,
                  command=command, source=record(__file__), summary=admission["summary"],
                  reports={n: record(directory / (n + ".json")) for n in reports},
                  full_dataset_reproduced=False, publication_ready=False,
                  limitations=["Same-host 16-gene fixture, not a full benchmark or portability test.",
                               "Historical template and external MAFFT/FastTree remain required."])
    save(directory / "result.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--directory", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(verify(args.repo.resolve(), args.directory.resolve()), indent=2))
