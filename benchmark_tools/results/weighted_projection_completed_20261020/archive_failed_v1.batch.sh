#!/bin/sh
# This script was created by sbatch --wrap.

env CUDA_VISIBLE_DEVICES= OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 PYTHONHASHSEED=0 /home/bizon/anaconda3/bin/python -c 'from pathlib import Path
import csv, hashlib, io, json, os, shutil, subprocess
import xml.etree.ElementTree as ET
from benchmark_tools import plan_weighted_projection_followup as planner
from benchmark_tools import report_hmm_evidence_campaign as campaign
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.report_hmm_transfer_campaign import resource_table
from benchmark_tools.run_membership_hillclimb import write_json
from benchmark_tools.run_profile_weighting_followup import successful_job
from benchmark_tools.weighted_projection_experiment import point
repo=Path.cwd()
assert subprocess.check_output(["git","rev-parse","HEAD"],text=True).strip()=="fa34236e4667cba6fd9ba80f20a81b8ad5eff805"
report, plan, chosen, checked, provenance=planner.authenticate(repo)
assert report["advancing_labels"]==chosen["advancing_candidates"]==chosen["execution_arms"]==[]
assert chosen["stopped_candidates"]==["repair_henikoff","repair_henikoff_cluster"]
assert len(report["cells"])==6 and len(report["uncertainty"]["contrasts"])==6 and report["uncertainty"]["multiplicity_endpoints"]==18
launch_revision="c3b697d72516fc952ca84b1f0d6f03c863b9a83b"
blob=subprocess.check_output(["git","show",launch_revision+":"+campaign.LEDGER])
entry=json.loads(blob)["weighted_projection_followup_selection"]
assert entry["runner_commit"]=="70063b31a9ff5c62705b6f4f2ab06d743842ba55" and entry["job"]["job_id"]=="24540" and entry["job"]["dependency"]=="afterany:24539"
selection_job=successful_job("24540",1,"16G")
selection_pin=record(repo/entry["job"]["selection"])
selection=json.loads(Path(selection_pin["path"]).read_text())
assert selection["status"]=="no_weighted_projection_followup_required" and selection["selection"]==chosen and selection["datasets"]==[]
assert selection["pilot_readback"]==provenance["reader"]["reader"] and selection["pilot_plan"]==report["execution_plan"]
assert selection["planner_revision"]==entry["runner_commit"] and selection["defaults_promoted"] is False and selection["goal_completed"] is False
source=(repo/entry["source"]).resolve()
assert subprocess.check_output(["git","rev-parse","HEAD"],cwd=source,text=True).strip()==selection["planner_revision"]
expected_files=planner.sources(source)+[record(source/"benchmark_tools"/planner.experiment.PROTOCOL),record(source/"benchmark_tools/plan_weighted_projection_followup.slurm")]
assert expected_files==selection["planner_files"]
assert campaign.unique_pins([*checked,*expected_files])==selection["inputs"]
for pin in expected_files: check(pin)
test_root=repo/"benchmarks/results/hmm_campaign_final_regression_v1_20261020"
test_submission=record(test_root/"submission.json")
test=json.loads(Path(test_submission["path"]).read_text())
assert test["job_id"]=="24544" and test["source_revision"]=="fa34236e4667cba6fd9ba80f20a81b8ad5eff805" and test["modules"]==43 and test["expected_tests"]==1647
successful_job("24544",1,"16G")
for pin in [*test["source_files"],*test["test_files"],test["previous_junit"]]: check(pin)
suite=ET.parse(test_root/"tests.xml").getroot().find("testsuite")
assert suite is not None and suite.get("tests")=="1647" and all(suite.get(k)=="0" for k in ("errors","failures","skipped"))
log_path=test_root/"regression-24544.out"
affinity_lines=[line.split("=",1)[1] for line in log_path.read_text().splitlines() if line.startswith("TEST_CPU_AFFINITY=")]
assert len(affinity_lines)==1
affinity=json.loads(affinity_lines[0]); assert len(affinity)==1 and type(affinity[0]) is int
preservation={"status":hashlib.sha256(subprocess.check_output(["git","status","--porcelain=v1","--untracked-files=normal","--","tests/samples"])).hexdigest(),"diff":hashlib.sha256(subprocess.check_output(["git","diff","--binary","--","tests/samples"])).hexdigest()}
assert preservation==test["preservation_before"]
command=["sacct","-j","24538,24539,24540,24544","-X","-P","-o","JobIDRaw,State,ExitCode,Elapsed,AllocCPUS,ReqMem,NodeList"]
raw=subprocess.check_output(command,text=True); jobs=list(csv.DictReader(io.StringIO(raw),delimiter="|"))
assert len(jobs)==4 and {r["JobIDRaw"] for r in jobs}=={"24538","24539","24540","24544"}
for job in jobs:
 assert job["State"]=="COMPLETED" and job["ExitCode"]=="0:0" and job["NodeList"]=="bizon"
 assert (job["AllocCPUS"],job["ReqMem"])==(("16","128G") if job["JobIDRaw"]=="24538" else ("1","16G"))
output=repo/"benchmarks/results/weighted_projection_completed_archive_v1_20261020"
assert not output.exists(); output.mkdir()
check(report["score_table"])
shutil.copyfile(report["score_table"]["path"],output/"scores.tsv")
assert record(output/"scores.tsv")["sha256"]==report["score_table"]["sha256"]
methods=[{"key":"weighted_projection:"+cell["label"],"cells":[cell]} for cell in report["cells"]]
costs,cost_pins=resource_table(output,methods)
assert len(costs)==24
for pin in cost_pins: check(pin)
individual=record(repo/"benchmark_tools/results/hmm_individual_effects_20261019/summary.json")
assert individual["sha256"]=="37630ffea7abf8f8fb087c4a0db2f14ab3a3fdfefbd0f77a4caf9fb8924c7c61"
cells=[{key:cell[key] for key in ("label","dataset","status","source_revision","arm_parameters","coverage","activation","candidate","products","rooted_species_tree","graph_stages","native_pair_turnover")} | {"point":point(cell),"raw_checked_records":len(cell["checked"])} for cell in report["cells"]]
result={"schema_version":1,"status":"complete_six_arm_interaction_admitted_no_combination_advances","archive_job_id":os.environ["SLURM_JOB_ID"],"archive_source_revision":"fa34236e4667cba6fd9ba80f20a81b8ad5eff805","scientific_revision":planner.SCIENCE,"source_admission":provenance,"source_report":provenance["reader"]["reader"],"scientific_plan":report["execution_plan"],"scientific_producer":report["producer"],"planner_submission":{"revision":launch_revision,"path":campaign.LEDGER,"git_blob_sha256":hashlib.sha256(blob).hexdigest()},"selection_source":selection_pin,"selection":chosen,"cells":cells,"matched_gates":report["matched_gates"],"uncertainty":report["uncertainty"],"source_content_records_rechecked":len(checked),"planner_source_records_rechecked":len(expected_files),"all_individual_effects":individual,"score_table":record(output/"scores.tsv"),"original_score_table":report["score_table"],"resource_table":record(output/"resources.tsv"),"resource_rows":len(costs),"scheduler":{"command":command,"raw":raw,"jobs":jobs},"tests":{"submission":test_submission,"junit":record(test_root/"tests.xml"),"log":record(log_path),"passed":1647,"modules":43,"errors":0,"failures":0,"skipped":0,"cpu_affinity":affinity,"source_records_rechecked":len(test["source_files"]),"test_source_records_rechecked":len(test["test_files"])},"preservation":preservation,"interpretation":"Both combinations improve F1 over canonical and repair solo but score below their corresponding weighting-only constituent, so both fixed stronger-constituent gates fail. Preserve the bounded negative interaction result without retuning or QfO/Three Kingdoms extension; no statistical inferiority or independent generalization claim.","scope":"Complete inherited separate frozen source/native/official-score admission, with original immutable launches, producer/plan, all source content fingerprints and selection gates freshly authenticated. No new biological inference, scoring or functional sampling. Four controls reproduced exact original graph/candidate/native/rooted-tree products and P/R/F. All15901 actual profiles in each candidate used the declared joint policies. This archive closes only the interaction branch, not the overall goal.","resource_scope":"Shared-host cached incremental phases, with different constituent/candidate cache reuse. Native child CPU/GNU-time maximum process RSS are not simultaneous aggregate memory. Candidate/arm totals are retained producer metadata, not independently retimed; do not sum totals and components or infer end-to-end tool speed ratios.","required_next_actions":["Freeze final retained experimental candidate scope before post-freeze confirmation or documented data limits.","Complete the updated all-method comparison/recommendation with endpoint-specific precision/recall, coverage, assessed populations and costs.","Finish requirement-by-requirement completion audit and verified commits/pushes; defaults remain unchanged unless justified."],"accuracy_result_admitted":True,"defaults_promoted":False,"goal_completed":False}
write_json(output/"summary.json",result)
print(json.dumps({"status":result["status"],"summary":record(output/"summary.json"),"tests_passed":1647,"resource_rows":24,"advancing_candidates":[]},indent=2))
'
