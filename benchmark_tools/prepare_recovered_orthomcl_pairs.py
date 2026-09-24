"""Convert admitted recovered native groups to unscored QfO group-clique pairs."""

import argparse
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.admit_blast_recovery_bpo import execution_identity
from benchmark_tools.admit_qfo_corrected_orthomcl import completed, frozen, unique_records
from benchmark_tools.audit_orthomcl_native_groups import audit
from benchmark_tools.orthomcl_matrix_to_pairwise import load_species
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.prepare_qfo_corrected_comparator_pairs import ENV_SHA
from benchmark_tools.prepare_qfo_corrected_orthomcl_pairs import write_pairs, SEMANTICS
from benchmark_tools.qfo_filter_pairs import load_mapping, filter_pairs
from benchmark_tools.run_blast_recovery_batch import save_status
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.verify_orthomcl_python_runtime import verify_runtime

ADMITTER_SHA = "3c02b1a701751cf98578107004a85b095f31aeb90041d17459f2caab20212821"
CONVERTER_SHA = "87b229469466cf330259db7f1ec5cd03609001bbb5d594ef49583bd4f832b23e"


def validate_admission(admission, job):
    if (admission["status"] != "recovered_orthomcl_native_outputs_admitted"
            or admission["conversion_authorized"] is not True
            or admission["admission_job_id"] != str(job) or admission["node"] != "bizon"
            or admission["allocated_cpus"] != 2 or admission["memory_mib"] != 65536
            or admission["accuracy_admitted"] is not False or admission["publication_ready"] is not False
            or admission["pair_semantics"] != SEMANTICS):
        raise ValueError("Require independently admitted recovered native groups")
    content = admission["content"]
    required = dict(State="COMPLETED", ExitCode="0:0", NodeList="bizon", AllocCPUS="180", ReqMem="900G")
    if (content["input_proteins"] != 984137 or content["input_species"] != 78
            or any(admission["scheduler"].get(k) != v for k, v in required.items())):
        raise ValueError("Wrong recovered native scope/allocation")
    for key in ("cross_species_clique_pairs", "final_groups", "grouped_proteins", "ungrouped_input_proteins"):
        if type(content[key]) is not int or content[key] < 0:
            raise ValueError("Invalid recovered group/pair counts")
    if content["grouped_proteins"] + content["ungrouped_input_proteins"] != 984137:
        raise ValueError("Inconsistent recovered group coverage")
    return content


def prepare(root, path, digest, job, executor, commit):
    identity, runtime = execution_identity(), verify_runtime(root)
    scheduler, accounting = completed(job, 2, "64G")
    expected = root / "benchmarks/results/qfo_blast_recovery_native_admission_v1/report.json"
    if path != expected or not executor.resolve().is_relative_to(root / "benchmarks/work"):
        raise ValueError("Wrong recovered native admission path/executor")
    frozen(executor, commit)
    admission = read_frozen(path, digest)
    content = validate_admission(admission, job)
    source = record(executor / "benchmark_tools/admit_recovered_orthomcl.py")
    converter = record(Path(__file__).with_name("prepare_qfo_corrected_orthomcl_pairs.py"))
    if (source["sha256"] != ADMITTER_SHA or source != admission["source"]
            or converter["sha256"] != CONVERTER_SHA):
        raise ValueError("Unreviewed native admission or pair converter source")
    native_scheduler, _ = completed(int(admission["scheduler"]["JobIDRaw"]), 180, "900G")
    if native_scheduler != admission["scheduler"]:
        raise ValueError("Native execution accounting changed")
    env_path = root / "benchmark_tools/results/qfo_assessment_environment_20260917.json"
    environment = read_frozen(env_path, ENV_SHA)
    mappings = [r for r in environment["reference_files"] if Path(r["path"]).name == "mapping.json.gz"]
    if len(mappings) != 1:
        raise ValueError("Require one frozen QfO mapping")
    mapping = mappings[0]
    checked = unique_records([record(path), source, converter, *admission["checked_records"],
        *admission["outputs"], admission["native_groups"], record(env_path), mapping,
        *[record(p) for p in sorted(Path(__file__).parent.glob("*.py"))]])
    for item in checked:
        check(item)
    base = root / "benchmarks/results/qfo_blast_recovery_native_v1"
    native = Path(admission["native_groups"]["path"]).parent
    if native.parent != base / "native_tool" or admission["native_groups"]["path"] != str(native / "all_orthomcl.out"):
        raise ValueError("Unexpected recovered native group path")
    gg = base / "inputs/all.gg"
    if record(gg) not in checked:
        raise ValueError("Missing admitted species mapping identity")
    owners = load_species(gg)
    if len(owners) != 984137 or len(set(owners.values())) != 78:
        raise ValueError("Wrong recovered species mapping scope")
    directory = root / "benchmarks/results/qfo_blast_recovery_pairs_v1"
    if directory.exists() or directory.is_symlink():
        raise FileExistsError(directory)
    directory.mkdir(parents=True, exist_ok=False)
    report = dict(status="recovered_pairs_preparing", source=record(__file__), method="orthomcl",
        participant="qfo_corrected_orthomcl_recovered", admission=record(path), admission_scheduler=scheduler,
        admission_accounting=accounting, mapping=mapping, checked_records=checked, runtime_before=runtime,
        query_coverage=admission["query_coverage"], accuracy_evaluated=False, publication_ready=False,
        job_id=identity["admission_job_id"], started_epoch=time.time(), semantics=SEMANTICS)
    save_status(directory / "results.json", report)
    try:
        validation = audit(native / "all_orthomcl.out", native / "tmp/all_ortho.mcl", native / "tmp/all_ortho.idx",
                           gg, directory / "groups.json")
        if validation["content"] != content:
            raise ValueError("Independent group audit differs from recovered admission")
        partial, mapped = directory / "pairs.partial.tsv", directory / "pairs.qfo.partial.tsv"
        with partial.open("x") as stream:
            counts = write_pairs(native / "all_orthomcl.out", owners, stream)
        expected = {key: content[key] for key in ("final_groups", "grouped_proteins", "ungrouped_input_proteins")}
        expected["total_pairs"] = content["cross_species_clique_pairs"]
        if counts != expected:
            raise ValueError("Recovered pair conversion counts differ")
        observed, retained = filter_pairs(partial, mapped, load_mapping(Path(mapping["path"])))
        if observed != retained or observed != counts["total_pairs"]:
            raise ValueError("Unexpected recovered mapping loss/count mismatch")
        for item in [*checked, *validation["checked_records"]]:
            check(item)
        runtime_after = verify_runtime(root)
        pairs, filtered = directory / "pairs.tsv", directory / "pairs.qfo.tsv"
        partial.rename(pairs)
        mapped.rename(filtered)
        report.update(status="recovered_orthomcl_pairs_prepared_unscored", content=counts,
            pairs=record(pairs), filtered_pairs=record(filtered), group_audit=record(directory / "groups.json"),
            total_pairs=observed, retained_pairs=retained, removed_mapping_pairs=0, native_duplicate_relations=0,
            runtime_after=runtime_after, limitations=[
                "Cross-species final-group cliques, not native phylogenetic ortholog pairs or graph edges.",
                "Ungrouped proteins and failed queries remain accounted for; no missing hits are repaired.",
                "Conversion remains unscored and requires terminal verification before assessment."])
    except BaseException as error:
        report.update(status="recovered_pair_conversion_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        report["finished_epoch"] = time.time()
        save_status(directory / "results.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("root", "admission", "executor"):
        parser.add_argument("--" + flag, type=Path, required=True)
    parser.add_argument("--admission-sha256", required=True)
    parser.add_argument("--job", type=int, required=True)
    parser.add_argument("--commit", required=True)
    args = parser.parse_args()
    print(prepare(args.root.resolve(), args.admission.resolve(), args.admission_sha256,
                  args.job, args.executor.resolve(), args.commit)["status"])
