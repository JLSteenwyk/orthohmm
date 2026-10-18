from dataclasses import replace
import json
from pathlib import Path
import shutil
import subprocess

from Bio import AlignIO
import numpy as np
import pytest

from orthohmm.search import selective_hmmer as hmmer
from orthohmm.search.matrices import get_background_freqs, get_matrix
from orthohmm.search.msa_profile import build_msa_profile


def source_fixture():
    rng = np.random.default_rng(18102026)
    alphabet = np.array(list("ACDEFGHIKLMNPQRSTVWY"))
    ancestor = rng.choice(alphabet, 200)
    members = {}
    for name in ("z_source|1", "a_source|2", "m_source|3"):
        sequence = ancestor.copy()
        positions = rng.choice(200, 5, replace=False)
        sequence[positions] = rng.choice(alphabet, len(positions))
        members[name] = "".join(sequence)
    profile = build_msa_profile(list(members.values()), list(members), get_matrix("BLOSUM62"), get_background_freqs("BLOSUM62"))
    query = {"full|candidate": members["z_source|1"], "fragment|candidate": members["a_source|2"][40:130],
             "tandem|candidate": members["z_source|1"] + "G" * 20 + members["z_source|1"],
             "unknown|candidate": "X" * 200, "unrelated|candidate": "".join(rng.choice(alphabet, 200))}
    return profile, members, query


@pytest.fixture(scope="module")
def tools():
    paths = {role: shutil.which(role) for role in ("hmmbuild", "hmmsearch")}
    if not all(paths.values()):
        pytest.skip("Requires real HMMER3.4 tools")
    for role, path in paths.items():
        output = subprocess.check_output([path, "-h"], text=True)
        if "# HMMER 3.4 " not in output:
            pytest.skip("Requires the frozen HMMER3.4 backend, not another version")
    return paths


@pytest.fixture(scope="module")
def scored(tmp_path_factory, tools):
    directory = tmp_path_factory.mktemp("selective_hmmer") / "score"
    profile, members, query = source_fixture()
    result = hmmer.score_family_hmmer(profile, members, query, 10000, directory, **tools)
    return directory, result, profile, members, query


def test_real_hmmer_control_preserves_source_columns_members_and_full_database_size(scored):
    directory, result, profile, members, query = scored
    assert result["status"] == "hmmer_profile_scored"
    assert not result["orthology_inferred"] and not result["edges_changed"]
    assert result["calibration_complete"] and not result["missing_members"]
    assert result["model_length"] == profile.length == 200
    assert float(result["model_header"]["EFFN"]) == 3
    assert result["parameters"]["database_genes"] == 10000
    assert result["parameters"]["searched_distinct_targets"] == 8
    build, search = result["commands"]["build"], result["commands"]["search"]
    assert all(flag in build for flag in ("--hand", "--amino", "--wnone", "--enone"))
    assert build[build.index("--cpu") + 1] == search[search.index("--cpu") + 1] == "0"
    assert search[search.index("-Z") + 1] == "10000"
    assert "--max" in search and "--nonull2" not in search
    original = AlignIO.read(directory / "source.sto", "stockholm")
    assert original.column_annotations["reference_annotation"].count("x") == profile.length
    assert [result["token_to_gene"][row.id] for row in original] == list(members)
    assert all(row["reported"] for gene, row in result["measurements"].items() if gene in members)
    full = result["measurements"]["full|candidate"]
    assert full["reported"] and full["sequence_evalue"] < 1e-4
    assert full["full_sequence_score"] >= result["minimum_member_score"]
    assert full["domains"] and full["domains"][0]["coordinate_system"] == "zero_based_half_open"
    assert result["measurements"]["fragment|candidate"]["reported"]
    assert set(result["measurements"]) == set(members) | set(query)
    assert json.loads((directory / "measurements.json").read_text()) == result
    assert all(hmmer._record(pin["path"]) == pin for pin in result["files"])
    with pytest.raises(FileExistsError):
        hmmer.score_family_hmmer(profile, members, query, 10000, directory)


def test_real_hmmer_full_database_scaling_changes_evalues_not_bit_scores(scored, tmp_path, tools):
    _, old, profile, members, query = scored
    new = hmmer.score_family_hmmer(profile, members, query, 1000000, tmp_path / "larger", **tools)
    for gene, before in old["measurements"].items():
        after = new["measurements"][gene]
        assert after["reported"] == before["reported"]
        if before["reported"]:
            assert after["full_sequence_score"] == before["full_sequence_score"]
            assert after["sequence_evalue"] == pytest.approx(100 * before["sequence_evalue"], rel=.06, abs=0)


def test_real_domain_coordinates_use_model_and_target_not_swapped(scored):
    _, result, _, _, query = scored
    fragment = result["measurements"]["fragment|candidate"]
    domain = fragment["domains"][0]
    assert 25 <= domain["model_start"] <= 60
    assert domain["model_end"] > 100
    assert 0 <= domain["target_start"] < domain["target_end"] <= len(query["fragment|candidate"])
    assert domain["envelope_start"] <= domain["target_start"] < domain["target_end"] <= domain["envelope_end"]
    tandem = result["measurements"]["tandem|candidate"]
    assert len(tandem["domains"]) == tandem["domain_counts"]["domain_reported_num"] == 2
    assert [row["index"] for row in tandem["domains"]] == [1, 2]


def table_copy(scored, tmp_path):
    directory, result, profile, members, query = scored
    sequence, domain = tmp_path / "sequence.tbl", tmp_path / "domain.tbl"
    shutil.copyfile(directory / "sequence.tbl", sequence)
    shutil.copyfile(directory / "domain.tbl", domain)
    lengths = {token: len({**members, **query}[gene]) for token, gene in result["token_to_gene"].items()}
    return sequence, domain, result["token_to_gene"], lengths, int(profile.length)


def mutate_first_data(path, column, value):
    lines = path.read_text().splitlines()
    index = next(index for index, line in enumerate(lines) if not line.startswith("#") and line.strip())
    fields = lines[index].split()
    fields[column] = value
    lines[index] = " ".join(fields)
    path.write_text("\n".join(lines) + "\n")


@pytest.mark.parametrize("kind", ["footer", "version", "program", "model", "token", "nan", "negative_evalue", "duplicates", "interrupted"])
def test_bad_or_incomplete_sequence_tables_are_failures_not_missing_hits(scored, tmp_path, kind):
    args = table_copy(scored, tmp_path)
    sequence = args[0]
    if kind == "footer":
        sequence.write_text(sequence.read_text().replace("# [ok]", "# interrupted"))
    elif kind == "version":
        sequence.write_text(sequence.read_text().replace("# Version:", "# old-version:"))
    elif kind == "program":
        sequence.write_text(sequence.read_text().replace("# Program:", "# old-program:"))
    elif kind in ("model", "token", "nan", "negative_evalue"):
        column, value = {"model": (2, "WrongProfile"), "token": (0, "g99999999"),
                         "nan": (5, "nan"), "negative_evalue": (4, "-1")}[kind]
        mutate_first_data(sequence, column, value)
    else:
        lines = sequence.read_text().splitlines()
        index = next(i for i, line in enumerate(lines) if not line.startswith("#") and line.strip())
        lines.insert(index if kind == "duplicates" else index + 1,
                     lines[index] if kind == "duplicates" else "# premature footer")
        sequence.write_text("\n".join(lines) + "\n")
    with pytest.raises(ValueError):
        hmmer.read_hmmer_tables(*args)


@pytest.mark.parametrize("column,value", [(2, "999"), (5, "999"), (9, "0"), (10, "999"), (15, "0"),
    (18, "99999"), (19, "999"), (21, "nan"), (21, "1.5"), (13, "inf"), (11, "-1")])
def test_invalid_domain_identity_lengths_coordinates_and_values_rejected(scored, tmp_path, column, value):
    args = table_copy(scored, tmp_path)
    mutate_first_data(args[1], column, value)
    with pytest.raises(ValueError):
        hmmer.read_hmmer_tables(*args)


def test_genuinely_unreported_target_remains_explicitly_missing(scored, tmp_path):
    args = list(table_copy(scored, tmp_path))
    token = next(token for token, gene in args[2].items() if gene == "full|candidate")
    for path in args[:2]:
        path.write_text("\n".join(line for line in path.read_text().splitlines()
                        if line.startswith("#") or not line.strip() or line.split()[0] != token) + "\n")
    result = hmmer.read_hmmer_tables(*args)
    assert result["full|candidate"] == {"reported": False, "full_sequence_score": None, "sequence_evalue": None, "domains": []}
    assert any(row["reported"] for row in result.values())


def test_real_unreported_source_member_cannot_be_silently_removed_from_calibration(tmp_path, tools):
    _, members, query = source_fixture()
    members["m_source|3"] = "X" * 200
    profile = build_msa_profile(list(members.values()), list(members), get_matrix("BLOSUM62"), get_background_freqs("BLOSUM62"))
    result = hmmer.score_family_hmmer(profile, members, query, 10000, tmp_path / "incomplete", **tools)
    assert result["status"] == "hmmer_profile_scored" and not result["calibration_complete"]
    assert "m_source|3" in result["missing_members"]
    assert result["minimum_member_score"] is None
    assert not result["measurements"]["m_source|3"]["reported"]
    assert result["source_alignment_unknowns"]["m_source|3"]["encoded_unknown_residues"] == 200
    assert result["source_alignment_unknowns"]["m_source|3"]["aligned_ungapped_length"] == 0


def test_overlapping_query_and_member_are_one_target_with_both_roles(tmp_path, tools):
    profile, members, _ = source_fixture()
    query = {"z_source|1": members["z_source|1"]}
    result = hmmer.score_family_hmmer(profile, members, query, 10000, tmp_path / "deduplicated", **tools)
    assert result["parameters"]["searched_distinct_targets"] == 3
    assert set(result["measurements"]) == set(members)
    assert "z_source|1" in result["queries"] and "z_source|1" in result["members"]


def test_scoring_preserves_readonly_native_parameters_and_source_order(tmp_path, tools):
    profile, members, query = source_fixture()
    arrays = [profile.match_emissions, profile.insert_emissions, profile.transitions, profile.consensus]
    before = [array.tobytes() for array in arrays]
    sources_before, queries_before, ids_before = dict(members), dict(query), list(profile.member_ids)
    for array in arrays:
        array.setflags(write=False)
    result = hmmer.score_family_hmmer(profile, members, query, 10000, tmp_path / "readonly", **tools)
    assert result["status"] == "hmmer_profile_scored"
    assert [array.tobytes() for array in arrays] == before
    assert members == sources_before and query == queries_before and profile.member_ids == ids_before


@pytest.mark.parametrize("kind", ["emissions", "consensus", "transitions", "weighted", "source_ownership", "query_ownership", "database_count"])
def test_invalid_native_source_or_scope_fails_before_external_tools(tmp_path, monkeypatch, kind):
    profile, members, query = source_fixture()
    database_count = 10000
    if kind == "emissions":
        profile.match_emissions[0, 0] += 1
    elif kind == "consensus":
        profile.consensus[0] = (profile.consensus[0] + 1) % 20
    elif kind == "transitions":
        profile.transitions[1] = -13
    elif kind == "weighted":
        profile = replace(profile, sequence_weighting="henikoff")
    elif kind == "source_ownership":
        members.pop(next(iter(members)))
    elif kind == "query_ownership":
        query[next(iter(members))] = "A" * 200
    else:
        database_count = 1
    def forbidden(*args, **kwargs):
        raise AssertionError("Invalid source must not run external scoring")
    monkeypatch.setattr(hmmer, "_run_tool", forbidden)
    directory = tmp_path / "failed"
    with pytest.raises(ValueError):
        hmmer.score_family_hmmer(profile, members, query, database_count, directory)
    receipt = json.loads((directory / "measurements.json").read_text())
    assert receipt["status"] == "failed" and receipt["error_type"] == "ValueError"
    assert not receipt["edges_changed"] and not receipt["orthology_inferred"]


def executable(path, source):
    path.write_text("#!/usr/bin/env python3\n" + source)
    path.chmod(0o755)
    return str(path)


def test_wrong_backend_version_is_rejected_and_retained(tmp_path):
    tool = executable(tmp_path / "wrong", "print('# hmmbuild :: wrong version')\nprint('# HMMER 3.5 (future)')\n")
    profile, members, query = source_fixture()
    with pytest.raises(ValueError, match="HMMER3.4"):
        hmmer.score_family_hmmer(profile, members, query, 10000, tmp_path / "failed", hmmbuild=tool)
    result = json.loads((tmp_path / "failed/measurements.json").read_text())
    assert result["status"] == "failed" and result["tools"]["hmmbuild"]["path"] == str(Path(tool).resolve())
    assert any(Path(pin["path"]).name == "hmmbuild_help.stdout.txt" for pin in result["files"])


def test_failed_build_preserves_exact_command_stderr_and_terminal_receipt(tmp_path, tools):
    tool = executable(tmp_path / "build_failure", "import sys\nif '-h' in sys.argv:\n print('# hmmbuild :: fixture')\n print('# HMMER 3.4 (fixture)')\nelse:\n print('controlled build failure', file=sys.stderr)\n sys.exit(7)\n")
    profile, members, query = source_fixture()
    with pytest.raises(RuntimeError, match="exited 7"):
        hmmer.score_family_hmmer(profile, members, query, 10000, tmp_path / "failed", hmmbuild=tool, hmmsearch=tools["hmmsearch"])
    result = json.loads((tmp_path / "failed/measurements.json").read_text())
    assert result["status"] == "failed" and result["error_type"] == "RuntimeError"
    assert "--hand" in result["commands"]["build"]
    assert (tmp_path / "failed/build.stderr.txt").read_text().strip() == "controlled build failure"
    assert not result["edges_changed"]


def test_tool_timeout_is_a_failure_not_an_unreported_hit(tmp_path, monkeypatch):
    def timeout(command, **kwargs):
        raise subprocess.TimeoutExpired(command, kwargs["timeout"])
    monkeypatch.setattr(hmmer.subprocess, "run", timeout)
    profile, members, query = source_fixture()
    with pytest.raises(subprocess.TimeoutExpired):
        hmmer.score_family_hmmer(profile, members, query, 10000, tmp_path / "failed", hmmbuild=__file__)
    result = json.loads((tmp_path / "failed/measurements.json").read_text())
    assert result["status"] == "failed" and result["error_type"] == "TimeoutExpired"
    assert not result.get("measurements")
