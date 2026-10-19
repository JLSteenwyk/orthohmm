"""Controlled, opt-in HMMER measurements; no edge or orthology decisions.

Uses the original native source alignment and match-column mask, uniform
relative weights and actual member count. This is an explicit HMMER control,
not the tool's default weighting or a numerically identical native model.
"""

import json
from io import StringIO
import math
from numbers import Integral
import os
from pathlib import Path
import shutil
import subprocess
import time

import numpy as np
from Bio import AlignIO, SearchIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ..phylogeny_pipeline import _atomic_write, _sha256_file
from .matrices import ALPHABET, get_background_freqs, get_matrix
from .msa_center_star import center_star_msa
from .msa_profile import compute_pssm, encode_alignment
from .sequences import encode_sequence


MODEL_NAME = "SelectiveProfile"


def _record(path):
    path = Path(path).resolve()
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": _sha256_file(path)}


def _normalized_sequence(sequence):
    if not isinstance(sequence, str) or not sequence:
        raise ValueError("Require nonempty protein sequence strings")
    return "".join(ALPHABET[int(aa)] if aa < 20 else "X" for aa in encode_sequence(sequence))


def _run_tool(command, directory, label, environment, timeout):
    before = time.monotonic()
    with (directory / (label + ".stdout.txt")).open("x") as out, (directory / (label + ".stderr.txt")).open("x") as err:
        completed = subprocess.run(command, stdin=subprocess.DEVNULL, stdout=out, stderr=err,
                                   env=environment, timeout=timeout, check=False)
    if completed.returncode:
        raise RuntimeError(f"HMMER {label} exited {completed.returncode}; retained stderr: {directory / (label + '.stderr.txt')}")
    return {"command": command, "returncode": completed.returncode, "wall_seconds": time.monotonic() - before}


def _checked_table(path):
    """Check completion and row boundaries before the structured parser runs."""
    lines = path.read_text().splitlines()
    if not lines or lines[-1].strip() != "# [ok]":
        raise ValueError("Incomplete HMMER table footer: " + str(path))
    if not any(line.startswith("# Program:") and line.split(":", 1)[1].strip() == "hmmsearch" for line in lines):
        raise ValueError("Unexpected HMMER table program")
    if not any(line.startswith("# Version:") and line.split(":", 1)[1].strip().startswith("3.4 ") for line in lines):
        raise ValueError("Unexpected HMMER table version")
    rows, footer = [], False
    for line in lines:
        if not line.strip() or line.startswith("#"):
            footer |= bool(rows)
        else:
            if footer:
                raise ValueError("HMMER data row after an interrupted table body")
            rows.append(line)
    return rows


def _parsed_rows(rows, format_name):
    # Grouped SearchIO output retains only the last row's repeated metadata.
    # Parsing each row also checks lengths/identities in earlier domain rows.
    results = []
    for row in rows:
        parsed = list(SearchIO.parse(StringIO(row + "\n"), format_name))
        if len(parsed) != 1 or len(parsed[0]) != 1:
            raise ValueError("Unexpected HMMER single-row structure")
        results.append(parsed[0])
    return results


def _finite(value, name, *, nonnegative=False):
    value = float(value)
    if not math.isfinite(value) or (nonnegative and value < 0):
        raise ValueError("Invalid HMMER " + name)
    return value


def read_hmmer_tables(sequence_table, domain_table, token_to_gene, target_lengths, model_length):
    """Read complete HMMER3.4 tables using Biopython's HMMER parsers.

    A missing reported target remains missing; it is not assigned a synthetic
    score. Coordinates are normalized by SearchIO to zero-based half-open.
    """
    sequence_table, domain_table = Path(sequence_table), Path(domain_table)
    if (set(target_lengths) != set(token_to_gene) or len(set(token_to_gene.values())) != len(token_to_gene)
            or any(not isinstance(gene, str) or not gene for gene in token_to_gene.values())
            or type(model_length) is not int or model_length <= 0
            or any(type(length) is not int or length <= 0 for length in target_lengths.values())):
        raise ValueError("Invalid HMMER model/target identity or lengths")
    sequence_rows, domain_rows = _checked_table(sequence_table), _checked_table(domain_table)
    sequences = _parsed_rows(sequence_rows, "hmmer3-tab")
    domains = _parsed_rows(domain_rows, "hmmsearch3-domtab")
    if any(result.id != MODEL_NAME for result in (*sequences, *domains)):
        raise ValueError("Unexpected HMMER model identity")
    measured = {}
    for result in sequences:
        hit = result[0]
        if hit.id not in token_to_gene or hit.id in measured or len(hit.hsps) != 1:
            raise ValueError("Unexpected or duplicate HMMER target identity")
        counts = {key: getattr(hit, key) for key in ("region_num", "cluster_num", "overlap_num", "env_num",
                  "domain_obs_num", "domain_reported_num", "domain_included_num")}
        if any(type(value) is not int or value < 0 for value in counts.values()):
            raise ValueError("Invalid HMMER domain counts")
        if not counts["domain_included_num"] <= counts["domain_reported_num"] <= counts["domain_obs_num"]:
            raise ValueError("Inconsistent HMMER observed/reported/included domain counts")
        best = hit.hsps[0]
        measured[hit.id] = {
            "reported": True, "full_sequence_score": _finite(hit.bitscore, "sequence bit score"),
            "sequence_evalue": _finite(hit.evalue, "sequence E-value", nonnegative=True),
            "sequence_bias": _finite(hit.bias, "sequence bias"),
            "best_domain_score": _finite(best.bitscore, "best-domain bit score"),
            "best_domain_evalue": _finite(best.evalue, "best-domain E-value", nonnegative=True),
            "best_domain_bias": _finite(best.bias, "best-domain bias"),
            "expected_domains": _finite(hit.domain_exp_num, "expected domains", nonnegative=True),
            "domain_counts": counts, "domains": [],
        }
    if len(measured) != len(sequence_rows):
        raise ValueError("Parsed HMMER sequence-row count differs")
    for result, raw in zip(domains, domain_rows):
        if result.seq_len != model_length:
            raise ValueError("HMMER domain table model length differs")
        for hit in result:
            if hit.id not in measured or hit.seq_len != target_lengths[hit.id] or len(hit.hsps) != 1:
                raise ValueError("HMMER domain target identity or length differs")
            summary = measured[hit.id]
            # SearchIO omits the redundant per-row "of" count.
            if int(raw.split()[10]) != summary["domain_counts"]["domain_reported_num"]:
                raise ValueError("HMMER domain total differs from sequence summary")
            for key, actual in (("full_sequence_score", hit.bitscore), ("sequence_evalue", hit.evalue), ("sequence_bias", hit.bias)):
                if _finite(actual, key) != summary[key]:
                    raise ValueError("HMMER sequence/domain summaries differ")
            for hsp in hit.hsps:
                if not (0 <= hsp.query_start < hsp.query_end <= model_length
                        and 0 <= hsp.env_start <= hsp.hit_start < hsp.hit_end <= hsp.env_end <= hit.seq_len):
                    raise ValueError("Invalid HMMER domain alignment/envelope coordinates")
                accuracy = _finite(hsp.acc_avg, "domain posterior accuracy")
                if not 0 <= accuracy <= 1:
                    raise ValueError("Invalid HMMER domain posterior accuracy")
                summary["domains"].append({
                    "index": hsp.domain_index, "model_start": hsp.query_start, "model_end": hsp.query_end,
                    "target_start": hsp.hit_start, "target_end": hsp.hit_end,
                    "envelope_start": hsp.env_start, "envelope_end": hsp.env_end,
                    "score": _finite(hsp.bitscore, "domain bit score"), "bias": _finite(hsp.bias, "domain bias"),
                    "independent_evalue": _finite(hsp.evalue, "domain independent E-value", nonnegative=True),
                    "conditional_evalue": _finite(hsp.evalue_cond, "domain conditional E-value", nonnegative=True),
                    "posterior_accuracy": accuracy, "coordinate_system": "zero_based_half_open",
                })
    if (sum(len(row["domains"]) for row in measured.values()) != len(domain_rows)
            or any(len(row["domains"]) != row["domain_counts"]["domain_reported_num"] for row in measured.values())):
        raise ValueError("Parsed HMMER domain-row count differs")
    if any(sorted(domain["index"] for domain in row["domains"]) != list(range(1, len(row["domains"]) + 1))
           for row in measured.values()):
        raise ValueError("HMMER domain indices are duplicated or incomplete")
    return {gene: measured.get(token, {"reported": False, "full_sequence_score": None,
                                      "sequence_evalue": None, "domains": []})
            for token, gene in token_to_gene.items()}


def score_family_hmmer(profile, member_sequences, query_sequences, full_database_genes,
                      output_directory, *, hmmbuild="hmmbuild", hmmsearch="hmmsearch"):
    """Measure one original family and its selected query targets with HMMER.

    Every source member participates in calibration. Failure and missing output
    remain distinct; neither this routine nor its tables establish orthology.
    """
    directory = Path(output_directory)
    directory.mkdir(parents=True, exist_ok=False)
    receipt = directory / "measurements.json"
    result = {"schema_version": 1, "status": "running", "commands": {}, "tools": {},
              "orthology_inferred": False, "edges_changed": False, "files": []}
    _atomic_write(receipt, json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
    try:
        members = list(profile.member_ids)
        if (not 3 <= len(members) <= 200 or len(set(members)) != len(members)
                or any(not isinstance(gene, str) or not gene for gene in (*members, *query_sequences))
                or set(member_sequences) != set(members) or not 1 <= len(query_sequences) <= 128
                or profile.sequence_weighting != "none"):
            raise ValueError("Require an original unweighted family and bounded unique query set")
        sources = {gene: _normalized_sequence(member_sequences[gene]) for gene in members}
        queries = {gene: _normalized_sequence(sequence) for gene, sequence in query_sequences.items()}
        if any(sources[gene] != queries[gene] for gene in set(sources) & set(queries)):
            raise ValueError("Query/source sequence ownership differs")
        targets = {**sources, **queries}
        if (isinstance(full_database_genes, (bool, np.bool_)) or not isinstance(full_database_genes, Integral)
                or not len(targets) <= full_database_genes <= np.iinfo(np.int32).max):
            raise ValueError("Full original database gene count must cover all distinct targets")
        matrix, background = get_matrix("BLOSUM62"), get_background_freqs("BLOSUM62")
        aligned = center_star_msa([sources[gene] for gene in members], matrix,
                                  gap_open=-12, gap_extend=-3, n_threads=1,
                                  **({"projection_policy": profile.projection_policy}
                                     if profile.projection_policy != "legacy" else {}))
        if aligned is None:
            raise ValueError("Native source alignment could not be reconstructed")
        aligned_residues = {gene: sequence.replace("-", "").replace(".", "").upper()
                            for gene, sequence in zip(members, aligned)}
        if any(aligned_residues[gene].replace("X", "") != sources[gene].replace("X", "") for gene in members):
            raise ValueError("Native alignment changed known source residues")
        emissions, mask, consensus = compute_pssm(encode_alignment(aligned), matrix, background)
        transitions = np.array([0, -12, -12, -1, -3, -1, -3], dtype=np.int32)
        if (profile.length <= 0 or profile.length != len(emissions)
                or not np.array_equal(emissions, profile.match_emissions)
                or not np.array_equal(consensus, profile.consensus)
                or not np.array_equal(profile.transitions, transitions)
                or not np.array_equal(profile.insert_emissions, np.full(20, -1, dtype=np.int8))):
            raise ValueError("Reconstructed native profile differs from original parameters")
        token_to_gene = {f"g{index:08d}": gene for index, gene in enumerate(sorted(targets))}
        gene_to_token = {gene: token for token, gene in token_to_gene.items()}
        alignment = MultipleSeqAlignment([SeqRecord(Seq(sequence), id=gene_to_token[gene], description="")
                                          for gene, sequence in zip(members, aligned)])
        alignment.column_annotations["reference_annotation"] = "".join("x" if match else "." for match in mask)
        AlignIO.write(alignment, directory / "source.sto", "stockholm")
        SeqIO.write([SeqRecord(Seq(targets[gene]), id=token, description="") for token, gene in token_to_gene.items()],
                    directory / "targets.faa", "fasta")
        result.update({"model_name": MODEL_NAME, "members": members, "queries": sorted(queries),
            "token_to_gene": token_to_gene, "model_length": int(profile.length),
            "source_alignment_unknowns": {gene: {"encoded_unknown_residues": sources[gene].count("X"),
                "aligned_unknown_residues": aligned_residues[gene].count("X"),
                "encoded_length": len(sources[gene]), "aligned_ungapped_length": len(aligned_residues[gene])}
                for gene in members},
            "source_alignment_scope": "Preserves the existing native representation, where unknown residues may become gaps. Calibration FASTA retains original native-encoded residues; this is not a correction of that behavior.",
            "parameters": {"database_genes": int(full_database_genes), "searched_distinct_targets": len(targets),
                "sequence_normalization": "native_20_amino_acids_else_X", "sequence_weighting": "none",
                "effective_sequence_number": len(members), "seed": 42, "worker_threads": 0,
                "per_tool_timeout_seconds": 600, "score_precision": "HMMER CLI reported values; ties may reflect rounding"}})
        if profile.projection_policy != "legacy":
            result["parameters"]["profile_projection_policy"] = profile.projection_policy
            result["source_alignment_scope"] = (
                "Residue-preserving native projection; every normalized source residue, including X, is retained."
            )
        environment = {**os.environ, "LC_ALL": "C", "HMMER_NCPU": "0", "OMP_NUM_THREADS": "1",
                       "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
        result["environment_overrides"] = {key: environment[key] for key in (
            "LC_ALL", "HMMER_NCPU", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}
        for role, executable in (("hmmbuild", hmmbuild), ("hmmsearch", hmmsearch)):
            binary = Path(shutil.which(str(executable)) or str(executable)).resolve(strict=True)
            result["tools"][role] = _record(binary)
            help_command = [str(binary), "-h"]
            result["commands"][role + "_help"] = help_command
            _run_tool(help_command, directory, role + "_help", environment, 10)
            help_text = (directory / (role + "_help.stdout.txt")).read_text()
            if (not help_text.startswith("# " + role + " ::")
                    or not any(line.startswith("# HMMER 3.4 ") for line in help_text.splitlines())):
                raise ValueError("Require the explicit HMMER3.4 " + role + " backend")
        build = [result["tools"]["hmmbuild"]["path"], "--amino", "--hand", "--wnone", "--enone",
                 "--seed", "42", "--cpu", "0", "-n", MODEL_NAME, "-O", str(directory / "hmmer_source.sto"),
                 str(directory / "profile.hmm"), str(directory / "source.sto")]
        search = [result["tools"]["hmmsearch"]["path"], "--max", "--seed", "42", "--cpu", "0", "--noali",
                  "-T", "-1000000000", "--domT", "-1000000000", "-Z", str(int(full_database_genes)),
                  "--tblout", str(directory / "sequence.tbl"), "--domtblout", str(directory / "domain.tbl"),
                  str(directory / "profile.hmm"), str(directory / "targets.faa")]
        result["commands"].update(build=build, search=search)
        _atomic_write(receipt, json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
        result["build_execution"] = _run_tool(build, directory, "build", environment, 600)
        annotated = AlignIO.read(directory / "hmmer_source.sto", "stockholm")
        if (len(annotated) != len(members) or {record.id for record in annotated} != {gene_to_token[gene] for gene in members}
                or any(str(record.seq).replace("-", "").replace(".", "").upper() != aligned_residues[token_to_gene[record.id]]
                       for record in annotated)):
            raise ValueError("HMMER annotated source alignment changed sequence ownership")
        header = {}
        model_text = (directory / "profile.hmm").read_text()
        for line in model_text.splitlines():
            if line.startswith("HMM "):
                break
            if not line.strip():
                continue
            key, *values = line.split()
            if key in ("NAME", "LENG", "NSEQ", "EFFN"):
                if key in header or len(values) != 1:
                    raise ValueError("Invalid HMMER model header")
                header[key] = values[0]
        if (not {"NAME", "LENG", "NSEQ", "EFFN"} <= header.keys()
                or not model_text.startswith("HMMER3/f") or not model_text.rstrip().endswith("//")
                or header.get("NAME") != MODEL_NAME or int(header["LENG"]) != profile.length
                or int(header["NSEQ"]) != len(members) or float(header["EFFN"]) != len(members)):
            raise ValueError("HMMER model identity, length or effective count differs")
        result["model_header"] = header
        result["search_execution"] = _run_tool(search, directory, "search", environment, 600)
        measurements = read_hmmer_tables(directory / "sequence.tbl", directory / "domain.tbl", token_to_gene,
            {token: len(targets[gene]) for token, gene in token_to_gene.items()}, int(profile.length))
        if any(_record(pin["path"]) != pin for pin in result["tools"].values()):
            raise ValueError("HMMER binary changed during scoring")
        missing_members = [gene for gene in members if not measurements[gene]["reported"]]
        result.update({"measurements": measurements, "missing_members": missing_members,
            "missing_queries": sorted(gene for gene in queries if not measurements[gene]["reported"]),
            "calibration_complete": not missing_members,
            "minimum_member_score": None if missing_members else min(measurements[gene]["full_sequence_score"] for gene in members),
            "status": "hmmer_profile_scored"})
    except Exception as error:
        result.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        result["files"] = [_record(path) for path in sorted(directory.iterdir()) if path.is_file() and path != receipt]
        _atomic_write(receipt, json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
    return result
