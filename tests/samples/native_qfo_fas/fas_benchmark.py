#!/usr/bin/env python3
import concurrent.futures
import csv
import itertools
import json
import logging
import multiprocessing
import subprocess
import tempfile
import random

import numpy
import sqlite3
from pathlib import Path
from typing import TextIO
from tqdm import tqdm

from JSON_templates import write_assessment_dataset
from helpers import auto_open, load_json_file

logger = logging.getLogger("FAS-Benchmark")
MAX_PAIRS_COMPUTE = 9_000

def load_precomputed_fas_scores(f: Path):
    dat = load_json_file(str(f))
    data = {}
    for pair, vals in dat.items():
        p1, p2 = pair.split('_')
        if p1 > p2: p1, p2 = p2, p1
        try:
            data[(p1, p2)] = float(numpy.array(vals, dtype="float").mean())
        except ValueError as e:
            logger.debug("FAS score for (%s,%s) is not parseable: %s", p1, p2, e)
    return data

def generate_prot_to_annoationfile_map(annotations: Path):
    def map_one_file(json_file:Path):
        """ Read annotation json file and return a dictionary
        where IDs are the protein IDs and values are their corresponding filename
        """
        tax = json_file.stem
        with json_file.open(mode='rb') as jf:
            anno_dict = json.load(jf)
        return {id_: tax for id_ in anno_dict['feature']}

    taxa_dict = {}
    with concurrent.futures.ThreadPoolExecutor() as ex:
        for res in ex.map(map_one_file, annotations.glob("*.json")):
            taxa_dict.update(res)
    return taxa_dict


def compute_fas_scores_for_pairs(pairs, prot2tax, annotations, nr_cpus):
    logger.info("For %d orthologous pairs we don't have precomputed FAS scores", len(pairs))
    if len(pairs) > MAX_PAIRS_COMPUTE:
        logger.warning("Too many pairs to analyse, will sub-sample to %d", MAX_PAIRS_COMPUTE)
        random.shuffle(pairs)
        pairs = pairs[:MAX_PAIRS_COMPUTE]

    with tempfile.TemporaryDirectory() as tmp:
        missing_fn = Path(tmp) / "missing.txt"
        with missing_fn.open(mode='wt') as fh:
            for p1, p2 in pairs:
                fh.write(f"{p1}\t{prot2tax[p1]}\t{p2}\t{prot2tax[p2]}\n")

        fas_cmds = ["fas.runMultiTaxa", "--input", missing_fn, "-a", annotations,
                    "-o", Path(tmp)/"fas_res", "--bidirectional", "--tsv", "--domain", "--no_config", "--json",
                    "--mergeJson", "--outName", "computed_results", "--max_cardinality", "40", "--paths_limit", "15",
                    "--pairLimit", "30000", "--cpus", str(nr_cpus)]
        logger.info("running subprocess: %s", fas_cmds)
        res = subprocess.run(fas_cmds) # , capture_output=True)
        if res.returncode != 0:
            logger.error("Computing fas.runMultiTaxa failed: %s", res.stderr)
            scores = {}
        else:
            scores = load_precomputed_fas_scores(Path(tmp) / "fas_res" / "computed_results.json")
    return scores


def compute_fas_benchmark(precomputed_scores: Path, annotations: Path, db_path: Path, nr_cpus: int, raw_out: TextIO, limited_species=False):
    def iter_all_orthologs(species=None):
        cur = con.cursor()
        query = "SELECT DISTINCT p1.uniprot_id, p2.uniprot_id FROM orthologs JOIN proteomes as p1 ON orthologs.prot_nr1 = p1.prot_nr JOIN proteomes as p2 ON orthologs.prot_nr2 = p2.prot_nr WHERE p1.uniprot_id < p2.uniprot_id"
        if species is not None:
            sp_lst = ','.join(['"%s"' % v for v in species])
            query += " AND p1.species IN (%s) AND p2.species IN (%s)" % (sp_lst, sp_lst)
        cur.execute(query)
        cur.arraysize = 50000
        while True:
            chunk = cur.fetchmany()
            if len(chunk) == 0:
                break
            yield from chunk

    con = sqlite3.connect(db_path)
    scores_lookup = load_precomputed_fas_scores(precomputed_scores)
    prot_2_tax_map = generate_prot_to_annoationfile_map(annotations)
    logger.info("feature annotations available for %d proteins", len(prot_2_tax_map))
    logger.debug(list(itertools.islice(prot_2_tax_map.items(), 30)))

    scores = []; missing_pairs = []; no_fa = 0
    species = ("HUMAN", "MOUSE", "RATNO", "YEAST", "ECOLI", "ARATH") if limited_species else None
    for p1, p2 in tqdm(iter_all_orthologs(species)):
        if '_' in p1 or '_' in p2:
            # skipping uniprot_ids, only uniprot accessions!
            continue
        if (p1, p2) in scores_lookup:
            scores.append((p1, p2))
        elif p1 in prot_2_tax_map and p2 in prot_2_tax_map:
            missing_pairs.append((p1, p2))
        else:
            logger.info("No feature architecture found for relation %s/%s. Skipping", p1, p2)
            no_fa += 1

    nr_orthologs = len(scores) + len(missing_pairs)
    logger.info("%d pairs precomputed, %d missing (will compute); %d no feature annotations",
                len(scores), len(missing_pairs), no_fa)
    frac_precomputed = len(scores) / (len(scores) + len(missing_pairs))
    logger.info("ratio of precomputed vs missing pairs: ~%.0f:%.0f", 100*frac_precomputed, 100*(1-frac_precomputed))
    if len(missing_pairs) > 0:
        compute_nr = min(len(missing_pairs), MAX_PAIRS_COMPUTE)
        nr_precomp_maintain_frac = round( compute_nr * frac_precomputed / (1-frac_precomputed))
        logger.info("to maintain ratio of precomputed vs missing, we will compute %d new pairs "
                    "and sample %d precomputed pairs", compute_nr, nr_precomp_maintain_frac)
        random.shuffle(scores)
        scores = scores[:nr_precomp_maintain_frac]

        score2 = compute_fas_scores_for_pairs(missing_pairs, prot_2_tax_map, annotations, nr_cpus)
        scores_lookup.update(score2)
    csv_writer = csv.writer(raw_out, dialect="excel-tab")
    csv_writer.writerow(("Acc1", "Acc2", "FAS"))
    scores_list = []
    for part, pairs in zip(("precomputed", "missing"), (scores, missing_pairs)):
        score_part = []
        for pair in pairs:
            try:
                score = scores_lookup[pair]
                csv_writer.writerow((pair[0], pair[1], score))
                score_part.append(score)
            except KeyError:
                pass
        scores_list.extend(score_part)
        score_part = numpy.array(score_part, dtype="float")
        logger.info("FAS score[%s]: %f +- %f [N=%d]", part, score_part.mean(),
                    score_part.std(ddof=1) / numpy.sqrt(numpy.size(score_part)),
                    len(score_part))

    fas_scores = numpy.array(scores_list, dtype="float")
    fas_mean = fas_scores.mean()
    fas_sem = fas_scores.std(ddof=1) / numpy.sqrt(numpy.size(fas_scores))

    logger.info(f"FAS_mean: {fas_mean} +- {fas_sem}; nr_orthologs: {nr_orthologs}; sample_size: {len(fas_scores)} vs {numpy.size(fas_scores)}")

    metrics = [{"name": "FAS", "value": float(fas_mean), "stderr": float(fas_sem)},
               {"name": "NR_ORTHOLOGS", "value": nr_orthologs, "stderr": 0},
               ]
    return metrics




def write_assessment_json_stub(fn, community, participant, result, challenge):
    stubs = []
    for metric in result:
        id_ = "{}:{}_{}_{}_A".format(community, challenge, metric['name'],
                                     participant.replace(' ', '-').replace('_', '-'))
        stubs.append(write_assessment_dataset(id_, community, challenge, participant, metric['name'], metric['value'],
                                              metric.get('stderr', 0)))
    with auto_open(fn, 'wt') as fout:
        json.dump(stubs, fout, sort_keys=True, indent=4, separators=(',', ': '))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="FAS benchmark assessment")
    parser.add_argument('--assessment-out', required=True, help="Path where assessment json file will be stored")
    parser.add_argument('--outdir', required=True, help="Folder to store the raw output file in")
    parser.add_argument('--db', required=True, help="Path to sqlite database with pairwise predictions")
    parser.add_argument('--com', required=True, help="community id")
    parser.add_argument('--fas-precomputed-scores', required=True, help="Path to text file with VGNC asserted orthologs")
    parser.add_argument('--fas-data', help="Path to the input data of protein to feature architecture mapping. If not "
                                          "provided, missing pairs in the fas-precomputed-sores will be simply skipped "
                                          "and not computed on the fly")
    parser.add_argument('--limited-species', action="store_true", help="run on limited species set (6 species)")
    parser.add_argument('--participant', required=True, help="Name of participant method")
    parser.add_argument('--log', help="Path to log file. Defaults to stderr")
    parser.add_argument('--cpus', type=int, help="nr of cpus to use. defaults to all available cpus")
    parser.add_argument('-d', '--debug', action="store_true", help="Set logging to debug level")
    conf = parser.parse_args()

    log_conf = {'level': logging.INFO, 'format': "%(asctime)-15s %(levelname)-7s: %(message)s"}
    if conf.log is not None:
        log_conf['filename'] = conf.log
    if conf.debug:
        log_conf['level'] = logging.DEBUG
    if conf.cpus is None:
        conf.cpus = multiprocessing.cpu_count()
    logging.basicConfig(**log_conf)
    logger.info("running fas_benchmark with following arguments: {}".format(conf))
    challenge = "FAS"
    
    outdir = Path(conf.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    outfn_path = outdir / "{}_{}_raw.txt.gz".format(challenge, conf.participant.replace(' ', '-').replace('_', '-'))

    with auto_open(str(outfn_path), 'wt') as raw_out_fh:
        res = compute_fas_benchmark(Path(conf.fas_precomputed_scores), Path(conf.fas_data), Path(conf.db), conf.cpus, raw_out_fh, limited_species=conf.limited_species)
    write_assessment_json_stub(conf.assessment_out, conf.com, conf.participant, res, challenge)
