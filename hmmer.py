#!/usr/bin/env python3
import os
import re
import sys
from pathlib import Path
from typing import Iterable, Union
from functools import partial
from multiprocessing import Pool

from pyhmmer import hmmsearch
from pyhmmer.easel import SequenceFile
from pyhmmer.plan7 import HMM, HMMFile

QUERIES_DIR = Path(sys.argv[1])
OUT_FILE = Path(sys.argv[2])
GENOMES = sys.argv[3:]

GENOME_REGEX = re.compile(r"(GCF_\d+\.\d)\.faa$")

WORKERS = 12
TOP_HITS = False  # Include Top Hits Object


class HMMFiles(Iterable[HMM]):
    def __init__(self, *files: Union[str, bytes, os.PathLike]):
        self.files = files

    def __iter__(self):
        for file in self.files:
            with HMMFile(file) as hmm_file:
                yield from hmm_file


def get_hmms(queries_path):
    queries_path = Path(queries_path)
    hmms_files = HMMFiles(*queries_path.iterdir())
    return hmms_files


def parse_genome(genome_path):
    genome_path = str(genome_path)
    genome = re.search(GENOME_REGEX, genome_path).group(1)
    return genome


def parse_hit(hit):
    out = [
        (pid := hit.name.decode("utf-8")),
        (query := hit.hits.query_accession.decode("utf-8")),
        (score := hit.score),
        (start := hit.best_domain.env_from),
        (end := hit.best_domain.env_to),
        (included := hit.included),
        (reported := hit.reported),
        (pid_description := hit.description.decode("utf-8")),
        (query_description := hit.hits.query_name.decode("utf-8")),
    ]

    out = [str(i) for i in out]
    return "\t".join(out) + "\n"


def run_genome(genome_path, hmms_files):

    genome_id = parse_genome(genome)

    out = {}
    out[genome_id] = {"tsv": None, "top_hits": None}

    with SequenceFile(genome_path, digital=True) as genome_file:
        genome = genome_file.read_block()
        results = hmmsearch(hmms_files, genome)

    tsv = list()
    for top_hits in results:
        for hit in top_hits:
            parsed = parse_hit(hit)
            hit_tsv = f"{genome_id}\t{parsed}"
            tsv.append(hit_tsv)

    if TOP_HITS:
        out[genome_id]["top_hits"] = results
    else:
        del results

    out[genome_id]["tsv"] = tsv

    return out


if __name__ == "__main__":


    hmms_files = get_hmms(QUERIES_DIR)
    worker = partial(run_genome, hmms_files=hmms_files)

    with Pool(WORKERS) as pool:
        results = pool.map(worker, GENOMES)

    # merge dictionaries
    merged = {}
    for resultD in results:
        merged = merged | resultD

    # write down tsv
    with open(OUT_FILE, "w") as tsv:
        for genome_id in merged:
            hits = merged[genome_id]["tsv"]

            for hit in hits:
                tsv.write(hit)

        for genome in OUT:
            for hit_data in OUT[genome]:
                tsv.write(hit_data)
