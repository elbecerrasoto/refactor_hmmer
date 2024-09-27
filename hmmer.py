#!/usr/bin/env python3
import os
import re
import sys
import numpy as np
from pathlib import Path
from typing import Iterable, Union
from functools import partial
from multiprocessing import Pool

from pyhmmer import hmmsearch
from pyhmmer.easel import SequenceFile
from pyhmmer.plan7 import HMM, HMMFile

QUERIES_DIR = Path(sys.argv[1])
GENOMES_FILE = sys.argv[2]
OUT_FILE = Path(sys.argv[3])

GENOME_REGEX = re.compile(r"(GC[FA]_\d+\.\d)\.faa$")


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
        (evalue := hit.evalue),
        (start := hit.best_domain.env_from),
        (end := hit.best_domain.env_to),
        (pid_description := hit.description.decode("utf-8")),
        (query_description := hit.hits.query_name.decode("utf-8")),
    ]

    out = [str(i) for i in out]
    return "\t".join(out) + "\n"


def run_genomes(genome_paths, hmms_files):

    def run_genome(genome_path):
        genome_id = parse_genome(genome_path)

        with SequenceFile(genome_path, digital=True) as genome_file:
            genome = genome_file.read_block()
            results = hmmsearch(hmms_files, genome)

        tsv = []
        for top_hits in results:
            for hit in top_hits:
                parsed = parse_hit(hit)
                hit_tsv = f"{genome_id}\t{parsed}"
                tsv.append(hit_tsv)

        del results

        out = {}
        out[genome_id] = tsv

        return out

    merged = {}
    for resultD in map(run_genome, genome_paths):
        merged |= resultD

    return merged


if __name__ == "__main__":

    with open(GENOMES_FILE, "r") as genomes_file:
        genomes_paths = [Path(line.rstrip()) for line in genomes_file]
        n_batches = len(genomes_paths)
        batches = np.array_split(np.array(genomes_paths), n_batches)

    hmms_files = get_hmms(QUERIES_DIR)
    worker = partial(run_genomes, hmms_files=hmms_files)

    with Pool() as pool:
        results = pool.imap_unordered(worker, batches)
        pool.close()
        pool.join()

    merged = {}
    for result in results:
        merged |= result

    with open(OUT_FILE, "w") as tsv:
        for genome_id in merged:
            top_hits = merged[genome_id]
            for hit in top_hits:
                tsv.write(hit)
