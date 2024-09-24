#!/usr/bin/env python3
import os
import re
import sys
import numpy as np
from pathlib import Path
from typing import Iterable, Union
from functools import partial
from multiprocessing import Pool, cpu_count
from contextlib import ExitStack

from pyhmmer import hmmsearch
from pyhmmer.easel import SequenceFile
from pyhmmer.plan7 import HMM, HMMFile

QUERIES_DIR = Path(sys.argv[1])
OUT_FILE = Path(sys.argv[2])
BATCH_SIZE = int(sys.argv[3])
WORKERS = int(float(sys.argv[4]) * cpu_count())
TOP_HITS = sys.argv[5] == "True"  # Include Top Hits Object
GENOMES_FILE = sys.argv[6]

CHUNKSIZE = WORKERS

# QUERIES_DIR = Path("queries")
# OUT_FILE = Path("test0.tsv")
# BATCH_SIZE = 512
# WORKERS = int(2 * cpu_count())
# TOP_HITS = False  # Include Top Hits Object
# GENOMES_FILE = Path("genomes.txt")

GENOME_REGEX = re.compile(r"(GCF_\d+\.\d)\.faa$")


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


def run_genomes(genome_paths, hmms_files):

    genomes_ids = list(map(parse_genome, genome_paths))

    with ExitStack() as stack:
        genomes_files = [
            stack.enter_context(SequenceFile(genome_path, digital=True))
            for genome_path in genome_paths
        ]

        results = []
        genomes = [g.read_block() for g in genomes_files]

        for genome in genomes:
            results.append(hmmsearch(hmms_files, genome))

    merged = {}
    for idx, result in enumerate(results):

        genome_id = genomes_ids[idx]

        out = {}
        out[genome_id] = {"tsv": None, "top_hits": None}

        tsv = []
        for top_hits in result:
            for hit in top_hits:
                parsed = parse_hit(hit)
                hit_tsv = f"{genome_id}\t{parsed}"
                tsv.append(hit_tsv)

        if TOP_HITS:
            out[genome_id]["top_hits"] = results
        else:
            del result

        out[genome_id]["tsv"] = tsv

        merged |= out

    return merged


if __name__ == "__main__":

    with open(GENOMES_FILE, "r") as genomes_file:
        genomes_paths = [Path(line.rstrip()) for line in genomes_file]
        batches = np.array_split(np.array(genomes_paths), BATCH_SIZE)

    hmms_files = get_hmms(QUERIES_DIR)
    worker = partial(run_genomes, hmms_files=hmms_files)

    with Pool(WORKERS) as pool:
        results = pool.imap_unordered(worker, batches, chunksize=CHUNKSIZE)
        pool.close()
        pool.join()

    # merge dictionaries
    merged = {}
    for resultD in results:
        merged |= resultD

    # write down tsv
    with open(OUT_FILE, "w") as tsv:
        for genome_id in merged:
            hits = merged[genome_id]["tsv"]

            for hit in hits:
                tsv.write(hit)
