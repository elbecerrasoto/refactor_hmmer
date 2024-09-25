#!/usr/bin/env python3
from pathlib import Path
from multiprocessing import Pool
from contextlib import ExitStack

from pyhmmer import hmmsearch
from pyhmmer.easel import SequenceFile
from pyhmmer.plan7 import HMMFile

QUERY = Path("queries/YwqJ_PF14431.hmm")
GENOMES_FILE = Path("genomes.txt")

hmm_file = HMMFile(QUERY)

with open(GENOMES_FILE, "r") as genomes_file:
    genomes_paths = [Path(line.rstrip()) for line in genomes_file]

with ExitStack() as stack:

    genomes_files = [
        stack.enter_context(SequenceFile(genome_path, digital=True))
        for genome_path in genomes_paths
    ]

    def worker(genome_file):
        genome = genome_file.read_block()
        return hmmsearch(hmm_file, genome)

    with Pool() as pool:
        results = pool.map(worker, genomes_files)

hmm_file.close()
