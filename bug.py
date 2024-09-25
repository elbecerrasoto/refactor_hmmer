#!/usr/bin/env python3
from pathlib import Path
from multiprocessing import Pool
from contextlib import ExitStack

from pyhmmer import hmmsearch
from pyhmmer.easel import SequenceFile
from pyhmmer.plan7 import HMM, HMMFile

QUERY = Path("queries/YwqJ_PF14431.hmm")
GENOMES_FILE = Path("genomes.txt")

with open(GENOMES_FILE, "r") as genomes_file:
    genomes_paths = [Path(line.rstrip()) for line in genomes_file]

with ExitStack() as stack:

    genomes_files = [
        stack.enter_context(SequenceFile(genome_path, digital=True))
        for genome_path in genomes_paths
    ]

    genomes = [g.read_block() for g in genomes_files]

    with HMMFile(QUERY) as hmm_file:

        def worker(genome):
            return hmmsearch(hmm_file, genome)

        with Pool() as pool:
            results = pool.imap_unordered(worker, genomes)
            pool.close()
            pool.join()

print(results)
# for i in results:
#     print(i)
