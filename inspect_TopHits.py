#!/usr/bin/env python3
import re
import os
from typing import Iterable, Union
from pathlib import Path
from pyhmmer import hmmsearch
from pyhmmer.easel import SequenceFile
from pyhmmer.plan7 import HMM, HMMFile

QUERY = Path("queries13")
GENOMES_FILE = Path("genomes.txt")
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

with open(GENOMES_FILE, "r") as genomes_file:
    genomes_paths = [Path(line.rstrip()) for line in genomes_file]


HMM_FILES = get_hmms(QUERY)
RESULTS = {}

for genome_path in genomes_paths:
    genome_id = parse_genome(genome_path)
    with SequenceFile(genome_path, digital=True) as genome_file:
        genome = genome_file.read_block()
        result = hmmsearch(HMM_FILES, genome)
        RESULTS[genome_id] = result

def a_tophits():
    for genome in RESULTS:
        print(genome)
        result = RESULTS[genome]
        for tophits in result:
            if len(tophits) > 0:
                return tophits


OUT_tophits = a_tophits()
