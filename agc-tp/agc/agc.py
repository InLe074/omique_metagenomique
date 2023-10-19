#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
import textwrap
from pathlib import Path
from collections import Counter
from typing import Iterator, Dict, List
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
# import nwalign3 as nw


def read_fasta(amplicon_file: Path, minseqlen: int) -> Iterator[str]:
    """Read a compressed FASTA file and extract sequences with length >= minseqlen.

    :param amplicon_file: (Path) Path to the FASTA.gz file.
    :param minseqlen: (int) Minimum amplicon sequence length.
    :return: A generator object that yields valid FASTA sequences (str).
    """
    with gzip.open(amplicon_file, 'rt') as f:
        sequence = ""
        for line in f:
            if line.startswith('>'):
                # A new sequence is starting
                if len(sequence) >= minseqlen:
                    yield sequence
                sequence = ""
            else:
                # Append the current line to the sequence
                sequence += line.strip()
        # Check the last sequence in the file
        if len(sequence) >= minseqlen:
            yield sequence



def dereplication_fulllength(amplicon_file: Path, minseqlen: int, mincount: int) -> Iterator[List]:
    """Dereplicate the set of sequence

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length
    :param mincount: (int) Minimum amplicon count
    :return: A generator object that provides a (list)[sequences, count] of sequence with a count >= mincount and a length >= minseqlen.
    """
    pass

def dereplication_fulllength(amplicon_file: Path, minseqlen: int, mincount: int) -> Iterator[List]:
    """Dereplicate the set of sequences.

    :param amplicon_file: (Path) Path to the FASTA.gz file.
    :param minseqlen: (int) Minimum amplicon sequence length.
    :param mincount: (int) Minimum amplicon count.
    :return: A generator object that yields [sequence, count] of unique sequences with count >= mincount.
    """
    sequence_counter = Counter()

    # Use the read_fasta generator to read sequences from the FASTA.gz file
    for sequence in read_fasta(amplicon_file, minseqlen):
        sequence_counter[sequence] += 1

    # Yield unique sequences with counts >= mincount in descending order of count
    for sequence, count in sequence_counter.items():
        if count >= mincount:
            yield [sequence, count]

def abundance_greedy_clustering(amplicon_file: Path, minseqlen: int, mincount: int, chunk_size: int, kmer_size: int) -> List:
    """Compute an abundance greedy clustering regarding sequence count and identity.
    Identify OTU sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length.
    :param mincount: (int) Minimum amplicon count.
    :param chunk_size: (int) A fournir mais non utilise cette annee
    :param kmer_size: (int) A fournir mais non utilise cette annee
    :return: (list) A list of all the [OTU (str), count (int)] .
    """
    pass


def write_OTU(OTU_list: List, output_file: Path) -> None:
    """Write the OTU sequence in fasta format.

    :param OTU_list: (list) A list of OTU sequences
    :param output_file: (Path) Path to the output file
    """
    pass


#==============================================================
# Main program
#==============================================================


if __name__ == '__main__':
    # Example usage:
    amplicon_file = Path("C:/Users/Bismilah/Desktop/M1-IPFB/M2_BI/Omique_1_projet_par_jour/Omique_2nd_jour/agc-tp/data/amplicon.fasta.gz")
    minseqlen = 50  # Replace with your desired minimum sequence length

    # for sequence in read_fasta(amplicon_file, minseqlen):
        # print(sequence)

    minseqlen = 50  # Replace with your desired minimum sequence length
    mincount = 5  # Replace with your desired minimum count

    for sequence, count in dereplication_fulllength(amplicon_file, minseqlen, mincount):
        print(f"Count: {count}, Sequence: {sequence}")