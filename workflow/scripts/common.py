#!/usr/bin/env python

import os
import glob
import pandas as pd

# function to get genome names and sequences path
def parse_names_and_sequences(metadata, seq_dir=""):
    sequence_map = {}
    df = pd.read_table(metadata).set_index("genome", drop=False)
    genomes = df.index.tolist()
    for genome in genomes:
      genome_seq = df.loc[genome, "seq_name"]
      sequence_map[genome] = os.path.join(seq_dir, genome_seq)
    return genomes, sequence_map

