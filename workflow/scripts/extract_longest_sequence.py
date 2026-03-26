#!/usr/bin/env python3

import click
from Bio import SeqIO

@click.command(
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=150),
    help='Usage:\n python extract_longest_sequence.py -i <input fasta file> -o <output fasta file>'
)
@click.option('-i',
    '--input',
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    help='Input fasta'
)
@click.option('-n',
    '--name',
    required=True,
    help='Name for sequence header'
)
@click.option('-o',
    '--output',
    default="longest_sequence.fasta",
    show_default=True,
    help='Output fasta'
)

def extract_longest_sequence(input, name, output):

    longest_seq = None
    
    for record in SeqIO.parse(input, "fasta"):
        if longest_seq is None or len(record.seq) > len(longest_seq.seq):
            longest_seq = record

    # rename the sequence header
    longest_seq.id = name
    longest_seq.description = "" 

    # write out the longest sequence
    with open(output, "w") as output_handle:
        SeqIO.write(longest_seq, output_handle, "fasta")

if __name__ == '__main__':
    extract_longest_sequence()
