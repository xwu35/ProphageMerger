#!/usr/bin/env python3

import os
import click
import pandas as pd
import numpy as np

@click.command(
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=150),
    help='Usage:\n python fill_in_nonprophage_coordinates.py -i <Prophage region coordinates> -t <Trimmed read counts> -o <Output file name>'
)
@click.option('-i',
    '--input',
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    help='Prophage region coordinates'
)
@click.option('-g',
    '--genome',
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    help='Genome length'
)
@click.option('-o',
    '--output',
    default="whole_genome_coordinates.txt",
    type=click.File("w"),
    show_default=True,
    help=('Output file name')
)

def process_data(input, genome, output):
    """
    fill in the gaps between prophage regions as nonprophage coordinates
    """
    
    # read in genome length
    genome = pd.read_csv(genome, sep='\t')
    genome_length_dict = dict(zip(genome['#name'], genome['length']))

    # read in prophage region coordinates
    samplenames = os.path.basename(input).split('_final_coordinates.txt')[0]
    prophage = pd.read_csv(input, sep='\t')[['chr', 'start', 'end', 'region']]
    prophage = prophage.rename(columns={'region': 'name'})
    prophage['gieStain'] = 'acen'
    
    # find gaps
    min_len = 0
    max_len = genome_length_dict[samplenames]
    chr_name = prophage['chr'].unique()
    # create a previous end column
    prophage['prev_end'] = prophage['end'].shift(1).fillna(min_len-1).astype(int)
    gaps_start = prophage['prev_end'] + 1
    gaps_end = prophage['start'] - 1
    # keep gaps where start < end
    repeated_chr = np.repeat(chr_name, len(gaps_start))
    gaps = pd.DataFrame({'chr': repeated_chr, 'start': gaps_start, 'end': gaps_end})
    gaps = gaps[gaps['start'] < gaps['end']]
    # handle the final gap after the last row
    if prophage['end'].iloc[-1] < max_len:
        final_gap = pd.DataFrame({'chr': chr_name, 'start': [prophage['end'].iloc[-1] + 1], 'end': [max_len]})
        gaps = pd.concat([gaps, final_gap])
    gaps['name'] = ['n' + str(i) for i in range(1, len(gaps) + 1)]
    gaps['gieStain'] = 'gneg'
    
    # merge with prophage region
    prophage = prophage.drop('prev_end', axis=1)
    whole_genome = pd.concat([gaps, prophage]).sort_values(by='start')

    # save the filtered df
    whole_genome.to_csv(output, sep='\t', index=False, header=True)

if __name__ == '__main__':
    process_data()
