#!/usr/bin/env python3

import click
import pandas as pd
from Bio import SeqIO
from io import StringIO

@click.command(
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=150),
    help='Usage:\n python combine_prophage_coordinates_arguments.py --virsorter2 <VirSorter2 boundary file> '
    '--checkv <CheckV provirus fasta> --genomad <GeNomad provirus file> --vibrant <Vibrant prophage file> '
    '--all_coordinates <Output all coordinates> '
    '--final_coordinates <Output maximized coordiantes> '
    '--bed_file <Output maximized 0-based coordiantes>'
)
@click.option(
    '--virsorter2',
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    help='VirSorter2 boundary file'
)
@click.option(
    '--checkv',
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    help='CheckV provirus fasta'
)
@click.option(
    '--genomad',
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    help='GeNomad provirus file'
)
@click.option(
    '--vibrant',
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    help='Vibrant prophage file'
)
@click.option(
    '--all_coordinates',
    default="all_coordinates.txt",
    type=click.File("w"),
    show_default=True,
    help=('Output all coordinates')
)
@click.option(
    '--final_coordinates',
    default="final_coordinates.txt",
    type=click.File("w"),
    show_default=True,
    help=('Output maximized coordinates')
)
@click.option(
    '--bed_file',
    default="final_coordinates_0-based.bed",
    type=click.File("w"),
    show_default=True,
    help=('Output maximized 0-based coordinates')
)

def process_data(virsorter2, checkv, genomad, vibrant, 
                all_coordinates, final_coordinates, bed_file):

    # VirSorter2 + CheckV
    # read in virsorter2 boundary table
    df_vir = pd.read_table(virsorter2)

    # parse headers from checkv proviruses.fna
    with open(checkv, "r") as fasta_handle:
        records = list(SeqIO.parse(fasta_handle, "fasta"))

    data = []

    # iterate through the sequence records 
    if records:
        for seq_record in records:
            header = seq_record.description
            # split the header by space
            header_parts = header.split()
            # split the second parts by / to get the coordinates and length
            header_pos = header_parts[1].split('/')
            # assign column names to the split parts
            chr_name=header_parts[0].rsplit('||', 1)[0].strip()
            seq_name = header_parts[0].rsplit('_', 1)[0].strip()
            seq_id = header_parts[0].rsplit('_', 1)[1].strip()
            input_seq_len = header_pos[1].strip()
            start = header_pos[0].split('-')[0].strip()
            end = header_pos[0].split('-')[1].strip()
            # append the extracted data as a row to the list
            data.append({'chr_name': chr_name, 
                        'seq_name': seq_name, 
                        'seq_id': seq_id, 
                        'input_seq_len': input_seq_len, 
                        'start': start, 
                        'end': end})
            # create dataframe from the list
            df_checkv = pd.DataFrame(data)
    else:
        df_checkv = pd.DataFrame(columns=['seq_name'])

    # merge the two dfs
    merged_vir_checkv = pd.merge(df_checkv, df_vir, left_on = 'seq_name', right_on = 'seqname_new', how = 'left')

    # calculate the coordinates and save it as a new df
    if not merged_vir_checkv.empty:
        df_vir_checkv = pd.DataFrame({'chr_name': merged_vir_checkv['chr_name'].astype(str),
                                'seq_name': merged_vir_checkv['seq_name'],
                                'start': merged_vir_checkv['trim_bp_start'].astype(int) + merged_vir_checkv['start'].astype(int) - 1,
                                'end': merged_vir_checkv['trim_bp_start'].astype(int) + merged_vir_checkv['end'].astype(int) - 1,
                                'method': 'VirSorter2 + CheckV'})
    else:
        column_names = ['seq_name','start','end','method']
        df_vir_checkv = pd.DataFrame(columns=column_names)

    # GeNomad
    df_genomad = pd.read_table(genomad)[['source_seq', 'seq_name','start','end']]
    df_genomad.rename(columns={'source_seq': 'chr_name'}, inplace=True)
    df_genomad[['chr_name']] = df_genomad[['chr_name']].astype(str)
    df_genomad['method'] = 'GeNomad'

    # Vibrant
    df_vibrant = pd.read_table(vibrant)[['scaffold', 'fragment', 'nucleotide start','nucleotide stop']]
    # unlike the others, vibrant doesn't remove the space and texts after space in the seq headers
    df_vibrant['scaffold'] = df_vibrant['scaffold'].str.split().str[0]
    df_vibrant['method'] = 'Vibrant'
    df_vibrant.columns = ['chr_name','seq_name','start','end', 'method']
    df_vibrant[['chr_name']] = df_vibrant[['chr_name']].astype(str)

    # Combine all coordinates
    combined_df = pd.concat([df_vir_checkv, df_genomad, df_vibrant]).sort_values(by=['chr_name','start'])

    if combined_df.empty:
        empty_df = pd.read_csv(StringIO("""NO PROPHAGE REGIONS FOUND"""))
        empty_df.to_csv(all_coordinates, index=False)
        empty_df.to_csv(final_coordinates, index=False)
        empty_df.to_csv(bed_file, index=False)
        
    else:
        # save the combined results
        combined_df.to_csv(all_coordinates, sep='\t', index=False, header=True)

        # merge overlapping regions
        # get the max of the previous rows within each chromosome; the first row will have NaN value, so fill it with 0
        combined_df['previous_end_max'] = combined_df.groupby('chr_name')['end'].shift(1).groupby(combined_df['chr_name']).cummax().fillna(0)
        # calculate the difference between row's start corrdinates and the previous rows max end coordinates to see if they overlap
        combined_df['diff'] = combined_df['start'] - combined_df['previous_end_max']
        # add a condition column, if the diff is > 0, then output True (a new region), otherwise False
        combined_df['condition'] = combined_df['diff'] > 0
        # use cumsum to create region identifiers
        combined_df['region_id'] = combined_df.groupby('chr_name')['condition'].cumsum()
        # find the minimum value of start coordinates and maximum value of end coordiantes for each region
        results = combined_df.groupby(['chr_name','region_id']).agg(
            min_start=('start', 'min'),
            max_end=('end', 'max')
        ).reset_index()
        # rename the columns
        results.columns = ['chr','region','start','end']
        # add string 'region' to column region
        results['region'] = 'region_' + results['region'].astype(str)
        # calculate the prophage length
        results['length'] = results['end'] - results['start'] + 1
        # save the results
        results.to_csv(final_coordinates, sep='\t', index=False, header=True)

        # output another 0-based coordinates for samtools
        results_0based = pd.DataFrame({'chr': results['chr'],
                                    'start': results['start'] - 1,
                                    'end': results['end'],
                                    'region': results['region']})

        results_0based.to_csv(bed_file, sep='\t', index=False, header=False)

if __name__ == '__main__':
    process_data()
