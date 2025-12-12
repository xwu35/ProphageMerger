#!/usr/bin/env python3

import sys
import os
import subprocess
import click
from datetime import datetime

version = "1.1.0"
@click.version_option(version, "--version", "-v")

def validate_test_run(ctx, param, value):
    """
    A callback to make --genome_sequence required if --test_run is not specified.
    """
    if not ctx.params.get('test_run') and value is None:
        raise click.BadParameter(f"Option '{param.name}' is required unless --test_run is used.")
    return value

@click.command(
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=150),
    help='Usage:\n ProphageMerger.py -g <chromosome genome sequence> -o <output directory>'
)
@click.option("-g",
    '--genome_sequence',
    callback=validate_test_run,
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
    help='Path to chromosome genome sequence (required unless --test_run is used)'
)
@click.option("-o",
    '--output_dir',
    default="OUTPUT",
    type=click.Path(dir_okay=True, resolve_path=True),
    show_default=True,
    help=('Output directory')
)
@click.option(
    '--test_run',
    is_flag=True,
    default=False,
    show_default=True,
    help='ProphageMerger test run'
)
@click.option(
    '--dryrun',
    is_flag=True,
    default=False,
    show_default=True,
    help='Check rules to run and files to produce'
)
@click.option(
    '--profile',
    default='slurm',
    show_default=True,
    help='Snakemake profile for cluster execution'
)

def run_prophagemerger(genome_sequence, output_dir, test_run, dryrun, profile):

    # get snakefile and conda_envs path
    script_dir=os.path.dirname(os.path.abspath(__file__))
    snakefile=os.path.join(script_dir, "workflow", "Snakefile")
    conda_envs=os.path.join(script_dir, "conda_envs")

    if test_run:
        genome_sequence=os.path.join(script_dir, "test_data", "genome_sequences.fa")
        output_dir="test_output"

    # write run log if it is not a dry run
    if not dryrun:
        os.makedirs(output_dir, exist_ok=True)
        logfile = os.path.join(output_dir, f"{os.path.basename(output_dir)}_run.log")
        with open(logfile, "w") as log:
            log.write("================ProphageMerger run log==============\n")
            log.write(f"Start time: {datetime.now()}\n")
            log.write(f"ProphageMerger version: {version}\n")
            log.write(f"Genome sequence: {genome_sequence}")

    cmd = (
        'snakemake --snakefile {snakefile} '
        '--use-conda --conda-frontend mamba '
        '--conda-prefix {conda_envs} '
        '--profile {profile} --rerun-incomplete ' 
        '--printshellcmds --nolock --show-failed-logs '
        '{dryrun} '
        '--config genome_sequence={seq} '
        'results_dir={results}'
        ).format(
            snakefile=snakefile,
            conda_envs=conda_envs,
            profile=profile,
            dryrun='--dryrun' if dryrun else '',
            seq=genome_sequence,
            results=output_dir
            )

    # run snakemake with command-line config
    try:
        subprocess.run(cmd, check=True, shell=True)
    except subprocess.CalledProcessError:
        print("Snakemake failed. see log for details.", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    run_prophagemerger()
