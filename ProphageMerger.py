#!/usr/bin/env python3

import sys
import os
import subprocess
import click
from datetime import datetime

version = "1.0.0"
@click.version_option(version, "--version", "-v")

@click.command(
    context_settings=dict(help_option_names=['-h', '--help'], max_content_width=150),
    help='Usage:\n python ProphageMerger.py -g <genome chromosome sequence> '
    '-d <database directory> -o <output directory>'
)
@click.option("-g",
    '--genome_sequence',
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
    help='Path to genome chromosome sequence'
)
@click.option("-d",
    '--database',
    default="databases",
    type=click.Path(dir_okay=True, resolve_path=True),
    show_default=True,
    help='Directory to save the downloaded databases'
)
@click.option("-o",
    '--output_dir',
    default="OUTPUT",
    type=click.Path(dir_okay=True, resolve_path=True),
    show_default=True,
    help=('Output directory')
)
@click.option(
    '--dryrun',
    is_flag=True,
    default=False,
    show_default=True,
    help='Check rules to run and files to produce'
)
@click.option(
    '--conda_envs',
    default='',
    show_default=True,
    help='Directory to store conda environments.'
    ' By default, the ".snakemake" directory relative to the invocation directory is used'
)
@click.option(
    '--profile',
    default='slurm',
    show_default=True,
    help='Snakemake profile for cluster execution'
)

def run_prophagemerger(genome_sequence, database, output_dir, dryrun, conda_envs, profile):

    # get snakefile path
    script_dir=os.path.dirname(os.path.abspath(__file__))
    snakefile=os.path.join(script_dir, "workflow", "Snakefile")

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
        '{conda_envs} '
        '--profile {profile} --rerun-incomplete ' 
        '--printshellcmds --nolock --show-failed-logs '
        '{dryrun} '
        '--config genome_sequence={seq} '
        'database={db} '
        'results_dir={results}'
        ).format(
            snakefile=snakefile,
            conda_envs='' if conda_envs=='' else '--conda-prefix {}'.format(conda_envs),
            profile=profile,
            dryrun='--dryrun' if dryrun else '',
            seq=genome_sequence,
            db=database,
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
