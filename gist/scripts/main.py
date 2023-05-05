import click
from gist.gist import GetSubSamplingByState, GetSimilarGenomes
from gist.util import read_get_states_input, check_dir, check_file
import os
import sys


@click.group()
def cli():
    "This tool englobes a set of functions to perform subsampling of gisaid information using augur."
    pass


@cli.command()
@click.option("--sequences", help="path to gisaid sequences tar.xz", required=True)
@click.option("--metadata", help="path to gisaid metadata tar.xz", required=True)
@click.option("--ncov_dir", help="path to ncov directory", required=True)
@click.option("--output_dir", help="path to output directory", default=".")
@click.option("--threads", help="threads", default=2, type=int)
@click.argument("subsampling_json_schema")
def get_states(
    sequences, metadata, ncov_dir, output_dir, threads, subsampling_json_schema
):
    "Get brazilian state sequences based on json input"
    GetSubSamplingByState(
        sequences, metadata, ncov_dir, output_dir, threads, subsampling_json_schema
    )


@cli.command()
@click.option("--input", help="path to gisaid subsampled sequences", required=True)
@click.option("--sequences", help="path to gisaid genomes fasta file", required=True)
@click.option("--metadata", help="path to gisaid metadata tsv file", required=True)
@click.option("--output_dir", help="path to output directory", default=".")
@click.option("--threads", help="threads", default=2, type=int)
@click.argument("get_similar_genomes_sampling_schema")
def get_genomes(
    input, sequences, metadata, output_dir, threads, get_similar_genomes_sampling_schema
):
    "Get gisaid similar genomes based on json input"
    GetSimilarGenomes(
        input,
        sequences,
        metadata,
        output_dir,
        threads,
        get_similar_genomes_sampling_schema,
    )
