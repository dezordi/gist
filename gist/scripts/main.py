import click
from gist.gist import GetSubSamplingByState
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
@click.option("--ncov_dir", help="path to ncov dir", required=True)
@click.argument("subsampling_json_schema")
def get_states(sequences, metadata, ncov_dir, subsampling_json_schema):
    "Get brazilian state sequences based on json input"
    GetSubSamplingByState(sequences, metadata, ncov_dir, subsampling_json_schema)
