## GIsaid Subsampling Toolkit

Just a set of scripts organized as a toolkit for conducting sub-sample analyses based on the Augur tool, focused on Brazilian states.

### Install

> Requirements: python 3.9+

```bash
git clone https://github.com/dezordiPhD/gist.git
git add gist
conda env create -f env/gis_ubuntu.yml
conda activate gist
pip install .

## check installation
gist --help

## clone ncov nextstrain repository
git clone https://github.com/nextstrain/ncov.git
```

> Mac users should install ncbi+blast 

### Run

#### get-states

This mode get subsampling data based on specific lineages on specific brazilian states and other countries. The input json file should be configure as the template present on `templates/get_by_states.json`

```bash
gist get-states --ncov_dir ncov --sequences <gisaid_genomes.tar.xz> --metadata <gisaid_metadata.tar.xz> --threads <number_of_threads> templates/get_by_states.json              
```

#### get-genomes

Get gisaid genomes based on blast analysis. The input json file should be configured as the template present on `templates/get_similar_genomes.json`

```bash
gist get-genomes --input <query_genomes.fasta> --sequences <gisaid_genomes.fasta> --metadata <gisaid_metadata.tsv> templates/get_similar_genomes.json
```
