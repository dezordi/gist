from gist.util import (
    read_get_states_input,
    read_get_similar_genomes_input,
    check_dir,
    check_file,
    check_mask_pos_file,
)
import sys
import click
import os
import subprocess
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
from Bio.Blast.Applications import NcbimakeblastdbCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO


class GetSubSamplingByState:
    def __init__(
        self,
        sequences: str,
        metadata: str,
        ncov_dir: str,
        output_dir: str,
        threads: int,
        subsampling_json_schema: str,
    ) -> dict:
        try:
            for file in sequences, metadata, subsampling_json_schema:
                check_file(file)
            for directory in ncov_dir, output_dir:
                check_dir(directory)

        except Exception as err:
            click.echo(f"Failed to validate input files: {err}", err=True)
            sys.exit(1)
        try:
            valid_subsampling_schema = read_get_states_input(subsampling_json_schema)
        except Exception as err:
            click.echo(f"Invalid json input: {err}", err=True)
            sys.exit(1)

        self.sequences = sequences
        self.metadata = metadata
        self.ncov_dir_scripts = os.path.join(ncov_dir, "scripts")
        self.threads = threads
        self.output_job_dir = os.path.join(
            output_dir, valid_subsampling_schema.job_name
        )
        self.valid_subsampling_schema = valid_subsampling_schema
        self.sanitize_sequences()
        self.sanitize_metadata()
        self.get_global()

    def _augur_index(self) -> None:
        augur_index_cmd = f"augur index \
                            --sequences {self.output_job_dir}/sequences_gisaid.fasta.gz \
                            --output {self.output_job_dir}/sequence_index_gisaid.tsv.gz"
        run_augur_index_cmd = subprocess.run(augur_index_cmd, shell=True)

    def sanitize_sequences(self) -> None:
        print("Sanitize and index sequences")
        if os.path.exists(self.output_job_dir) == False:
            print(f"Creating {self.output_job_dir}")
            os.mkdir(self.output_job_dir)
        sanitize_sequences_cmd = f"python {self.ncov_dir_scripts}/sanitize_sequences.py \
                                    --sequences {self.sequences} \
                                    --strip-prefixes hCoV-19/ SARS-CoV-2/ \
                                    --output {self.output_job_dir}/sequences_gisaid.fasta.gz"
        run_sanitize_sequences_cmd = subprocess.run(sanitize_sequences_cmd, shell=True)
        self._augur_index()

    def sanitize_metadata(self) -> None:
        print("Sanitize metadata")
        sanitize_metadata_cmd = f"python {self.ncov_dir_scripts}/sanitize_metadata.py \
                                    --metadata {self.metadata} \
                                    --metadata-id-columns strain name 'Virus name' \
                                    --database-id-columns 'Accession ID' gisaid_epi_isl genbank_accession \
                                    --parse-location-field Location \
                                    --rename-fields 'Virus name=strain' Type=type 'Accession ID=gisaid_epi_isl' 'Collection date=date' 'Sequence length=length' Host=host 'Pango lineage=pango_lineage' 'Host=host'\
                                    --strip-prefixes hCoV-19/ SARS-CoV-2/ \
                                    --output {self.output_job_dir}/metadata_gisaid.tsv.gz"
        run_sanitize_metadata_cmd = subprocess.run(sanitize_metadata_cmd, shell=True)

    def _filter_state(self, state: object) -> str:
        print(f"filter {state.name} genomes")
        output = os.path.join(self.output_job_dir, f"state_{state.sigla}_sub.txt")
        filter_by_state_cmd = f"""augur filter \
                                --metadata {self.output_job_dir}/metadata_gisaid.tsv.gz \
                                --sequence-index {self.output_job_dir}/sequence_index_gisaid.tsv.gz \
                                --min-date {self.valid_subsampling_schema.min_date} \
                                --max-date {self.valid_subsampling_schema.max_date} \
                                --min-length {str(self.valid_subsampling_schema.min_genome_len)} \
                                --query "(country == 'Brazil') & (division == '{state.name}') & (pango_lineage == {self.valid_subsampling_schema.target_lineages}) & (host == 'Human') " \
                                --subsample-max-sequences {str(state.max_genomes)} \
                                --exclude-ambiguous-dates-by any \
                                --group-by pango_lineage year month \
                                --output-strains {output}"""
        run_filter_by_state_cmd = subprocess.run(filter_by_state_cmd, shell=True)

        return output

    def _filter_country(self, country: object) -> str:
        print(f"filter {country.name} genomes")
        output = os.path.join(self.output_job_dir, f"country_{country.sigla}_sub.txt")
        filter_by_country_cmd = f"""augur filter \
                                --metadata {self.output_job_dir}/metadata_gisaid.tsv.gz \
                                --sequence-index {self.output_job_dir}/sequence_index_gisaid.tsv.gz \
                                --min-date {self.valid_subsampling_schema.min_date} \
                                --max-date {self.valid_subsampling_schema.max_date} \
                                --min-length {str(self.valid_subsampling_schema.min_genome_len)} \
                                --query "(country == '{country.name}') & (pango_lineage == {self.valid_subsampling_schema.target_lineages}) & (host == 'Human')" \
                                --subsample-max-sequences {str(country.max_genomes)} \
                                --exclude-ambiguous-dates-by any \
                                --group-by pango_lineage year month \
                                --output-strains {output}"""
        run_filter_by_country_cmd = subprocess.run(filter_by_country_cmd, shell=True)

        return output

    def _filter_outgroup(self) -> str:
        print(f"filter outgroup genomes")
        output = os.path.join(self.output_job_dir, "outgroup_sub.txt")
        filter_outgroup_cmd = f"""augur filter \
                                --metadata {self.output_job_dir}/metadata_gisaid.tsv.gz \
                                --sequence-index {self.output_job_dir}/sequence_index_gisaid.tsv.gz \
                                --max-date {self.valid_subsampling_schema.min_date} \
                                --min-length {str(self.valid_subsampling_schema.min_genome_len)} \
                                --query "pango_lineage == {self.valid_subsampling_schema.outgroup_lineages} & (host == 'Human')" \
                                --subsample-max-sequences 30 \
                                --exclude-ambiguous-dates-by any \
                                --group-by pango_lineage year month \
                                --output-strains {output}"""
        run_filter_outgroup_cmd = subprocess.run(filter_outgroup_cmd, shell=True)

        return output

    def _filter_gisaid(self) -> str:
        print(f"filter gisaid")
        subsampling_files = self.get_samples()
        output = os.path.join(self.output_job_dir, f"gisaid_sub.txt")
        filter_gisaid_cmd = f"""augur filter \
                                --metadata {self.output_job_dir}/metadata_gisaid.tsv.gz \
                                --sequence-index {self.output_job_dir}/sequence_index_gisaid.tsv.gz \
                                --min-date {self.valid_subsampling_schema.min_date} \
                                --max-date {self.valid_subsampling_schema.max_date} \
                                --min-length {str(self.valid_subsampling_schema.min_genome_len)} \
                                --exclude {subsampling_files} \
                                --query "(country != 'Brazil') & (pango_lineage == {self.valid_subsampling_schema.target_lineages}) & (host == 'Human') " \
                                --exclude-ambiguous-dates-by any \
                                --output-strains {output}"""

        run_filter_gisaid_cmd = subprocess.run(filter_gisaid_cmd, shell=True)

        return output

    def get_samples(self) -> str:
        sub_sampling_target_states = []
        sub_sampling_files = []

        with ThreadPoolExecutor(max_workers=self.threads) as executor:
            job_pool = []
            for state in self.valid_subsampling_schema.states:
                job = executor.submit(self._filter_state, state)
                job.state_sigla = state.sigla
                job_pool.append(job)

            for job in job_pool:
                subsampled_file = job.result()
                if job.state_sigla in self.valid_subsampling_schema.target_states:
                    sub_sampling_target_states.append(subsampled_file)
                else:
                    sub_sampling_files.append(subsampled_file)

        if self.valid_subsampling_schema.countries:
            with ThreadPoolExecutor(max_workers=self.threads) as executor:
                job_pool = []
                for country in self.valid_subsampling_schema.countries:
                    job = executor.submit(self._filter_country, country)
                    job_pool.append(job)
                for job in job_pool:
                    subsampled_file = job.result()
                    sub_sampling_files.append(subsampled_file)

        sub_sampling_files.append(self._filter_outgroup())

        if self.valid_subsampling_schema.target_states:
            sub_sampling_target_states = " ".join(sub_sampling_target_states)
        sub_sampling_files = " ".join(sub_sampling_files)

        if self.valid_subsampling_schema.target_states:
            get_samples_target_states_cmd = f"augur filter \
                                --metadata {self.output_job_dir}/metadata_gisaid.tsv.gz \
                                --sequence-index {self.output_job_dir}/sequence_index_gisaid.tsv.gz \
                                --sequences {self.output_job_dir}/sequences_gisaid.fasta.gz \
                                --exclude-all \
                                --include {sub_sampling_target_states} \
                                --output-metadata {self.output_job_dir}/subsampled_target_states_metadata_gisaid.tsv \
                                --output-sequences {self.output_job_dir}/subsampled_target_states_sequences_gisaid.fasta"
            run_get_samples_target_states_cmd = subprocess.run(
                get_samples_target_states_cmd, shell=True
            )

        get_samples_cmd = f"augur filter \
                                --metadata {self.output_job_dir}/metadata_gisaid.tsv.gz \
                                --sequence-index {self.output_job_dir}/sequence_index_gisaid.tsv.gz \
                                --sequences {self.output_job_dir}/sequences_gisaid.fasta.gz \
                                --exclude-all \
                                --include {sub_sampling_files} \
                                --output-metadata {self.output_job_dir}/subsampled_metadata_gisaid.tsv \
                                --output-sequences {self.output_job_dir}/subsampled_sequences_gisaid.fasta"
        run_get_samples_cmd = subprocess.run(get_samples_cmd, shell=True)

        if self.valid_subsampling_schema.target_states:
            return sub_sampling_target_states + " " + sub_sampling_files

        return sub_sampling_files

    def get_global(self) -> None:
        print(f"creating gisaid base without sampling sequences")
        gisaid_samples = self._filter_gisaid()

        filter_by_state_cmd = f"""augur filter \
                                --metadata {self.output_job_dir}/metadata_gisaid.tsv.gz \
                                --sequence-index {self.output_job_dir}/sequence_index_gisaid.tsv.gz \
                                --sequences {self.output_job_dir}/sequences_gisaid.fasta.gz \
                                --exclude-all \
                                --include {gisaid_samples} \
                                --output-metadata {self.output_job_dir}/gisaid_filtered_metadata_gisaid.tsv \
                                --output-sequences {self.output_job_dir}/gisaid_filtered_sequences_gisaid.fasta"""
        run_filter_by_state_cmd = subprocess.run(filter_by_state_cmd, shell=True)


class GetSimilarGenomes:
    def __init__(
        self,
        input_file: str,
        sequences: str,
        metadata: str,
        output_dir: str,
        threads: int,
        get_similar_genomes_sampling_schema: dict,
    ) -> None:
        try:
            for file in (
                input_file,
                sequences,
                metadata,
                get_similar_genomes_sampling_schema,
            ):
                check_file(file)
            check_dir(output_dir)

        except Exception as err:
            click.echo(f"Failed to validate input files: {err}", err=True)
            sys.exit(1)
        try:
            valid_subsampling_schema = read_get_similar_genomes_input(
                get_similar_genomes_sampling_schema
            )
        except Exception as err:
            click.echo(f"Invalid json input: {err}", err=True)
            sys.exit(1)
        self.input_file = input_file
        self.sequences = sequences
        self.metadata = metadata
        self.output_job_dir = os.path.join(
            output_dir, valid_subsampling_schema.job_name
        )
        self.threads = threads
        self.valid_subsampling_schema = valid_subsampling_schema
        self.get_gisaid_filtered_genomes()

    def create_blast_database(self) -> None:
        make_db = NcbimakeblastdbCommandline(dbtype="nucl", input_file=self.sequences)
        stdout, stderr = make_db()

    def perform_blast(self) -> None:
        if os.path.exists(self.output_job_dir) == False:
            print(f"Creating {self.output_job_dir}")
            os.mkdir(self.output_job_dir)
        output = os.path.join(self.output_job_dir, "gisaid_blastn.tsv")
        blast = NcbiblastnCommandline(
            query=self.input_file,
            db=self.sequences,
            out=output,
            outfmt="6 qseqid sseqid pident evalue bitscore qcovhsp",
            task="megablast",
            evalue=0.00001,
            num_threads=self.threads,
        )
        stdout, stderr = blast()

    def filter_blast_results(self) -> list:
        blast_results = os.path.join(self.output_job_dir, "gisaid_blastn.tsv")
        output = os.path.join(self.output_job_dir, "gisaid_blastn_filtered.tsv")
        blast_columns = ["qseqid", "sseqid", "pident", "evalue", "bitscore", "qcovhsp"]
        blast_results_df = pd.read_csv(
            blast_results, names=blast_columns, delimiter="\t"
        )
        blast_results_df = blast_results_df[blast_results_df["qcovhsp"] >= 99.9]
        blast_results_df = blast_results_df[
            (blast_results_df["pident"] >= self.valid_subsampling_schema.min_id)
            & (blast_results_df["pident"] <= self.valid_subsampling_schema.max_id)
        ]
        blast_results_df = blast_results_df.sort_values(by="bitscore", ascending=False)
        blast_results_df = blast_results_df.drop_duplicates(subset=["sseqid"])
        blast_results_df = blast_results_df.groupby("qseqid").head(valid_subsampling_schema.max_number_of_genomes_per_query)
        blast_results_df = blast_results_df.head(
            self.valid_subsampling_schema.max_number_of_similar_genomes
        )
        gisaid_filtered_matches = blast_results_df["sseqid"].tolist()
        blast_results_df.to_csv(output, sep="\t", index=False)

        return gisaid_filtered_matches

    def get_gisaid_filtered_genomes(self):
        self.create_blast_database()
        self.perform_blast()
        gisaid_filtered_matches = self.filter_blast_results()
        output_sequences = os.path.join(
            self.output_job_dir, "gisaid_similar_genomes.fasta"
        )
        output_metadata = os.path.join(
            self.output_job_dir, "gisaid_similar_genomes.tsv"
        )
        selected_records = []
        with open(self.sequences, "r") as handle:
            records = SeqIO.parse(handle, "fasta")
            for record in records:
                if record.id in gisaid_filtered_matches:
                    selected_records.append(record)
        with open(output_sequences, "w") as output:
            SeqIO.write(selected_records, output, "fasta")

        metadata_df = pd.read_csv(self.metadata, sep="\t")
        selected_df = metadata_df[metadata_df["strain"].isin(gisaid_filtered_matches)]
        selected_df.to_csv(output_metadata, sep="\t", index=False)


class GetAlignment:
    def __init__(
        self,
        input_file: str,
        reference: str,
        output_dir: str,
        threads: int,
        mask_pos: str,
    ) -> None:
        try:
            for file in (input_file, reference):
                check_file(file)
            check_dir(output_dir)
            if mask_pos:
                check_file(mask_pos)
                check_mask_pos_file(mask_pos)

        except Exception as err:
            click.echo(f"Failed to validate input files: {err}", err=True)
            sys.exit(1)

        self.input_file = input_file
        self.reference = reference
        self.output_dir = output_dir
        self.threads = threads
        self.mask_pos = mask_pos
        self.perform_alignment()

    def fix_sequence_names(self, alignment: str, renamed_alignment: str) -> None:
        with open(alignment, "r") as handle_in, open(
            renamed_alignment, "w"
        ) as output_handle:
            for record in SeqIO.parse(handle_in, "fasta"):
                if "|" in record.id:
                    record.id = record.id.split("|")[1]
                output_handle.write(
                    ">" + str(record.id) + "\n" + str(record.seq) + "\n"
                )

    def mask_alignment(self, alignment: str) -> None:
        alignment_seqs = SeqIO.parse(alignment, "fasta")
        positions_file = self.mask_pos
        positions = []
        with open(positions_file) as f:
            for line in f:
                start, end = line.strip().split("\t")
                positions.append((int(start), int(end)))

        output_path = os.path.join(self.output_dir, f"{os.path.basename(self.input_file)}.mask.algn")
        with open(output_path, "w") as output_handle:
            for record in alignment_seqs:
                sequence = record.seq
                for start, end in positions:
                    sequence = (
                        sequence[: start - 1] + "N" * (end - start + 1) + sequence[end:]
                    )
                record.seq = sequence
                output_handle.write(
                    ">" + str(record.id) + "\n" + str(record.seq) + "\n"
                )

    def perform_alignment(self):
        output = os.path.join(self.output_dir, f"sequences.algn.fa")
        renamed_alignment = os.path.join(self.output_dir, f"{os.path.basename(self.input_file)}.algn")
        mafft_cmd = f"mafft --inputorder \
                     --keeplength \
                     --compactmapout \
                     --anysymbol \
                     --kimura 1 \
                     --add {self.input_file} \
                     --6merpair \
                     --thread {self.threads} {self.reference} > {output}"
        run_mafft_cmd = subprocess.run(mafft_cmd, shell=True)
        self.fix_sequence_names(output, renamed_alignment)

        if self.mask_pos:
            self.mask_alignment(renamed_alignment)

        os.remove(output)
