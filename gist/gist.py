from gist.util import read_get_states_input, check_dir, check_file
import sys
import click
import os
import subprocess


class GetSubSamplingByState:
    def __init__(
        self, sequences: str, metadata: str, ncov_dir: str, subsampling_json_schema: str
    ) -> dict:
        try:
            for file in sequences, metadata, subsampling_json_schema:
                check_file(file)
            check_dir(ncov_dir)

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
        self.ncov_dir_data_job = os.path.join(
            ncov_dir, "data", valid_subsampling_schema.job_name
        )
        self.valid_subsampling_schema = valid_subsampling_schema
        self.sanitize_sequences()
        self.sanitize_metadata()
        self.get_samples()

    def _augur_index(self) -> None:
        augur_index_cmd = f"augur index \
                            --sequences {self.ncov_dir_data_job}/sequences_gisaid.fasta.gz \
                            --output {self.ncov_dir_data_job}/sequence_index_gisaid.tsv.gz"
        run_augur_index_cmd = subprocess.run(augur_index_cmd, shell=True)

    def sanitize_sequences(self) -> None:
        print("Sanitize and index sequences")
        if os.path.exists(self.ncov_dir_data_job) == False:
            print(f"Creating {self.ncov_dir_data_job}")
            os.mkdir(self.ncov_dir_data_job)
        sanitize_sequences_cmd = f"python {self.ncov_dir_scripts}/sanitize_sequences.py \
                                    --sequences {self.sequences} \
                                    --strip-prefixes hCoV-19/ SARS-CoV-2/ \
                                    --output {self.ncov_dir_data_job}/sequences_gisaid.fasta.gz"
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
                                    --output {self.ncov_dir_data_job}/metadata_gisaid.tsv.gz"
        run_sanitize_metadata_cmd = subprocess.run(sanitize_metadata_cmd, shell=True)
    
    def _filter_state(self, state: object) -> str:
        print(f"filter {state.name} genomes")
        output = os.path.join(self.ncov_dir_data_job,f"state_{state.sigla}_sub.txt")
        filter_by_state_cmd = f'''augur filter \
                                --metadata {self.ncov_dir_data_job}/metadata_gisaid.tsv.gz \
                                --sequence-index {self.ncov_dir_data_job}/sequence_index_gisaid.tsv.gz \
                                --min-date {self.valid_subsampling_schema.min_date} \
                                --max-date {self.valid_subsampling_schema.max_date} \
                                --min-length {str(self.valid_subsampling_schema.min_genome_len)} \
                                --query "(country == 'Brazil') & (division == '{state.name}') & (pango_lineage == {self.valid_subsampling_schema.target_lineages}) & (host == 'Human') " \
                                --subsample-max-sequences {str(state.max_genomes)} \
                                --exclude-ambiguous-dates-by any \
                                --group-by pango_lineage year month \
                                --output-strains {output}'''
        run_filter_by_state_cmd = subprocess.run(filter_by_state_cmd, shell=True)

        return output

    def _filter_country(self, country: object) -> str:
        print(f"filter {country.name} genomes")
        output = os.path.join(self.ncov_dir_data_job,f"country_{country.sigla}_sub.txt")
        filter_by_country_cmd = f'''augur filter \
                                --metadata {self.ncov_dir_data_job}/metadata_gisaid.tsv.gz \
                                --sequence-index {self.ncov_dir_data_job}/sequence_index_gisaid.tsv.gz \
                                --min-date {self.valid_subsampling_schema.min_date} \
                                --max-date {self.valid_subsampling_schema.max_date} \
                                --min-length {str(self.valid_subsampling_schema.min_genome_len)} \
                                --query "(country == '{country.name}') & (pango_lineage == {self.valid_subsampling_schema.target_lineages}) & (host == 'Human')" \
                                --subsample-max-sequences {str(country.max_genomes)} \
                                --exclude-ambiguous-dates-by any \
                                --group-by pango_lineage year month \
                                --output-strains {output}'''
        run_filter_by_country_cmd = subprocess.run(filter_by_country_cmd, shell=True)

        return output
    
    def _filter_outgroup(self) -> str:
        print(f"filter outgroup genomes")
        output = os.path.join(self.ncov_dir_data_job,"outgroup_sub.txt")
        filter_outgroup_cmd = f'''augur filter \
                                --metadata {self.ncov_dir_data_job}/metadata_gisaid.tsv.gz \
                                --sequence-index {self.ncov_dir_data_job}/sequence_index_gisaid.tsv.gz \
                                --max-date {self.valid_subsampling_schema.min_date} \
                                --min-length {str(self.valid_subsampling_schema.min_genome_len)} \
                                --query "pango_lineage == {self.valid_subsampling_schema.outgroup_lineages} & (host == 'Human')" \
                                --subsample-max-sequences 30 \
                                --exclude-ambiguous-dates-by any \
                                --group-by pango_lineage year month \
                                --output-strains {output}'''
        run_filter_outgroup_cmd = subprocess.run(filter_outgroup_cmd, shell=True)

        return output       

    def get_samples(self) -> None:
        sub_sampling_files = []
        for state in self.valid_subsampling_schema.states:
            sub_sampling_files.append(self._filter_state(state))
        if self.valid_subsampling_schema.countries:
            for country in self.valid_subsampling_schema.countries:
                sub_sampling_files.append(self._filter_country(country))
        sub_sampling_files.append(self._filter_outgroup())
        sub_sampling_files = ' '.join(sub_sampling_files)
        
        get_samples_cmd = f"augur filter \
                            --metadata {self.ncov_dir_data_job}/metadata_gisaid.tsv.gz \
                            --sequence-index {self.ncov_dir_data_job}/sequence_index_gisaid.tsv.gz \
                            --sequences {self.ncov_dir_data_job}/sequences_gisaid.fasta.gz \
                            --exclude-all \
                            --include {sub_sampling_files} \
                            --output-metadata {self.ncov_dir_data_job}/subsampled_metadata_gisaid.tsv \
                            --output-sequences {self.ncov_dir_data_job}/subsampled_sequences_gisaid.fasta"
        
        run_get_samples_cmd = subprocess.run(get_samples_cmd, shell=True)