import json
from dacite import from_dict
from gist.config import SubSamplingStateLineage, GetSimilarGenomes
from gist.error import InvalidInput
import os


def _valid_input_information(input_information: dict, configuration: object) -> object:
    return from_dict(data=input_information, data_class=configuration)


def read_get_states_input(input_file: str) -> object:
    with open(input_file, "r") as file_reader:
        validated_input = json.loads(file_reader.read())

        return _valid_input_information(validated_input, SubSamplingStateLineage)


def read_get_similar_genomes_input(input_file: str) -> object:
    with open(input_file, "r") as file_reader:
        validated_input = json.loads(file_reader.read())

        return _valid_input_information(validated_input, GetSimilarGenomes)


def check_dir(directory: str) -> Exception:
    if os.path.exists(directory) == False:
        raise InvalidInput(f"{directory} don't exist")
    if os.path.isdir(directory) == False:
        raise InvalidInput(f"{directory} isn't a valid directory")


def check_file(file: str) -> Exception:
    if os.path.exists(file) == False:
        raise InvalidInput(f"{file} don't exist")
    if os.path.isfile(file) == False:
        raise InvalidInput(f"{file} isn't a valid file")
