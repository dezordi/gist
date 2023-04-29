from dataclasses import dataclass, field
from typing import Optional, List
from gist.constants import BRAZILIAN_STATES
from datetime import datetime


@dataclass
class States:
    """
    Brazilian states information

    Attributes:
        name:
        max_genomes:
    """

    name: str
    sigla: str
    max_genomes: int

    def __init__(self, name: str, sigla:str, max_genomes: int) -> None:
        assert max_genomes > 0, "The state should have at least one genome"
        assert name in BRAZILIAN_STATES, f"{name} Invalid Brazilian state"
        self.name = name
        self.sigla = sigla
        self.max_genomes = max_genomes


@dataclass
class Countries:
    """
    Countries states information

    Attributes:
        name:
        max_genomes:
    """

    name: str
    sigla: str
    max_genomes: int

    def __init__(self, name: str, sigla: str, max_genomes: int) -> None:
        assert max_genomes > 0, "The country should have at least one genome"
        self.name = name
        self.sigla = sigla
        self.max_genomes = max_genomes


@dataclass
class SubSamplingStateLineage:
    """
    Information for subsampling

    Attributes:
        name:
        max_genomes:
    """

    job_name: str
    max_date: Optional[str]
    target_lineages: list
    states: List[States]
    countries: Optional[List[Countries]]
    min_date: str = "2019-12-26"
    min_genome_len: int = 28400
    outgroup_lineages: List[str] = field(default_factory=lambda: ["B.1"])

    def is_a_valid_date(self, date: str) -> bool:
        try:
            datetime.strptime(date, "%Y-%m-%d")
        except:
            return False

        return True

    def __init__(
        self,
        job_name: str,
        min_date: str,
        max_date: Optional[str],
        min_genome_len: int,
        target_lineages: list,
        outgroup_lineages: list,
        states: list,
        countries: list,
    ):
        assert (
            self.is_a_valid_date(min_date) == True
        ), "The min_date should follow YYYY-MM-DD pattern."
        if max_date:
            assert (
                self.is_a_valid_date(max_date) == True
            ), "The max_date should follow YYYY-MM-DD pattern."
        assert min_genome_len > 0, "The min_genome_len should be grather than 0."
        assert len(target_lineages) > 0, "At least one pango lineage should be parsed."
        assert len(states) > 0, "At least one brazilian state should be parsed"

        self.job_name = job_name
        self.min_date = min_date
        self.max_date = max_date
        self.min_genome_len = min_genome_len
        self.target_lineages = target_lineages
        self.outgroup_lineages = outgroup_lineages
        self.states = states
        self.countries = countries
