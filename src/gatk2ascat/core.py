from collections import defaultdict

from typing import DefaultDict
from typing import Iterable
from typing import Iterator
from typing import List
from typing import NamedTuple
from typing import Optional
from typing import Tuple

from gatk2ascat.exceptions import UncoveredPositionError


class BAF(NamedTuple):
    chromosome: str
    position: int
    ref_count: int
    alt_count: int
    ref_nucleotide: str
    alt_nucleotide: str

    @property
    def frequency(self):
        return self.alt_count / (self.ref_count + self.alt_count)


class Segment(NamedTuple):
    chromosome: str
    start: int
    end: int
    logr: float


class ASCATDataPoint(NamedTuple):
    chromosome: str
    position: int
    value: float

    @property
    def name(self) -> str:
        return self.chromosome + '_' + str(self.position)

    def __str__(self) -> str:
        return '\t'.join([f'{self.chromosome}_{self.position}', self.chromosome, str(self.position), str(self.value)])


class Segmentation:

    def __init__(self, segments: Iterable[Segment]) -> None:
        self._segments: DefaultDict[str, List[Segment]] = defaultdict(list)
        for segment in segments:
            self._segments[segment.chromosome].append(segment)

    def logr(self, chromosome: str, position: int) -> float:
        for segment in self._segments[chromosome]:
            if segment.start <= position <= segment.end:
                return segment.logr
        raise UncoveredPositionError


def get_consensus_bafs(tumor_bafs: List[BAF], normal_bafs: List[BAF]) -> Iterator[Tuple[BAF, BAF]]:

    for tumor_baf, normal_baf in zip(tumor_bafs, normal_bafs):
        if all([tumor_baf.ref_nucleotide == normal_baf.ref_nucleotide,
                tumor_baf.alt_nucleotide == normal_baf.alt_nucleotide]):
            yield tumor_baf, normal_baf


def generate_ascat_input(bafs: List[BAF], segmentation: Optional[Segmentation] = None) -> Iterator[Tuple[ASCATDataPoint, ASCATDataPoint]]:
    """Takes as input a list of BAF objects and an optional Segmentation object. Yields pairs
    of ASCATDataPoint objects with BAFs and LOGRs. If a Segmentation object is provided, LOGR
    values are determined from that. Otherwise, all LOGR value are set to 0.0."""

    for baf in bafs:

        if segmentation:
            logr = segmentation.logr(chromosome=baf.chromosome, position=baf.position)
        else:
            logr = 0.0

        baf_data_point = ASCATDataPoint(chromosome=baf.chromosome, position=baf.position, value=baf.frequency)
        logr_data_point = ASCATDataPoint(chromosome=baf.chromosome, position=baf.position, value=logr)

        yield baf_data_point, logr_data_point
