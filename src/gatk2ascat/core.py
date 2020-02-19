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
    frequency: float


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


def generate_ascat_input(bafs: List[BAF], segmentation: Optional[Segmentation] = None) -> Iterator[Tuple[ASCATDataPoint, ASCATDataPoint]]:

    for baf in bafs:

        logr = segmentation.logr(chromosome=baf.chromosome, position=baf.position) if segmentation else 0.0

        yield ASCATDataPoint(chromosome=baf.chromosome, position=baf.position, value=baf.frequency), ASCATDataPoint(chromosome=baf.chromosome, position=baf.position, value=logr)
