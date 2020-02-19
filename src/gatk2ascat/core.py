from collections import defaultdict

from typing import DefaultDict
from typing import Iterable
from typing import List
from typing import NamedTuple

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

    def name(self) -> str:
        return self.chromosome + '_' + str(self.position)


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
