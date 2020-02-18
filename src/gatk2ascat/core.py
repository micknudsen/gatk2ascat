from collections import defaultdict

from typing import DefaultDict
from typing import Iterable
from typing import List
from typing import NamedTuple


class Segment(NamedTuple):
    chromosome: str
    start: int
    end: int
    logr: float


class Segmentation:

    def __init__(self, segments: Iterable[Segment]) -> None:
        self._segments: DefaultDict[str, List[Segment]] = defaultdict(list)
        for segment in segments:
            self._segments[segment.chromosome].append(segment)
