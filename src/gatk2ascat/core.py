from typing import Iterable
from typing import NamedTuple


class Segment(NamedTuple):
    chromosome: str
    start: int
    end: int
    logr: float


class Segmentation:

    def __init__(self, segments: Iterable[Segment]) -> None:
        pass
