from typing import NamedTuple


class Segment(NamedTuple):
    chromosome: str
    start: int
    end: int
    logr: float
