from typing import Iterable
from typing import List

from gatk2ascat.core import Segment
from gatk2ascat.core import Segmentation


def parse_segments(stream: Iterable[str]) -> Segmentation:

    segments: List[Segment] = []

    skipped_header = False
    for line in stream:
        if line.startswith('@'):
            continue
        if not skipped_header:
            skipped_header = True
            continue
        chromosome, start, end, logr = line.split('\t')
        segments.append(Segment(chromosome=chromosome,
                                start=int(start),
                                end=int(end),
                                logr=float(logr)))

    return Segmentation(segments=segments)
