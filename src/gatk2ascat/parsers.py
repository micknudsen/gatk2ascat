from typing import Iterable
from typing import List

from gatk2ascat.core import Segment
from gatk2ascat.core import Segmentation


def parse_segments(stream: Iterable[str]) -> Segmentation:
    """Parses output from GATK DenoiseReadCounts, which is a SAM-style
    header comprising lines starting with @ followed by single line
    with column names (CONTIG, START, STOP, and LOG2_COPY_RATIO)."""

    segments: List[Segment] = []

    skipped_column_names_line = False
    for line in stream:
        if line.startswith('@'):
            continue
        if not skipped_column_names_line:
            skipped_column_names_line = True
            continue
        chromosome, start, end, logr = line.split('\t')
        segments.append(Segment(chromosome=chromosome,
                                start=int(start),
                                end=int(end),
                                logr=float(logr)))

    return Segmentation(segments=segments)
