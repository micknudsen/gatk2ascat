from typing import Iterable
from typing import List

from gatk2ascat.core import Segment
from gatk2ascat.core import Segmentation


def parse_segments(stream: Iterable[str]) -> Segmentation:
    """Parses output from GATK DenoiseReadCounts, which is a SAM-style
    header comprising lines starting with @ followed by single line
    with column names (CONTIG, START, STOP, and LOG2_COPY_RATIO)."""

    segments: List[Segment] = []

    # Make sure not to try to parse column names line.
    skipped_column_names_line = False

    for line in stream:

        # Skip the SAM-style header lines.
        if line.startswith('@'):
            continue

        # This first line after the SAM-style header is the
        # column names lines. Skip that one, too.
        if not skipped_column_names_line:
            skipped_column_names_line = True
            continue

        # All remaining lines correspond to a segment.
        chromosome, start, end, logr = line.split('\t')
        segment = Segment(chromosome=chromosome, start=int(start), end=int(end), logr=float(logr))
        segments.append(segment)

    return Segmentation(segments=segments)
