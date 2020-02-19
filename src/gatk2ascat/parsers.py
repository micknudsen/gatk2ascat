from typing import Iterator
from typing import List

from gatk2ascat.core import Segment
from gatk2ascat.core import Segmentation
from gatk2ascat.core import SNP


def skip_header(stream: Iterator[str]) -> None:
    """Expecst an iterator of lines of GATK output. Fast-forwards the
    iterator by skipping the SAM-style header comprising lines starting
    with @ and the following line, which contains column names."""
    for line in stream:
        # Skip the SAM-style header lines.
        if line.startswith('@'):
            continue
        # Now we have also skipped the column names line.
        break


def parse_segments(stream: Iterator[str]) -> Segmentation:
    """Parses output from GATK DenoiseReadCounts, which is a SAM-style
    header comprising lines starting with @ followed by single line
    with column names (CONTIG, START, STOP, and LOG2_COPY_RATIO)."""

    segments: List[Segment] = []

    skip_header(stream)

    for line in stream:
        chromosome, start, end, logr = line.split('\t')
        segment = Segment(chromosome=chromosome, start=int(start), end=int(end), logr=float(logr))
        segments.append(segment)

    return Segmentation(segments=segments)


def parse_snps(stream: Iterator[str]) -> List[SNP]:
    """Parses allelic counts output from GATK ModelSegments, which is a SAM-style
     header comprising lines starting with @ followed by single line with column
     names (CONTIG, POSITION, REF_COUNT, ALT_COUNT, REF_NUCLEOTIDE, ALT_NUCLEOTIDE)."""

    snps: List[SNP] = []

    skip_header(stream)

    for line in stream:
        chromosome, position, ref_count, alt_count, *_ = line.split('\t')
        baf = float(alt_count) / (float(ref_count) + float(alt_count))
        snp = SNP(chromosome=chromosome, position=int(position), baf=baf)
        snps.append(snp)

    return snps
