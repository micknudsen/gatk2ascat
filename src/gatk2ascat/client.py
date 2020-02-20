import argparse

from typing import List
from typing import Optional

from gatk2ascat.core import BAF
from gatk2ascat.core import Segmentation

from gatk2ascat.core import generate_ascat_input

from gatk2ascat.parsers import get_sample_name
from gatk2ascat.parsers import parse_bafs
from gatk2ascat.parsers import parse_segments


def write_to_files(ascat_baf_file: str, ascat_logr_file: str, bafs: List[BAF], segmentation: Optional[Segmentation] = None):

    with open(ascat_baf_file, 'w') as baf_file, open(ascat_logr_file, 'w') as logr_file:

        header = '\t'.join(['', 'chromosome', 'position', 'SAMPLE_NAME'])
        print(header, file=baf_file)
        print(header, file=logr_file)

        for baf_entry, logr_entry in generate_ascat_input(bafs=bafs, segmentation=segmentation):
            print(baf_entry, file=baf_file)
            print(logr_entry, file=logr_file)


def main():

    parser = argparse.ArgumentParser()

    # Input files
    parser.add_argument('--denoised-copy-ratios', required=True)
    parser.add_argument('--allelic-counts-tumor', required=True)
    parser.add_argument('--allelic-counts-normal', required=True)

    # Output files
    parser.add_argument('--ascat-baf-tumor', required=True)
    parser.add_argument('--ascat-baf-normal', required=True)
    parser.add_argument('--ascat-logr-tumor', required=True)
    parser.add_argument('--ascat-logr-normal', required=True)

    args = parser.parse_args()

    with open(args.allelic_counts_tumor, 'r') as f:
        tumor_sample_name = get_sample_name(f)

    with open(args.allelic_counts_normal, 'r') as f:
        normal_sample_name = get_sample_name(f)

    with open(args.denoised_copy_ratios, 'r') as f:
        segmentation = parse_segments(stream=f)

    with open(args.allelic_counts_tumor, 'r') as f:
        tumor_bafs = parse_bafs(stream=f)

    with open(args.allelic_counts_normal, 'r') as f:
        normal_bafs = parse_bafs(stream=f)

    write_to_files(ascat_baf_file=args.ascat_baf_tumor, ascat_logr_file=args.ascat_logr_tumor, bafs=tumor_bafs segmentation=segmentation)
    write_to_files(ascat_baf_file=args.ascat_baf_normal, ascat_logr_file=args.ascat_logr_normal, bafs=normal_bafs)
