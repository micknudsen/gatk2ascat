import argparse

from gatk2ascat.core import generate_ascat_input

from gatk2ascat.parsers import parse_bafs
from gatk2ascat.parsers import parse_segments


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

    with open(args.denoised_copy_ratios, 'r') as f:
        segmentation = parse_segments(stream=f)

    with open(args.allelic_counts_tumor, 'r') as f:
        tumor_bafs = parse_bafs(stream=f)

    with open(args.allelic_counts_normal, 'r') as f:
        normal_bafs = parse_bafs(stream=f)

    with open(args.ascat_baf_tumor, 'w') as baf_file, open(args.ascat_logr_tumor, 'w') as logr_file:

        print('', 'chromosome', 'position', 'tumor', sep='\t', file=baf_file)
        print('', 'chromosome', 'position', 'tumor', sep='\t', file=logr_file)

        for baf_entry, logr_entry in generate_ascat_input(bafs=tumor_bafs, segmentation=segmentation):
            print(baf_entry, file=baf_file)
            print(logr_entry, file=logr_file)

    with open(args.ascat_baf_normal, 'w') as baf_file, open(args.ascat_logr_normal, 'w') as logr_file:

        print('', 'chromosome', 'position', 'tumor', sep='\t', file=baf_file)
        print('', 'chromosome', 'position', 'tumor', sep='\t', file=logr_file)

        for baf_entry, logr_entry in generate_ascat_input(bafs=normal_bafs):
            print(baf_entry, file=baf_file)
            print(logr_entry, file=logr_file)
