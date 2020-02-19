import argparse

from gatk2ascat.parsers import parse_segments
from gatk2ascat.parsers import parse_snps


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('--denoised-copy-ratios', required=True)
    parser.add_argument('--allelic-counts-tumor', required=True)
    parser.add_argument('--allelic-counts-normal', required=True)

    args = parser.parse_args()

    with open(args.denoised_copy_ratios, 'r') as f:
        segmentation = parse_segments(stream=f)

    with open(args.allelic_counts_tumor, 'r') as f:
        tumor_baf = parse_snps(stream=f)

    with open(args.allelic_counts_normal, 'r') as f:
        normal_baf = parse_snps(stream=f)

    # Just to silence flake8 while developing!
    print(segmentation)
    print(tumor_baf)
    print(normal_baf)
