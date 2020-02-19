import argparse

from gatk2ascat.parsers import parse_segments


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('--denoised-copy-ratios', required=True)
    parser.add_argument('--allelic-counts-tumor', required=True)
    parser.add_argument('--allelic-counts-normal', required=True)

    args = parser.parse_args()

    with open(args.denoised_copy_ratios, 'r') as f:
        segmentation = parse_segments(stream=f)
