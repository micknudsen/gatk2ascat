import argparse

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
        for baf in tumor_bafs:
            print(f'{baf.chromosome}_{baf.position}', baf.chromosome, baf.position, baf.frequency, sep='\t', file=baf_file)
            print(f'{baf.chromosome}_{baf.position}', baf.chromosome, baf.position, segmentation.logr(chromosome=baf.chromosome, position=baf.position), sep='\t', file=logr_file)

    with open(args.ascat_baf_normal, 'w') as baf_file, open(args.ascat_logr_normal, 'w') as logr_file:

        print('', 'chromosome', 'position', 'tumor', sep='\t', file=baf_file)
        for baf in normal_bafs:
            print(f'{baf.chromosome}_{baf.position}', baf.chromosome, baf.position, baf.frequency, sep='\t', file=baf_file)
            print(f'{baf.chromosome}_{baf.position}', baf.chromosome, baf.position, '0.0', sep='\t', file=logr_file)
