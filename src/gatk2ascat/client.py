import argparse

from typing import List
from typing import Optional

from gatk2ascat.core import BAF
from gatk2ascat.core import Segmentation

from gatk2ascat.core import generate_ascat_input

from gatk2ascat.parsers import get_sample_name
from gatk2ascat.parsers import parse_bafs
from gatk2ascat.parsers import parse_segments


def sample_name_from_gatk_output(file: str) -> str:
    with open(file, 'r') as stream:
        return get_sample_name(stream=stream)


def bafs_from_gatk_output(file: str) -> List[BAF]:
    with open(file, 'r') as stream:
        return parse_bafs(stream=stream)


def write_to_files(ascat_baf_file: str, ascat_logr_file: str, bafs: List[BAF], sample_name: str, segmentation: Optional[Segmentation] = None) -> None:

    with open(ascat_baf_file, 'w') as baf_file, open(ascat_logr_file, 'w') as logr_file:

        header = '\t'.join(['', 'chromosome', 'position', sample_name])
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

    tumor_sample_name = sample_name_from_gatk_output(file=args.allelic_counts_tumor)
    normal_sample_name = sample_name_from_gatk_output(file=args.allelic_counts_normal)

    tumor_bafs = bafs_from_gatk_output(file=args.allelic_counts_tumor)
    normal_bafs = bafs_from_gatk_output(file=args.allelic_counts_normal)

    with open(args.denoised_copy_ratios, 'r') as stream:
        segmentation = parse_segments(stream=stream)

    write_to_files(ascat_baf_file=args.ascat_baf_tumor, ascat_logr_file=args.ascat_logr_tumor, bafs=tumor_bafs, sample_name=tumor_sample_name, segmentation=segmentation)
    write_to_files(ascat_baf_file=args.ascat_baf_normal, ascat_logr_file=args.ascat_logr_normal, bafs=normal_bafs, sample_name=normal_sample_name)
