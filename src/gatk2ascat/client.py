import argparse


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('--denoised-copy-ratios', required=True)
    parser.add_argument('--allelic-counts-tumor', required=True)
    parser.add_argument('--allelic-counts-normal', required=True)

    args = parser.parse_args()
