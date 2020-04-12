[![install with conda](https://img.shields.io/badge/install%20with-conda-brightgreen.svg?style=flat)](https://conda.anaconda.org/micknudsen) ![CI](https://github.com/micknudsen/gatk2ascat/workflows/CI/badge.svg?branch=master) [![Coverage Status](https://coveralls.io/repos/github/micknudsen/gatk2ascat/badge.svg?branch=master)](https://coveralls.io/github/micknudsen/gatk2ascat?branch=master)

# gatk2ascat

Just a simple tool for converting output from the [GATK Somatic Copy Number Workflow](https://gatk.broadinstitute.org/hc/en-us/articles/360035531092?id=11682) to something which can be used as input for [ASCAT](https://www.crick.ac.uk/research/labs/peter-van-loo/software).

## Input

Three files produced by GATK are needed as input. These are the allelic counts at heterozygous sites for the tumor and normal samples, respectively, as well as the denoised copy number ratios. In the official GATK workflow, these files are named `SAMPLE.hets.tsv`, `SAMPLE.hets.normal.tsv`, and `SAMPLE.denoisedCR.tsv`.

## Running

The `gatk2ascat` command is then run as follows:

```
gatk2ascat --denoised-copy-ratios SAMPLE.denoisedCR.tsv
           --allelic-counts-tumor SAMPLE.hets.tsv
           --allelic-counts-normal SAMPLE.hets.normal.tsv
           --ascat-baf-tumor Tumor_BAF.txt
           --ascat-baf-normal Germline_BAF.txt
           --ascat-logr-tumor Tumor_LogR.txt
           --ascat-logr-normal Germline_LogR.txt
```

and produces four output files (`Tumor_BAF.txt`, `Germline_BAF.txt`, `Tumor_LogR.txt`, and `Germline_LogR.txt`) ready for use as input for ASCAT.

## Install

The simplest way to install `gatk2ascat` is by using conda:

```$ conda install -c micknudsen gatk2ascat```
