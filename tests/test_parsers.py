import unittest

from gatk2ascat.core import BAF
from gatk2ascat.core import Segment

from gatk2ascat.parsers import parse_bafs
from gatk2ascat.parsers import parse_segments
from gatk2ascat.parsers import skip_header


class TestParsers(unittest.TestCase):

    def test_skip_header(self):

        stream = iter(['\t'.join(['@HD', 'VN:1.6']),
                       '\t'.join(['@SQ', 'SN:chr1 LN:248956422']),
                       '\t'.join(['@SQ', 'SN:chr2 LN:242193529']),
                       '\t'.join(['@RG', 'ID:GATKCopyNumber', 'SM:TESTSAMPLE']),
                       '\t'.join(['CONTIG', 'START', 'END', 'LOG2_COPY_RATIO']),
                       '\t'.join(['chr1', '100', '200', '0.1']),
                       '\t'.join(['chr1', '500', '750', '3.4']),
                       '\t'.join(['chr2', '100', '200', '-1.7'])])

        skip_header(stream)
        self.assertEqual(list(stream), ['\t'.join(['chr1', '100', '200', '0.1']),
                                        '\t'.join(['chr1', '500', '750', '3.4']),
                                        '\t'.join(['chr2', '100', '200', '-1.7'])])

    def test_parse_segments(self):

        stream = iter(['\t'.join(['@HD', 'VN:1.6']),
                       '\t'.join(['@SQ', 'SN:chr1 LN:248956422']),
                       '\t'.join(['@SQ', 'SN:chr2 LN:242193529']),
                       '\t'.join(['@RG', 'ID:GATKCopyNumber', 'SM:TESTSAMPLE']),
                       '\t'.join(['CONTIG', 'START', 'END', 'LOG2_COPY_RATIO']),
                       '\t'.join(['chr1', '100', '200', '0.1']),
                       '\t'.join(['chr1', '500', '750', '3.4']),
                       '\t'.join(['chr2', '100', '200', '-1.7'])])

        segementation = parse_segments(stream=stream)

        self.assertEqual(segementation._segments['chr1'],
                         [Segment(chromosome='chr1', start=100, end=200, logr=0.1),
                          Segment(chromosome='chr1', start=500, end=750, logr=3.4)])

        self.assertEqual(segementation._segments['chr2'],
                         [Segment(chromosome='chr2', start=100, end=200, logr=-1.7)])

    def test_baf_parser(self):

        stream = iter(['\t'.join(['@HD', 'VN:1.6']),
                       '\t'.join(['@SQ', 'SN:chr1 LN:248956422']),
                       '\t'.join(['@SQ', 'SN:chr2 LN:242193529']),
                       '\t'.join(['@RG', 'ID:GATKCopyNumber', 'SM:TESTSAMPLE']),
                       '\t'.join(['CONTIG', 'POSITION', 'REF_COUNT', 'ALT_COUNT', 'REF_NUCLEOTIDE', 'ALT_NUCLEOTIDE']),
                       '\t'.join(['chr1', '150', '10', '90', 'A', 'G']),
                       '\t'.join(['chr1', '700', '0', '30', 'T', 'C']),
                       '\t'.join(['chr2', '100', '30', '50', 'G', 'A'])])

        bafs = parse_bafs(stream=stream)

        self.assertEqual(bafs, [BAF(chromosome='chr1', position=150, frequency=0.9),
                                BAF(chromosome='chr1', position=700, frequency=1.0),
                                BAF(chromosome='chr2', position=100, frequency=0.625)])
