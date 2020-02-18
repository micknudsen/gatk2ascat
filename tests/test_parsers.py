import unittest

from gatk2ascat.core import Segment
from gatk2ascat.parsers import parse_segments


class TestSegmentsParser(unittest.TestCase):

    def test_parse_segments(self):

        stream = ['\t'.join(['@HD', 'VN:1.6']),
                  '\t'.join(['@SQ', 'SN:chr1 LN:248956422']),
                  '\t'.join(['@SQ', 'SN:chr2 LN:242193529']),
                  '\t'.join(['@RG', 'ID:GATKCopyNumber', 'SM:TESTSAMPLE']),
                  '\t'.join(['CONTIG', 'START', 'END', 'LOG2_COPY_RATIO']),
                  '\t'.join(['chr1', '100', '200', '0.1']),
                  '\t'.join(['chr1', '500', '750', '3.4']),
                  '\t'.join(['chr2', '100', '200', '-1.7'])]

        segementation = parse_segments(stream=stream)

        self.assertEqual(segementation._segments['chr1'],
                         [Segment(chromosome='chr1', start=100, end=200, logr=0.1),
                          Segment(chromosome='chr1', start=500, end=750, logr=3.4)])

        self.assertEqual(segementation._segments['chr2'],
                         [Segment(chromosome='chr2', start=100, end=200, logr=-1.7)])
