import unittest

from gatk2ascat.core import Segment
from gatk2ascat.core import Segmentation
from gatk2ascat.core import SNP

from gatk2ascat.exceptions import UncoveredPositionError


class TestSegment(unittest.TestCase):

    def test_segment(self):
        segment = Segment(chromosome='chr2', start=100, end=200, logr=0.5)
        self.assertEqual(segment.chromosome, 'chr2')
        self.assertEqual(segment.start, 100)
        self.assertEqual(segment.end, 200)
        self.assertEqual(segment.logr, 0.5)


class TestSegmentation(unittest.TestCase):

    def setUp(self):
        self.segment_1a = Segment(chromosome='chr1', start=100, end=200, logr=0.5)
        self.segment_1b = Segment(chromosome='chr1', start=201, end=350, logr=1.2)
        self.segment_2a = Segment(chromosome='chr2', start=150, end=400, logr=-0.7)
        self.segmentation = Segmentation(segments=[self.segment_1a, self.segment_1b, self.segment_2a])

    def test_create_segmentation(self):
        self.assertEqual(self.segmentation._segments['chr1'], [self.segment_1a, self.segment_1b])
        self.assertEqual(self.segmentation._segments['chr2'], [self.segment_2a])

    def test_get_logr_for_covered_position(self):
        self.assertEqual(self.segmentation.logr(chromosome='chr1', position=100), 0.5)
        self.assertEqual(self.segmentation.logr(chromosome='chr1', position=200), 0.5)
        self.assertEqual(self.segmentation.logr(chromosome='chr1', position=300), 1.2)
        self.assertEqual(self.segmentation.logr(chromosome='chr2', position=400), -0.7)

    def test_get_logr_for_uncovered_position_raises_exception(self):
        with self.assertRaises(UncoveredPositionError):
            self.segmentation.logr(chromosome='chr2', position=100)


class TestSNP(unittest.TestCase):

    def test_snp(self):
        snp = SNP(chromosome='chr1', position=100, baf=0.45)
        self.assertEqual(snp.chromosome, 'chr1')
        self.assertEqual(snp.position, 100)
        self.assertEqual(snp.baf, 0.45)
