import unittest

from gatk2ascat.core import Segment
from gatk2ascat.core import Segmentation


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

    def test_get_logr(self):
        self.assertEqual(self.segmentation.logr(chromosome='chr1', position='100'), 0.5)
        self.assertEqual(self.segmentation.logr(chromosome='chr1', position='200'), 0.5)
        self.assertEqual(self.segmentation.logr(chromosome='chr1', position='300'), 1.2)
        self.assertEqual(self.segmentation.logr(chromosome='chr2', position='400'), -0.7)
