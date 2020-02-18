import unittest

from gatk2ascat.core import Segment


class TestSegment(unittest.TestCase):

    def test_segment(self):
        segment = Segment(chromosome='chr2', start=100, end=200, logr=0.5)
        self.assertEqual(segment.chromosome, 'chr2')
        self.assertEqual(segment.start, 100)
        self.assertEqual(segment.end, 200)
        self.assertEqual(segment.logr, 0.5)
