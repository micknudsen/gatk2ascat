import unittest

from gatk2ascat.core import BAF
from gatk2ascat.core import ASCATDataPoint
from gatk2ascat.core import Segment
from gatk2ascat.core import Segmentation

from gatk2ascat.core import generate_ascat_input
from gatk2ascat.core import get_consensus_bafs

from gatk2ascat.exceptions import UncoveredPositionError


class TestSimpleStructures(unittest.TestCase):

    def test_baf(self):
        baf = BAF(chromosome='chr1', position=100, ref_count=55, alt_count=45, ref_nucleotide='A', alt_nucleotide='T')
        self.assertEqual(baf.chromosome, 'chr1')
        self.assertEqual(baf.position, 100)
        self.assertEqual(baf.ref_count, 55)
        self.assertEqual(baf.alt_count, 45)
        self.assertEqual(baf.ref_nucleotide, 'A')
        self.assertEqual(baf.alt_nucleotide, 'T')
        self.assertEqual(baf.frequency, 0.45)

    def test_segment(self):
        segment = Segment(chromosome='chr2', start=100, end=200, logr=0.5)
        self.assertEqual(segment.chromosome, 'chr2')
        self.assertEqual(segment.start, 100)
        self.assertEqual(segment.end, 200)
        self.assertEqual(segment.logr, 0.5)

    def test_ascat_data_point(self):
        data_point = ASCATDataPoint(chromosome='chr1', position=100, value=3.2)
        self.assertEqual(data_point.chromosome, 'chr1')
        self.assertEqual(data_point.position, 100)
        self.assertEqual(data_point.value, 3.2)
        self.assertEqual(data_point.name, 'chr1_100')

    def test_ascat_data_point_string_representation(self):
        data_point = ASCATDataPoint(chromosome='chr1', position=100, value=3.2)
        self.assertEqual(data_point.__str__(), '\t'.join(['chr1_100', 'chr1', '100', '3.2']))


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


class TestConsensusBAFs(unittest.TestCase):

    def setUp(self):

        self.tumor_bafs = [BAF(chromosome='chr1', position=150, ref_count=80, alt_count=20, ref_nucleotide='T', alt_nucleotide='G'),
                           BAF(chromosome='chr1', position=175, ref_count=85, alt_count=15, ref_nucleotide='G', alt_nucleotide='C'),
                           BAF(chromosome='chr2', position=300, ref_count=10, alt_count=90, ref_nucleotide='C', alt_nucleotide='G')]

        self.normal_bafs = [BAF(chromosome='chr1', position=150, ref_count=50, alt_count=50, ref_nucleotide='T', alt_nucleotide='A'),
                            BAF(chromosome='chr1', position=175, ref_count=45, alt_count=55, ref_nucleotide='G', alt_nucleotide='C'),
                            BAF(chromosome='chr2', position=300, ref_count=57, alt_count=43, ref_nucleotide='A', alt_nucleotide='G')]

    def test_consensus_bafs(self):

        consensus_tumor_bafs, consensus_normal_bafs = get_consensus_bafs(tumor_bafs=self.tumor_bafs, normal_bafs=self.normal_bafs)

        self.assertEqual(list(consensus_tumor_bafs), [BAF(chromosome='chr1', position=175, ref_count=85, alt_count=15, ref_nucleotide='G', alt_nucleotide='C')])
        self.assertEqual(list(consensus_normal_bafs), [BAF(chromosome='chr1', position=175, ref_count=45, alt_count=55, ref_nucleotide='G', alt_nucleotide='C')])


class TestOutputGenerator(unittest.TestCase):

    def setUp(self):

        self.segmentation = Segmentation(segments=[Segment(chromosome='chr1', start=100, end=200, logr=0.5),
                                                   Segment(chromosome='chr1', start=201, end=350, logr=1.2),
                                                   Segment(chromosome='chr2', start=150, end=400, logr=-0.7)])

        self.tumor_bafs = [BAF(chromosome='chr1', position=150, ref_count=80, alt_count=20, ref_nucleotide='T', alt_nucleotide='A'),
                           BAF(chromosome='chr1', position=175, ref_count=85, alt_count=15, ref_nucleotide='G', alt_nucleotide='C'),
                           BAF(chromosome='chr2', position=300, ref_count=10, alt_count=90, ref_nucleotide='C', alt_nucleotide='G')]

        self.normal_bafs = [BAF(chromosome='chr1', position=150, ref_count=50, alt_count=50, ref_nucleotide='T', alt_nucleotide='A'),
                            BAF(chromosome='chr1', position=175, ref_count=45, alt_count=55, ref_nucleotide='G', alt_nucleotide='C'),
                            BAF(chromosome='chr2', position=300, ref_count=57, alt_count=43, ref_nucleotide='C', alt_nucleotide='G')]

    def test_generate_ascat_input_with_segmentation(self):

        ascat_baf, ascat_logr = zip(*generate_ascat_input(bafs=self.tumor_bafs, segmentation=self.segmentation))

        self.assertEqual(ascat_baf, (ASCATDataPoint(chromosome='chr1', position=150, value=0.2),
                                     ASCATDataPoint(chromosome='chr1', position=175, value=0.15),
                                     ASCATDataPoint(chromosome='chr2', position=300, value=0.9)))

        self.assertEqual(ascat_logr, (ASCATDataPoint(chromosome='chr1', position=150, value=0.5),
                                      ASCATDataPoint(chromosome='chr1', position=175, value=0.5),
                                      ASCATDataPoint(chromosome='chr2', position=300, value=-0.7)))

    def test_generate_ascat_input_without_segmentation(self):

        ascat_baf, ascat_logr = zip(*generate_ascat_input(bafs=self.normal_bafs))

        self.assertEqual(ascat_baf, (ASCATDataPoint(chromosome='chr1', position=150, value=0.5),
                                     ASCATDataPoint(chromosome='chr1', position=175, value=0.55),
                                     ASCATDataPoint(chromosome='chr2', position=300, value=0.43)))

        self.assertEqual(ascat_logr, (ASCATDataPoint(chromosome='chr1', position=150, value=0.0),
                                      ASCATDataPoint(chromosome='chr1', position=175, value=0.0),
                                      ASCATDataPoint(chromosome='chr2', position=300, value=-0.0)))
