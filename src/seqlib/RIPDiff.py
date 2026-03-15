"""Framework for RIP-Seq differential enrichment analysis.

Provides skeletal classes and functions for comparing RNA Immunoprecipitation
(RIP) sequencing data against an isotype control (IgG) or total RNA input to
identify transcript segments preferentially enriched in the RIP sample.

RIP-Seq (RNA Immunoprecipitation followed by Sequencing) is used to identify
RNA molecules bound by a specific RNA-binding protein.

Note: This module is largely unimplemented (placeholder pass statements) and
is retained as a design scaffold for future development.
"""
##################
#Imports
##################
from . import intervallib

##################
#Classes
##################

class RIPUnit(intervallib.Interval):
    """A genomic interval unit used as the basic unit of RIP-Seq differential analysis.

    Can represent an individual transcript or any other genomic region (e.g. a
    whole chromosome) that is to be tested for differential read enrichment
    between a RIP sample and its control.  Extends intervallib.Interval.

    Note: All methods are currently unimplemented placeholders.
    """
    def __init__(self,interval):
        """Initialises a RIPUnit from an existing Interval instance.

        Args:
            interval: An intervallib.Interval object to copy coordinates from.

        Raises:
            AssertionError: If interval is not an instance of
                intervallib.Interval.
        """
        assert isinstance(interval,intervallib.Interval)
        intervallib.Interval.__init__(interval)

    def scan(self):
        """Scans the interval for differential RIP peaks (not implemented)."""
        pass

    def makebins(self,binSize):
        """Divides the interval into bins of the given size (not implemented).

        Args:
            binSize: Size of each bin in base pairs.
        """
        pass

    def binBinom(self):
        """Applies a binomial test to each bin (not implemented)."""
        pass

    def binPois(self):
        """Applies a Poisson test to each bin (not implemented)."""
        pass

    def fetchReads(self,bamHandle):
        """Fetches aligned reads overlapping this interval from a BAM file (not implemented).

        Args:
            bamHandle: A pysam AlignmentFile handle.
        """
        pass


#################
#Functions
#################
def globalNorm(ripUnit,totReads):
    """Applies global normalisation to a RIPUnit based on total library size (not implemented).

    Args:
        ripUnit: A RIPUnit object representing the region to normalise.
        totReads: Total number of mapped reads in the library, used as the
            normalisation denominator.
    """
    pass

def localNorm(ripUnitA,ripUnitB):
    """Applies local normalisation between two RIPUnit objects (not implemented).

    Args:
        ripUnitA: A RIPUnit from the experimental (RIP) sample.
        ripUnitB: A RIPUnit from the control (IgG or input) sample for the
            same genomic region.
    """
    pass
