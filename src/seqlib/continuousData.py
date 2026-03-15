"""High-resolution genome-wide continuous data storage structures.

Provides ContinuousData for per-nucleotide or binned coverage arrays on a
single chromosome, and SimpleChIPData for loading, normalising, and scanning
multi-sample NimbleGen ChIP data.

Note: SimpleChIPData depends on rpy2 and tables (PyTables) as well as the
Chip module from this package.

First attempt at a data structure for high-resolution genome-wide data.

Originally created on Jun 30, 2009.

Author: lgoff
"""
import gzip
import sys

import numpy as np
import rpy2.robjects as rpy
from tables import *

from . import Chip, genomelib


class ContinuousData(object):
    """Per-chromosome continuous (coverage) data storage backed by numpy arrays.

    Stores strand-separated floating-point data at a configurable bin
    resolution. Supports interval-based data accumulation, range extraction,
    and gzipped binary serialisation.

    Attributes:
        name: Sample name string.
        chr: Chromosome name (must be in genomelib.chr_lengths).
        binSize: Resolution in base pairs per bin (default 1).
        fname: Default filename for binary output.
        data: Dict with "+" and "-" keys mapping to numpy float64 arrays.
    """

    def __init__(self,name,chr,binSize = 1,data = {}):
        """Construct a ContinuousData object for a single chromosome.

        If data is non-empty, it is used directly. Otherwise, two zero-filled
        numpy arrays of length chr_length // binSize are created.

        Args:
            name: Sample name string.
            chr: Chromosome name string (must be in genomelib.chr_lengths).
            binSize: Bin size in base pairs (default 1).
            data: Optional pre-existing dict with "+" and "-" numpy arrays.
        """
        self.name = name
        self.chr = chr
        self.binSize = int(binSize)
        self.fname = "%s_%s_%d.bin" % (self.name,self.chr,self.binSize)
        if len(data.keys())>0:
            self.data = data
        else:
            self.data = {
                         '+':np.zeros(genomelib.chr_lengths[chr]//binSize,'d'),
                         '-':np.zeros(genomelib.chr_lengths[chr]//binSize ,'d')
                         }

    def __len__(self):
        """Return the number of bins, equivalent to the chromosome length in bins."""
        return np.alen(self.data['+'])

    def __repr__(self):
        """Return the sample name string."""
        return self.name

    def __str__(self):
        """Return the sample name string."""
        return self.name

    def getMin(self,strand):
        """Return the minimum value in the data array for the given strand.

        Args:
            strand: "+" or "-".

        Returns:
            Minimum float value in self.data[strand].
        """
        return np.amin(self.data[strand])

    def getMax(self,strand):
        """Return the maximum value in the data array for the given strand.

        Args:
            strand: "+" or "-".

        Returns:
            Maximum float value in self.data[strand].
        """
        return np.amax(self.data[strand])

    def whichMax(self,strand):
        """Return the bin index of the maximum value for the given strand.

        Args:
            strand: "+" or "-".

        Returns:
            Integer index of the maximum element in self.data[strand].
        """
        return np.argmax(self.data[strand])

    def whichMin(self,strand):
        """Return the bin index of the minimum value for the given strand.

        Args:
            strand: "+" or "-".

        Returns:
            Integer index of the minimum element in self.data[strand].
        """
        return np.argmin(self.data[strand])

    def getDataRange(self,strand,start,end):
        """Return the data array slice corresponding to a genomic coordinate range.

        Args:
            strand: "+" or "-".
            start: Genomic start coordinate.
            end: Genomic end coordinate.

        Returns:
            Numpy array slice of self.data[strand] for the given range.
        """
        return self.data[strand][(start//self.binSize)-1:(end//self.binSize)-1]

    def addInterval(self,interval):
        """Accumulate an interval's count into the data arrays.

        Adds interval.count to each bin covered by the interval on its strand.
        Does nothing if the interval's chromosome does not match self.chr.

        Args:
            interval: An object with chr, strand, start, end, and count
                attributes.

        Returns:
            The string "Wrong data file" if interval.chr != self.chr,
            otherwise None.
        """
        if self.chr != interval.chr:
            return "Wrong data file"
        else:
            self.data[interval.strand][(interval.start//self.binSize)-1:(interval.end//self.binSize)-1]=self.data[interval.strand][(interval.start//self.binSize)-1:(interval.end//self.binSize)-1]+interval.count

    def write(self,fname=None):
        """Write data arrays to a gzipped binary file.

        Args:
            fname: Output file path. Defaults to self.fname if not provided.
        """
        if fname == None:
            fname = self.fname
        fd = gzip.open(fname,'wb')
        for s in self.data.keys():
            fd.write(self.data[s])
        fd.close()

    def read(self,fname):
        """Read data from a file (not yet implemented).

        Args:
            fname: Path to the file to read from.
        """
        pass

    def innerHeight(self,strand,start,end):
        """Return the maximum value (peak height) within a genomic range.

        Args:
            strand: "+" or "-".
            start: Genomic start coordinate.
            end: Genomic end coordinate.

        Returns:
            Maximum float value in the data range.
        """
        region = self.getDataRange(strand,start,end)
        return np.amax(region)

    def outerHeight(self,strand,start,end):
        """Return the total signal (sum) within a genomic range.

        Args:
            strand: "+" or "-".
            start: Genomic start coordinate.
            end: Genomic end coordinate.

        Returns:
            Sum of all values in the data range.
        """
        region = self.getDataRange(strand,start,end)
        return sum(region)

class SimpleChIPData(object):
    """Multi-sample NimbleGen ChIP-chip data container with normalisation and scanning.

    Loads NimbleGen GFF probe files, applies quantile normalisation via
    limma, joins probes into intervals, and scans intervals with a sliding
    window test.

    Attributes:
        data: Dict mapping sample name to list of probe Intervals.
        samples: List of sample name strings in load order.
        dataMatrix: 2D numpy float array of probe scores (set by makeMatrix).
        normMatrix: 2D numpy array of quantile-normalised scores (set by
            quantileNormalize).
        intervals: Dict mapping sample name to list of joined Intervals (set
            by joinProbes).
    """

    def __init__(self,files):
        """Load NimbleGen GFF files and initialise the data store.

        Args:
            files: List of GFF file paths to load. Each file's sample name is
                derived by stripping the ".gff" extension.
        """
        self.data = {}
        self.samples = []
        for fname in files:
            sampleName = fname.rstrip(".gff")
            self.samples.append(sampleName)
            sys.stderr.write("Parsing file '%s'...\n" % fname)
            self.data[sampleName] = Chip.parseNimblegen(fname)

    def doIt(self,permuted,windows=[5,6,7,8,9,10,11,12],threshold=0.05):
        """Run the full normalise-join-scan pipeline.

        Calls normalize(), joinProbes(), and then scan() for each window size.

        Args:
            permuted: Permuted score data passed to scan() for significance
                testing.
            windows: List of window sizes to scan (default [5..12]).
            threshold: Significance threshold for scanning (default 0.05).
        """
        self.normalize()
        self.joinProbes()
        for winSize in windows:
            self.scan(permuted,winSize,threshold)

    def makeMatrix(self):
        """Build self.dataMatrix from probe scores across all samples.

        Creates a 2D numpy float array of shape (n_probes, n_samples) where
        each column contains the scores for one sample in probe order.
        Writes a progress message to stderr on completion.
        """
        data_keys = list(self.data.keys())
        self.dataMatrix = np.empty((len(self.data[data_keys[0]]),len(self.samples)),'f')
        for i in range(0,len(data_keys)):
            self.dataMatrix[:,i]=[x.score for x in self.data[data_keys[i]]]
        sys.stderr.write("Created dataMatrix!\n")

    def quantileNormalize(self):
        """Apply quantile normalisation to self.dataMatrix using limma.

        Calls makeMatrix() first if dataMatrix is not yet set. Requires the
        R limma package. Stores the result in self.normMatrix.
        """
        if 'dataMatrix' not in self.__dict__: self.makeMatrix()
        rpy.r.library("limma")
        sys.stderr.write("Performing Quantile Normalization...\n")
        self.normMatrix = rpy.r.normalizeQuantiles(self.dataMatrix)

    def normalize(self):
        """Replace probe scores with quantile-normalised values.

        Calls quantileNormalize() first if normMatrix is not yet set. Updates
        the score attribute of every probe object in self.data in-place.
        """
        if 'normMatrix' not in self.__dict__: self.quantileNormalize()
        sys.stderr.write("Replacing values in data with normalized values...\n")
        data_keys = list(self.data.keys())
        for i in range(0,len(data_keys)):
            for j in range(0,np.shape(self.normMatrix)[0]):
                self.data[data_keys[i]][j].score = self.normMatrix[j,i]

    def joinProbes(self):
        """Join adjacent probes into contiguous intervals for each sample.

        Populates self.intervals dict via Chip.joinNimblegenIntervals().
        Writes per-sample progress messages to stderr.
        """
        sys.stderr.write("Joining Probes into intervals...\n")
        self.intervals = {}
        for sample in self.samples:
            sys.stderr.write("\t%s\n" % sample)
            self.intervals[sample] = Chip.joinNimblegenIntervals(self.data[sample])

    def scan(self,permuted,windowSize,threshold=0.05):
        """Scan all intervals with a sliding window test of the given size.

        Calls i.scan(permuted, windowSize, threshold) on every interval in
        every sample. Writes progress messages to stderr.

        Args:
            permuted: Permuted score data used by the interval scan method for
                significance testing.
            windowSize: Integer number of probes per sliding window.
            threshold: Significance threshold (default 0.05).
        """
        sys.stderr.write("Scanning with window of size %d..\n" % windowSize)
        for sample in self.samples:
            sys.stderr.write("\t%s\n" % sample)
            for i in self.intervals[sample]:
                i.scan(permuted,windowSize,threshold)
