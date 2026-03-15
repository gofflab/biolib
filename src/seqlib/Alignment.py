"""Short RNA read alignment data structure.

Provides the Alignment class for representing a single short-read alignment,
with methods for strand testing, BED output, and conversion to intervallib
Interval objects.

Originally created on Jun 30, 2009.

Author: lgoff
"""
from . import misc
from .intervallib import *


class Alignment(object):
    """Basic alignment class for short RNA reads.

    Can be bypassed in favour of aligner-specific implementations such as
    ShrimpRead or MAQRead. Supports score-based sorting (higher scores sort
    first) and conversion to BED or Interval format.

    Attributes:
        readname: Name/identifier of the aligned read.
        chr: Chromosome name.
        start: 0-based start coordinate.
        end: End coordinate.
        strand: Strand orientation ("+" or "-").
        score: Alignment score (float).
        readsequence: DNA sequence of the read.
        readcount: Integer read count (-1 if unset).
    """
    def __init__(self,readname,chr,start,end,strand,score=0,readcount = -1,readsequence=''):
        """Initialize an Alignment.

        Args:
            readname: Name/identifier of the read.
            chr: Chromosome name string.
            start: Start coordinate (converted to int).
            end: End coordinate (converted to int).
            strand: Strand string ("+" or "-").
            score: Alignment score (default 0, converted to float).
            readcount: Read count integer (default -1).
            readsequence: DNA sequence string of the read (default "").
        """
        self.readname = str(readname)
        self.chr = chr
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.score = float(score)
        self.readsequence = readsequence
        self.readcount = readcount

    def __lt__(self, b):
        """Compare by score in descending order (higher scores sort first)."""
        return self.score > b.score  # reversed because original was -cmp(self.score, b.score)

    def __eq__(self, b):
        """Return True if self and b have the same score."""
        return self.score == b.score

    def __str__(self):
        """Return a readname:chr:start:end string."""
        return "%s:%s:%d:%d" % (self.readname,self.chr,self.start,self.end)

    def __repr__(self):
        """Return a readname:chr:start:end string."""
        return "%s:%s:%d:%d" % (self.readname,self.chr,self.start,self.end)

    def __len__(self):
        """Return the length of the alignment in bases (end - start + 1)."""
        return self.end-self.start+1

    def isPlus(self):
        """Return True if the alignment is on the "+" strand.

        Returns:
            True if self.strand == "+", otherwise False.
        """
        if self.strand=="+":
            return True
        else:
            return False

    def isMinus(self):
        """Return True if the alignment is on the "-" strand.

        Returns:
            True if self.strand == "-", otherwise False.
        """
        if self.strand=="-":
            return True
        else:
            return False

    def toInterval(self):
        """Convert this alignment to an intervallib.Interval.

        Returns:
            An Interval with the same coordinates, score, readcount, and
            readname as this alignment.
        """
        return Interval(self.chr,self.start,self.end,self.strand,self.score,self.readcount,name=self.readname)

    def toBed(self):
        """Return a BED-formatted string for this alignment.

        The name field is encoded using misc.seq2nuID applied to the read
        sequence.

        Returns:
            Tab-delimited BED line string with a trailing newline.
        """
        return ("%s\t%d\t%d\t%s\t%d\t%s\n" % (self.chr,self.start,self.end,misc.seq2nuID(self.readsequence),self.readcount,self.strand))
