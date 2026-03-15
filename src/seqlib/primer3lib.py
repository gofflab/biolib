"""Primer3 output parsing and primer design helpers.

Provides data classes (Record, Primer) for representing primer3 output and
a generator function for parsing primer3 Boulder-IO output files.  Also
includes a convenience wrapper for running primer3_core directly from Python.

Requires primer3 >= v2.2.
"""
import subprocess
import sys

from RNASeq import sequencelib


class Record(object):
    """Represents the primer3 output for a single target sequence.

    Attributes:
        sequenceID: Value of the SEQUENCE_ID field from the primer3 record.
        sequence: Value of the SEQUENCE_TEMPLATE field.
        comments: Comment line(s) associated with the record.
        primers: List of Primer objects describing primer pairs designed for
            this target sequence.
        attributes: Dictionary of other global parameters in the primer3
            record that are not specific to an individual primer pair.
    """
    def __init__(self):
        """Initialises a Record with empty/default attribute values."""
        self.sequenceID = ""
        self.sequence = ""
        self.comments = ""
        self.primers = []
        self.attributes = {}

    def __iter__(self):
        """Iterates over the Primer objects in this record."""
        return iter(self.primers)

    def __repr__(self):
        """Returns a short string representation of the record."""
        return "%s: %d primer pair(s)" % (self.sequenceID,len(self.primers))

class Primer(object):
    """Represents a single primer pair designed by Primer3.

    Attributes:
        sequenceID: ID of the target sequence for which this primer was
            designed (matches the parent Record's sequenceID).
        number: 1-based rank of this primer pair within the record.
        size: Deprecated field; use product_size instead.
        forward_seq: Sequence of the forward (left) primer.
        forward_start: 0-based start position of the forward primer on the
            template.
        forward_length: Length of the forward primer in bases.
        forward_tm: Melting temperature of the forward primer in °C.
        forward_gc: GC content of the forward primer as a percentage.
        reverse_seq: Sequence of the reverse (right) primer.
        reverse_start: 0-based start position of the reverse primer on the
            template.
        reverse_length: Length of the reverse primer in bases.
        reverse_tm: Melting temperature of the reverse primer in °C.
        reverse_gc: GC content of the reverse primer as a percentage.
        product_size: Expected PCR product size in base pairs.
    """
    def __init__(self):
        """Initialises a Primer with zero/empty attribute values."""
        self.sequenceID=""
        self.number = 0
        self.size = 0
        self.forward_seq = ''
        self.forward_start = ''
        self.forward_length = ''
        self.forward_tm = 0.0
        self.forward_gc = 0.0
        self.reverse_seq = ''
        self.reverse_start = 0
        self.reverse_length = 0
        self.reverse_tm = 0.0
        self.reverse_gc = 0.0
        self.product_size = 0

    def __repr__(self):
        """Returns a short string representation showing the sequence ID, number, and primer sequences."""
        return "%s_%d\n\tFwd: %s\tRev: %s" % (self.sequenceID,self.number,self.forward_seq, self.reverse_seq)

def parse(handle):
    """Parses a primer3 Boulder-IO output file and yields Record objects.

    Reads lines from the file handle, accumulates them until a '=' record
    separator is encountered, then constructs a Record with its associated
    Primer objects and yields it.

    Args:
        handle: A readable file-like object containing primer3 output in
            Boulder-IO format (each record terminated by a line containing
            only '=').

    Yields:
        Record objects, one per primer3 sequence entry.  Each Record contains
        a list of Primer objects corresponding to the primer pairs returned
        by primer3 for that sequence.

    Raises:
        StopIteration: When the end of the file is reached.
        KeyError: If a required primer3 output field is missing from a record.
    """
    recordLines = []
    while True:
        line = handle.readline().rstrip()
        if not line: raise StopIteration
        if not line == "=":
            recordLines.append(line)
            continue
        else:
            recordLines = [x.split("=") for x in recordLines]
            recordDict = dict(zip([x[0] for x in recordLines],[x[1] for x in recordLines]))
            rdKeys = recordDict.keys()
            record = Record()
            record.sequenceID = recordDict['SEQUENCE_ID']
            record.sequence = recordDict['SEQUENCE_TEMPLATE']
            try:
                nPrimers = int(recordDict['PRIMER_PAIR_NUM_RETURNED'])
            except KeyError:
                nPrimers=0
            for i in xrange(nPrimers):
                primer = Primer()
                primer.sequenceID = record.sequenceID
                primer.number = i+1
                primer.size = int(recordDict['PRIMER_PAIR_%d_PRODUCT_SIZE' % i])
                primer.forward_seq = recordDict['PRIMER_LEFT_%d_SEQUENCE' % i]
                primer.forward_start = int(recordDict['PRIMER_LEFT_%d' % i].split(",")[0])
                primer.forward_length = int(recordDict['PRIMER_LEFT_%d' % i].split(",")[1])
                primer.forward_tm = float(recordDict['PRIMER_LEFT_%d_TM' % i])
                primer.forward_gc = float(recordDict['PRIMER_LEFT_%d_GC_PERCENT' % i])
                primer.reverse_seq = recordDict['PRIMER_RIGHT_%d_SEQUENCE' % i]
                primer.reverse_start = int(recordDict['PRIMER_RIGHT_%d' % i].split(",")[0])
                primer.reverse_length = int(recordDict['PRIMER_RIGHT_%d' % i].split(",")[1])
                primer.reverse_tm = float(recordDict['PRIMER_RIGHT_%d_TM' % i])
                primer.reverse_gc = float(recordDict['PRIMER_RIGHT_%d_GC_PERCENT' % i])
                primer.product_size = int(recordDict['PRIMER_PAIR_%d_PRODUCT_SIZE' % i])
                record.primers.append(primer)
            yield record
            recordLines = []

#######
#Context specific runs
#######
def runPrimer3(fastaFile,task="qpcr",p3CloneSetFile="/seq/compbio-hp/lgoff/lincRNAs/primer_design/P3_cloning_primer_settings.p3",p3PCRSetFile="/seq/compbio-hp/lgoff/lincRNAs/primer_design/P3_qPCR_primer_settings.p3"):
    """Runs primer3_core on a FASTA file to design qPCR or cloning primers.

    Converts the FASTA file to Boulder-IO format and launches a primer3_core
    subprocess with the appropriate settings file.  The output file path is
    returned; note that the subprocess is not waited on before returning.

    Args:
        fastaFile: Path to a FASTA file of sequences to design primers for.
        task: Either 'qpcr' (default) to design short amplicon primers, or
            'cloning' to design full-length amplification primers using a
            defined included region.
        p3CloneSetFile: Path to the primer3 settings file used for cloning
            primer design.
        p3PCRSetFile: Path to the primer3 settings file used for qPCR primer
            design.

    Returns:
        Path to the primer3 output file (baseName + '.p3out').
    """

    baseName = fastaFile.rstrip(".fa")
    iter = sequencelib.FastaIterator(open(fastaFile,'r'))
    tmpFname = baseName+".p3in"
    tmpHandle = open(tmpFname,'w')

    #Make Boulder-IO format...
    for i in iter:
        myString = "SEQUENCE_ID=%s\nSEQUENCE_TEMPLATE=%s\n" % (i['name'],i['sequence'])
        if task == "cloning":
            myString += "SEQUENCE_INCLUDED_REGION=1,%d\n" % len(i['sequence'])
        myString += "="
        print(myString, file=tmpHandle)
    tmpHandle.close()

    P3Command = "primer3_core -p3_settings_file=%s -output=%s.p3out %s"

    sys.stderr.write("Designing Primers...\n")
    if task == "qpcr":
        subprocess.Popen(P3Command % (p3PCRSetFile,baseName+"_qPCR",tmpFname),shell=True)
    elif task == "cloning":
        subprocess.Popen(P3Command % (p3CloneSetFile,baseName+"_cloning",tmpFname),shell=True)
    return baseName+".p3out"
