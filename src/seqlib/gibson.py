"""Tools for designing Gibson Assembly fragments from FASTA sequences.

Reads a FASTA file of sequences (e.g. cDNAs or genomic regions) and splits
each into overlapping fragments suitable for Gibson Assembly cloning.
Optionally prepends Gateway attB recombination sequences to the outermost
primers.  Fragments are written in a tab-delimited or pretty-printed format.

Usage::

    python gibson.py [options] <fastaFile.fa>
"""
#Imports
import getopt
import sys

from RNASeq import sequencelib

#Fixed attributes
attF = "GGGGACAAGTTTGTACAAAAAAGCAGGCT" #Sequence to be added to the forward primer for Gateway (TM) cloning
attR = "GGGGACCACTTTGTACAAGAAAGCTGGGT" #Sequence to be added to the reverse primer for Gateway (TM) cloning

#Error trapping
help_message = '''
usage:
python gibson.py [options] <fastaFile.fa>

options:
    -h or --help      Prints this helpful help message
    -o or --output    output file for pretty results (default = <fastaFile_primers.txt>
    -g                Add attB sites for gateway cloning
    -f                Fragment size (default: 500bp)
    -v                Verbose output
    -s                overhang size (default: 20bp)
    -t                tab-delimited output (more machine readable)
'''

class Usage(Exception):
    """Exception raised for command-line usage errors.

    Attributes:
        msg: Human-readable explanation of the error or the help message.
    """
    def __init__(self, msg):
        """Initialises a Usage exception with an error message.

        Args:
            msg: Human-readable error or help text.
        """
        self.msg = msg

def gibson(fname,gateway=True,fragSize=500,overhangSize=20):
    """Splits FASTA sequences into overlapping Gibson Assembly fragments.

    Reads each record from a FASTA file and divides its sequence into a series
    of fragments of approximately fragSize bp, with consecutive fragments
    overlapping by overhangSize bp.  When gateway is True, the Gateway attB
    forward site (attF) is prepended to the sequence and the reverse
    complement of the Gateway attB reverse site (attR) is appended before
    fragmentation.

    Args:
        fname: Path to a FASTA-format input file.
        gateway: If True, add Gateway attB recombination sequences flanking
            the insert before fragmentation (default: True).
        fragSize: Target size in base pairs for each Gibson fragment
            (default: 500).
        overhangSize: Length in base pairs of the overlap between adjacent
            fragments (default: 20).

    Returns:
        A dictionary mapping each FASTA record name to a list of fragment
        sequence strings in 5'-to-3' order.
    """
    res = {}

    #Fasta file handle
    handle = open(fname,'r')
    iter = sequencelib.FastaIterator(handle)

    #Iterate over records in input fasta file
    for i in iter:
        fragments = []
        seq = i['sequence'].upper()
        if gateway:
            seq = attF + seq + sequencelib.rcomp(attR)
        curpos = 0
        length = int(len(seq)-1)
        while curpos < length:
            if curpos < 0:
                curpos = 0
            fragStart = curpos
            fragEnd = min(curpos+fragSize,length)
            #print "%d\t%d" % (fragStart,fragEnd)
            fragSeq = seq[int(fragStart):int(fragEnd)]
            fragments.append(fragSeq)
            curpos = curpos+fragSize-overhangSize
        res[i['name']]=fragments

    return res

def printGibson(fragDict,outHandle):
    """Writes Gibson Assembly fragments to a file handle in tab-delimited format.

    For each sequence in fragDict, prints a header line with the sequence name
    followed by one line per fragment in the format:
        <name>_block<N>\\t<fragment_sequence>

    Args:
        fragDict: Dictionary mapping sequence names to lists of fragment
            sequence strings, as returned by gibson().
        outHandle: Writable file-like object to receive the output.
    """
    for k in fragDict.keys():
        print("%s:" % k, file=outHandle)
        blockCount = 0
        for fragment in fragDict[k]:
            blockCount += 1
            print("%s_block%d\t%s" % (k,blockCount,fragment), file=outHandle)
        print("\n", file=outHandle)



##############
# Main
##############
def main(argv=None):
    """Command-line entry point for the Gibson Assembly fragment designer.

    Parses command-line arguments, calls gibson() to generate fragments from
    the provided FASTA file, and writes the results with printGibson().

    Args:
        argv: List of command-line argument strings.  Defaults to sys.argv
            when None.

    Raises:
        SystemExit: On usage errors or when --help is requested.
    """
    if argv is None:
        argv = sys.argv
    verbose = False
    outFile = None
    gateway = False
    keepTmp = False
    tabDelim = False
    overhangSize = 20
    fragSize = 500
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hto:vs:gf:k", ["help", "output="])
        except getopt.error as msg:
            raise Usage(msg)
        # option processing
        for option, value in opts:
            if option == "-v":
                verbose = True
            if option == "-g":
                gateway = True
            if option == "-f":
                fragSize == value
            if option == "-k":
                keepTmp=True
            if option in ("-h", "--help"):
                raise Usage(help_message)
            if option in ("-o", "--output"):
                outFile = value
            if option == "-s":
                overhangSize=value
            if option == "-t":
                tabDelim = True
        try:
            assert len(args)==1
            fname=args[0]
        except:
            raise Usage(help_message)
        if outFile == None:
            outFile = fname.rstrip(".fa")+"_gibson.txt"
        outHandle = open(outFile,'w')

        #Put actual function call here...
        fragDict = gibson(fname,gateway=gateway,fragSize=fragSize,overhangSize=overhangSize)
        #pp(fragDict)
        printGibson(fragDict,outHandle)

    except Usage as err:
        print(sys.argv[0].split("/")[-1] + ": " + str(err.msg), file=sys.stderr)
        print("\t for help use --help", file=sys.stderr)
        sys.exit()

if __name__ == "__main__":
    sys.exit(main())
