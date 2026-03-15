'''
Created on Jun 23, 2011

@author: lgoff

NOTE: This module depends on 'pygr', an unmaintained Python 2-only library.
It is kept for reference only and is not functional in Python 3.
Do not import this module in production code.
'''

# NOTE: pygr is not available in Python 3. Imports are guarded below.
try:
    from pygr import annotation, mapping, worldbase
    _PYGR_AVAILABLE = True
except ImportError:
    _PYGR_AVAILABLE = False


###Classes
class MySliceInfo(object):
    """Stores coordinate information for a genomic slice in pygr convention.

    Holds the four fields required to identify a sequence slice: sequence ID,
    start, stop (exclusive), and orientation (+1 or -1).

    Attributes:
        id: Sequence (chromosome) identifier.
        start: 0-based start coordinate.
        stop: Exclusive end coordinate.
        orientation: Strand orientation; +1 for forward, -1 for reverse.
    """
    def __init__(self, seq_id, start, stop, orientation):
        """Initialises a MySliceInfo.

        Args:
            seq_id: Sequence (chromosome) identifier.
            start: 0-based start coordinate of the slice.
            stop: Exclusive end coordinate of the slice.
            orientation: Strand orientation; +1 for forward strand, -1 for
                reverse strand (pygr convention).
        """
        (self.id, self.start, self.stop, self.orientation) = \
            (seq_id, start, stop, orientation)


###GFF Futzing around

class GFF3Row(object):
    """Represents a single data row from a GFF3 annotation file.

    Parses one GFF3 line and stores the type, sequence ID, start/stop
    coordinates (converted to 0-based pygr convention), strand orientation,
    and all key=value attributes from column 9.

    Attributes:
        type: Feature type string from column 3 (e.g. 'gene', 'exon').
        id: Sequence (chromosome) ID from column 1.
        start: 0-based start coordinate (GFF3 1-based column 4 minus 1).
        stop: Exclusive end coordinate (GFF3 column 5).
        orientation: +1 for '+' strand, -1 for '-' strand.
        Additional attributes are set dynamically from column 9 key=value
        pairs; multi-value attributes (comma-separated) are stored as lists.
    """
    def __init__(self, line):
        """Parses a GFF3 line into a GFF3Row object.

        Args:
            line: A single tab-delimited GFF3 data line (not a comment).

        Raises:
            ValueError: If the strand character in column 7 is not '+' or '-'.
        """
        cols = line.split('\t')
        self.type = cols[2]
        self.id = cols[0]  # sequence ID
        self.start = int(cols[3]) - 1  # correct for 1-based coords
        self.stop = int(cols[4])
        if cols[6] == '+':  # convert to Pygr convention
            self.orientation = 1
        elif cols[6] == '-':
            self.orientation = -1
        else:
            raise ValueError('Bad strand: %s' % cols[6])
        for s in cols[8].split(';')[:-1]:  # parse attributes
            attr, val = s.strip().split(' ')
            if ',' in val:
                setattr(self, attr, val.split(','))
            else:
                setattr(self, attr, val)


def read_gff3(filename, genome):
    """Reads a GFF3 annotation file and builds pygr AnnotationDB objects.

    Parses a GFF3 file, groups features by type, and creates one pygr
    AnnotationDB per feature type, each associated with the provided genome
    sequence database.  Comment lines (starting with '#') are skipped.
    Features lacking a type or gene_id attribute are also skipped.

    Args:
        filename: Path to a GFF3-format annotation file.
        genome: A pygr sequence database object (e.g. a worldbase genome)
            used to associate annotation slices with genomic sequence.

    Returns:
        A dictionary mapping feature type strings to pygr AnnotationDB
        objects.

    Raises:
        ImportError: If the pygr library is not installed.
    """
    if not _PYGR_AVAILABLE:
        raise ImportError("pygr is required for read_gff3 but is not installed.")
    d = {}  # for different types of sliceDBs
    with open(filename) as ifile:
        for line in ifile:  # parse all the GFF3 lines
            if line.startswith('#'):  # ignore this line
                continue
            row = GFF3Row(line)
            try:
                d.setdefault(row.type, {})[row.gene_id] = row
            except AttributeError:
                pass  # no type or ID so ignore...
    annotations = {}
    for atype, sliceDB in d.items():  # create annotation DBs
        adb = annotation.AnnotationDB(sliceDB, genome)
        annotations[atype] = adb
    return annotations
