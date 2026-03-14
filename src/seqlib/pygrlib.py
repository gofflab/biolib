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
    def __init__(self, seq_id, start, stop, orientation):
        (self.id, self.start, self.stop, self.orientation) = \
            (seq_id, start, stop, orientation)


###GFF Futzing around

class GFF3Row(object):
    def __init__(self, line):
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
