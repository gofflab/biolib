"""Parsing and data structures for GTF (Gene Transfer Format) files.

All of this is very fragile and is absolutely dependent on a unique geneId
and unique transcriptId for any records.

Provides GTF_Entry, GTFTranscriptContainer, and GTFGeneContainer classes for
holding GTF data, along with iterator functions for streaming over transcripts
and genes, and utility functions for building attribute dictionaries and tables.

Originally created on Aug 31, 2010.

Author: lgoff
"""
###########
#Imports
###########
import sys

from . import intervallib
from .misc import uniqify

#import genomelib

#######################
#Error Handling
#######################
class Error(Exception):
    """Base class for exceptions in this module.

    Provides a message property with getter/setter so subclasses can store
    a human-readable error description.
    """
    def __str__(self):
        """Return the string representation of the error message."""
        return str(self.message)
    def _get_message(self, message): return self._message
    def _set_message(self, message): self._message = message
    message = property(_get_message, _set_message)

class ParsingError(Error):
    """
    Exception raised for errors in the input.

    Attributes:
        message -- explanation of the error
    """
    def __init__(self, message):
        self.message = message

#########################
#GTF Entry Class
#########################

class GTF_Entry:
    """Holds a single row's worth of GTF/GFF information.

    Attributes:
        contig: Sequence name / chromosome.
        source: Annotation source name.
        feature: Feature type (e.g. "exon", "CDS", "transcript").
        frame: Reading frame (".","0","1","2").
        start: 1-based start coordinate (integer).
        end: 1-based end coordinate (integer).
        score: Score value (float or ".").
        strand: Strand ("+" or "-" or ".").
        attributes: Dictionary of parsed key-value attribute pairs.
    """

    def __init__(self):
        """Construct a GTF_Entry with default empty/sentinel field values."""
        self.contig = "."
        self.source = "."
        self.feature = "."
        self.frame = "."
        self.start = 0
        self.end = 0
        self.score = "."
        self.strand = "."
        self.attributes = {}

    def __lt__(self, b):
        """Compare GTF entries by midpoint coordinate."""
        return (self.start + self.end) // 2 < (b.start + b.end) // 2

    def __eq__(self, b):
        """Return True if two GTF entries share the same midpoint coordinate."""
        return (self.start + self.end) // 2 == (b.start + b.end) // 2

    def __repr__(self):
        """Return a transcript_id:feature string representation."""
        return self.attributes['transcript_id']+":"+self.feature

    def addGTF_Entry(self,gtf_entry):
        """Copy all fields from another GTF_Entry into self.

        Args:
            gtf_entry: A GTF_Entry instance whose fields will be copied.
        """
        self.contig = gtf_entry.contig
        self.source = gtf_entry.source
        self.feature = gtf_entry.feature
        self.frame = gtf_entry.frame
        self.start = int(gtf_entry.start)
        self.end = int(gtf_entry.end)
        self.score = gtf_entry.score
        self.strand = gtf_entry.strand
        self.attributes = gtf_entry.attributes

    def read(self,line):
        """
        read gff entry from line.
        <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
        """
        data = line.rstrip().split("\t")

        try:
            (self.contig, self.source, self.feature,
             self.start, self.end, self.score, self.strand,
             self.frame ) = [x.strip() for x in data[:8]]
        except ValueError:
            raise ValueError( "parsing error in line `%s`" % line )

        ## note: frame might be "."
        (self.start, self.end) = map(int, (self.start, self.end))
        try:
            self.score = float(self.score)
        except:
            pass
        #TODO: This may only be necessary when I convert to an Interval object
        #self.start -= 1

        self.parseInfo( data[8], line )

    def parseInfo(self,myAttributes,line ):
        """
        Parse attributes.
        """
        # remove comments
        myAttributes = myAttributes.split( "#" )[0]
        # separate into fields
        fields = [x.strip() for x in myAttributes.split(";")[:-1]]
        self.attributes = {}

        for f in fields:
            d = [x.strip() for x in f.split(" ")]
            n,v = d[0], d[1]
            if len(d) > 2: v = d[1:]
            if v[0] == '"' and v[-1] == '"':
                v = v[1:-1]
            else:
                ## try to convert to a value
                try:
                    v = float( v )
                    v = int( v )
                except ValueError:
                    pass
                except TypeError:
                    pass
            self.attributes[n] = v

    def toGTF(self):
        """Serialize this entry back to a GTF-formatted string.

        Writes gene_id and transcript_id first (as required by the GTF spec),
        then all remaining attributes in arbitrary order.

        Returns:
            A GTF-formatted line string ending with a newline.
        """
        tmp = '%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t' % (self.contig,self.source,self.feature,self.start,self.end,str(self.score),self.strand,self.frame)
        #Print 'gene_id' and 'transcript_id' as first and second attributes (required by GTF spec.)
        for attr in ['gene_id','transcript_id']:
            try:
                tmp += '%s "%s"; ' % (attr,self.attributes[attr])
            except:
                pass
        #Print remainder of attributes in any order.
        for k,v in self.attributes.items():
            if k not in ['gene_id','transcript_id']:
                tmp += '%s "%s"; ' % (k,str(v))
        tmp += "\n"
        return tmp

############
#GTFTranscriptContainer
############
class GTFTranscriptContainer(object):
    """Container grouping all GTF_Entry instances sharing a transcript_id.

    Attributes:
        features: List of GTF_Entry objects belonging to this transcript.
        start: Minimum start coordinate across all features.
        end: Maximum end coordinate across all features.
        contig: Chromosome/contig name.
        strand: Strand orientation.
        transcriptId: transcript_id attribute value.
        geneId: gene_id attribute value.
    """

    def __init__(self):
        """Construct an empty GTFTranscriptContainer with sentinel values."""
        self.features = []
        self.start = -1
        self.end = -1
        self.contig = None
        self.strand = None
        self.transcriptId = ''
        self.geneId = ''

    def __len__(self):
        """Return the genomic span of the transcript (end - start + 1)."""
        return self.end-self.start+1

    def __lt__(self, b):
        """Compare transcript containers by midpoint coordinate."""
        return (self.start + self.end) // 2 < (b.start + b.end) // 2

    def __eq__(self, b):
        """Return True if two transcript containers share the same midpoint."""
        return (self.start + self.end) // 2 == (b.start + b.end) // 2

    def addFeature(self,gtf_entry):
        """Add a GTF_Entry to this transcript container.

        Initialises contig, strand, and transcriptId from the first feature
        added. Asserts that subsequent features share the same transcript_id.
        Updates self.start and self.end to span all features.

        Args:
            gtf_entry: A GTF_Entry instance to add.

        Raises:
            AssertionError: If gtf_entry has a different transcript_id.
        """
        if self.transcriptId == '':
            self.contig = gtf_entry.contig
            self.strand = gtf_entry.strand
            self.transcriptId = gtf_entry.attributes['transcript_id']
        assert self.transcriptId == gtf_entry.attributes['transcript_id']
        if 'gene_id' in gtf_entry.attributes.keys():
            self.geneId = gtf_entry.attributes['gene_id']
        self.features.append(gtf_entry)
        self.update()

    def update(self):
        """Recompute self.start and self.end from the current feature list."""
        self.start = min([x.start for x in self.features])
        self.end = max([x.end for x in self.features])

    def toSplicedInterval(self):
        """Convert this transcript container to a SplicedInterval.

        Extracts exon features, sorts them by exon_number, and constructs a
        SplicedInterval using their lengths and offsets.

        Returns:
            A SplicedInterval representing the spliced transcript.

        Raises:
            ValueError: If more than one distinct transcript_id is found
                in the feature list.
        """
        transcripts = uniqify([x.attributes['transcript_id'] for x in self.features])
        if len(transcripts) > 1:
            raise ValueError ("Something is wrong, there are too many different transcript_ids")
        for t in transcripts:
            exons = [x for x in self.features if (x.feature =='exon' and x.attributes['transcript_id']==t)]
            exons.sort(key=lambda x: x.attributes['exon_number'])
            transStart = min([x.start-1 for x in exons])
            myInt = intervallib.SplicedInterval(self.contig,transStart,max([x.end for x in exons]),self.strand,",".join([str(x.end-x.start+1) for x in exons]),",".join([str(x.start-1-transStart) for x in exons]),name=t)
            return myInt


############
#Gene Container
############

class GTFGeneContainer(object):
    """Container for all GTF_Entry instances sharing a common gene_id.

    Assumptions:
        - The gene_id field is unique to a gene locus (not shared among
          gene duplicates).
        - There is no guarantee that the row order is preserved during
          reading or when returning GTF output.

    Attributes:
        features: List of GTF_Entry objects for this gene.
        transcripts: List of GTFTranscriptContainer objects for this gene.
        start: Minimum start coordinate across all features/transcripts.
        end: Maximum end coordinate across all features/transcripts.
        contig: Chromosome/contig name.
        strand: Strand orientation.
        geneId: gene_id attribute value.
        sequence: DNA sequence string (empty by default).
    """

    def __init__(self):
        """Construct an empty GTFGeneContainer with sentinel values."""
        self.features = []
        self.transcripts = []
        self.start = -1
        self.end = -1
        self.contig = None
        self.strand = None
        self.geneId = ''
        self.sequence = ''

    def __len__(self):
        """Return the genomic span of the gene (end - start + 1)."""
        return self.end-self.start+1

    def __lt__(self, b):
        """Compare gene containers by midpoint coordinate."""
        return (self.start + self.end) // 2 < (b.start + b.end) // 2

    def __eq__(self, b):
        """Return True if two gene containers share the same midpoint."""
        return (self.start + self.end) // 2 == (b.start + b.end) // 2

    def addFeature(self,gtf_entry):
        """Add a GTF_Entry feature to this gene container.

        Initialises contig, strand, and geneId from the first feature added.
        Asserts that subsequent features share the same gene_id. Updates
        self.start and self.end.

        Args:
            gtf_entry: A GTF_Entry instance to add.

        Raises:
            AssertionError: If gtf_entry has a different gene_id.
        """
        if self.geneId == '':
            self.contig = gtf_entry.contig
            self.strand = gtf_entry.strand
            self.geneId = gtf_entry.attributes['gene_id']
        assert self.geneId == gtf_entry.attributes['gene_id']
        self.features.append(gtf_entry)
        self.update()

    def addGTFTranscript(self,gtf_transcript):
        """Add a GTFTranscriptContainer to this gene container.

        Initialises contig, strand, and geneId from the first transcript added.
        Asserts that subsequent transcripts share the same geneId, contig, and
        strand. Updates self.start and self.end via transcriptUpdate().

        Args:
            gtf_transcript: A GTFTranscriptContainer instance to add.

        Raises:
            AssertionError: If geneId, contig, or strand do not match.
        """
        if self.geneId == '':
            self.contig = gtf_transcript.contig
            self.strand = gtf_transcript.strand
            self.geneId = gtf_transcript.geneId
        assert self.geneId == gtf_transcript.geneId and self.contig == gtf_transcript.contig and self.strand == gtf_transcript.strand
        self.transcripts.append(gtf_transcript)
        self.transcriptUpdate()

    def update(self):
        """Recompute self.start and self.end from the current features list."""
        self.start = min([x.start for x in self.features])
        self.end = max([x.end for x in self.features])

    def transcriptUpdate(self):
        """Recompute self.start and self.end from the transcripts list."""
        self.start = min([x.start for x in self.transcripts])
        self.end = max([x.end for x in self.transcripts])


    def propogateLincName(self,lincName):
        """Set the linc_name attribute on all features, and gene_name if absent.

        Args:
            lincName: The lincRNA name string to propagate to all features.
        """
        for feat in self.features:
            feat.attributes['linc_name'] = lincName
            if 'gene_name' not in feat.attributes:
                feat.attributes['gene_name'] = lincName

    def addAttribute(self,key,value):
        """Add or overwrite an attribute key-value pair on all features.

        Args:
            key: Attribute name string.
            value: Attribute value to assign.
        """
        for feat in self.features:
            feat.attributes[key] = value

    def geneToBed(self):
        """This does not work yet"""
        raise Error ("This method does not work yet")
        return "%s\t%d\t%d\t%s\t0\t%s\t%s\t%s" % (self.contig,self.start,self.end,self.attributes['transcript_id'],self.strand,",".join(self.exonLengths),",".join(self.exonOffsets))

    def transcriptsToBed(self):
        """Placeholder for BED output of transcripts (not yet implemented)."""
        pass

    def getGTF(self):
        """Return a GTF string containing all features of this gene.

        Returns:
            Multi-line string of GTF-formatted rows for every feature.
        """
        tmp = ''
        for feat in self.features:
            tmp += feat.toGTF()
        return tmp

    def toInterval(self):
        """Convert this gene to an Interval spanning its genomic footprint.

        Returns:
            An Interval with 0-based start (start-1), end, strand, and the
            gene_id as its name.
        """
        return intervallib.Interval(self.contig,self.start-1,self.end,self.strand,name=self.geneId)

    # def fetchSequence(self,genome='hg19',connection=None):
    #     if connection == None:
    #         connection = genomelib.pygrConnect(genome)
    #     try:
    #         seq = connection[self.contig][self.start-1:self.end]
    #     except KeyError:
    #         seq = ''
    #     self.sequence=str(seq)
    #     return


#############
#lineIterator
#############
def lineIterator(gtfHandle):
    """Yield GTF_Entry objects for every non-comment line in gtfHandle.

    Skips lines starting with "#". Parses each remaining line into a
    GTF_Entry via GTF_Entry.read().

    Args:
        gtfHandle: An open file handle to a GTF file.

    Yields:
        GTF_Entry objects, one per data line.
    """
    while True:
        line = gtfHandle.readline()
        if not line: return
        if line.startswith("#"):continue
        gtf_entry = GTF_Entry()
        gtf_entry.read(line)
        yield gtf_entry

def GTFGeneIterator(gtfFile,verbose = False):
    """Iterate over genes in a GTF file, yielding one GTFGeneContainer per gene.

    Groups all GTF_Entry rows by gene_id and yields a fully-populated
    GTFGeneContainer for each unique gene_id found.

    Args:
        gtfFile: Path to the GTF file.
        verbose: If True, write progress messages to stderr (default False).

    Yields:
        GTFGeneContainer objects, one per unique gene_id.
    """
    handle = open(gtfFile,'r')
    iter = lineIterator(handle)
    res = {}
    if verbose:
        sys.stderr.write("Parsing GTF lines into genes...\n")
    for i in iter:
        res.setdefault(i.attributes['gene_id'],GTFGeneContainer())
        res[i.attributes['gene_id']].addFeature(i)
    for k in res.keys():
        yield res[k]

def GTFGeneIterator2(gtfFile,verbose=False):
    """Iterate over genes by grouping transcripts, yielding one GTFGeneContainer per gene.

    An alternative to GTFGeneIterator that builds genes from
    GTFTranscriptContainer objects rather than raw GTF_Entry rows.

    Args:
        gtfFile: Path to the GTF file.
        verbose: If True, write progress messages to stderr (default False).

    Yields:
        GTFGeneContainer objects, one per unique gene_id.
    """
    iter = GTFTranscriptIterator(gtfFile,verbose=verbose)
    res = {}
    for i in iter:
        res.setdefault(i.geneId,GTFGeneContainer())
        res[i.geneId].addGTFTranscript(i)
    for k in res.keys():
        yield res[k]

def GTFTranscriptIterator(gtfFile,verbose = False):
    """Iterate over transcripts in a GTF file, yielding one GTFTranscriptContainer per transcript.

    Groups all GTF_Entry rows by transcript_id and yields a fully-populated
    GTFTranscriptContainer for each unique transcript_id found.

    Args:
        gtfFile: Path to the GTF file.
        verbose: If True, write progress messages to stderr (default False).

    Yields:
        GTFTranscriptContainer objects, one per unique transcript_id.
    """
    handle = open(gtfFile,'r')
    iter = lineIterator(handle)
    res = {}
    if verbose:
        sys.stderr.write("Parsing GTF lines into transcripts...\n")
    for i in iter:
        res.setdefault(i.attributes['transcript_id'],GTFTranscriptContainer())
        res[i.attributes['transcript_id']].addFeature(i)
    for k in res.keys():
        yield res[k]

def GTFAttributeDict(gtfFile,idField='gene_id'):
    """Returns a dictionary of attributes for each unique gene_id"""
    handle = open(gtfFile)
    res = {}
    fields = set([])
    for line in handle:
        if line.startswith("#"):continue
        attributes = line.rstrip().split("\t")[8].split(";")[:-1]
        attrs = [ x.strip().split(" ")[0] for x in attributes]
        fields.update(attrs)
        values = [ x.strip().split(" ")[1].strip('"') for x in attributes]
        myDict = dict(zip(attrs,values))
        res.setdefault(myDict[idField],{})
        for k,v in myDict.items():
            res[myDict[idField]].setdefault(k,set([])).add(v)
    return res

def GTFAttributeTable(gtfFile,outfile,idField='gene_id'):
    """writes a table of attributes for each unique gene_id"""
    handle = open(gtfFile)
    outHandle = open(outfile,'w')
    res = {}
    fields = set([])
    for line in handle:
        if line.startswith("#"):continue
        attributes = line.rstrip().split("\t")[8].split(";")[:-1]
        attrs = [ x.strip().split(" ")[0] for x in attributes]
        fields.update(attrs)
        values = [ x.strip().split(" ")[1].strip('"') for x in attributes]
        myDict = dict(zip(attrs,values))
        res.setdefault(myDict[idField],{})
        for k,v in myDict.items():
            res[myDict[idField]].setdefault(k,set([])).add(v)

    #Print output to outHandle
    #Header
    print("%s\t%s" % (idField,"\t".join([str(x) for x in fields])), file=outHandle)

    for key in res.keys():
        outString = '%s\t' % key
        for field in fields:
            try:
                outString += ",".join(res[key][field]) + "\t"
            except KeyError:
                outString += "-\t"
        outString.rstrip("\t")
        print(outString, file=outHandle)
    return

def test():
    """Placeholder test function. No-op.

    Example usage (Python 2 style, for reference)::

        from RNASeq import GTFlib
        fname = 'linc_catalog.gtf'
        iter = GTFlib.GTFGeneIterator(fname)
        for i in iter:
            print i.getGTF(),
    """
    pass



