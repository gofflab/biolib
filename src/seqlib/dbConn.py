#!/usr/bin/env python
"""Database connection helpers and genomic data retrieval utilities.

Provides connection factories for several MySQL databases (Broad Institute
internal, UCSC Genome Browser public mirror, local UCSC mirror on 'valor', and
Ensembl) and a collection of query functions for fetching RefSeq transcripts,
wgRNA annotations, CpG islands, repeat overlaps, lincRNA records, and miRNA
seed sequences.

Most connection functions require network access to specific internal or public
servers and appropriate credentials.
"""
import sys
import time

import genomelib
import intervallib
import MySQLdb
import sequencelib


###################
#
#Connect to Broad MySQL Database
#
###################
def broadConnect():
    """Opens a DictCursor connection to the Broad Institute MySQL database.

    Connects to the lgoff_nextgen schema on mysql.broadinstitute.org.

    Returns:
        A MySQLdb DictCursor for the lgoff_nextgen database.
    """
    host="mysql.broadinstitute.org"
    user="lgoff"
    password=""
    db="lgoff_nextgen"
    broadDb=MySQLdb.connect(host=host,user=user,db=db,passwd=password)
    return broadDb.cursor(MySQLdb.cursors.DictCursor)
    
###################
#
#Connection to UCSC Genome Browser MySQL Database
#
###################
def gbdbConnect(gbdbname = "hg18"):
    """Opens a DictCursor connection to the UCSC Genome Browser public MySQL mirror.

    Args:
        gbdbname: UCSC genome database name (default: 'hg18').

    Returns:
        A MySQLdb DictCursor for the specified UCSC genome database.
    """
    gbHost = "genome-mysql.cse.ucsc.edu"
    gbUser = "genome"
    gbdb = MySQLdb.connect(host=gbHost,user=gbUser,db=gbdbname)
    return gbdb.cursor(MySQLdb.cursors.DictCursor)

###################
#
#Connection to Valor local UCSC Genome Browser MySQL Database
#
###################
def valorGbdbConnect(gbdbname='hg19'):
    """Opens a DictCursor connection to the local UCSC Genome Browser mirror on 'valor'.

    Connects to a locally hosted UCSC mirror database using the root account
    without a password.

    Args:
        gbdbname: Local UCSC genome database name (default: 'hg19').

    Returns:
        A MySQLdb DictCursor for the specified local genome database.
    """
    gbHost = 'localhost'
    gbUser = 'root'
    gbPass = ''
    gbdb = MySQLdb.connect(host=gbHost,user=gbUser,passwd=gbPass,db=gbdbname)
    return gbdb.cursor(MySQLdb.cursors.DictCursor)

###################
#
#Connection to Ensembl MySQL Database
#
####################
def ensemblConnect():
    """Opens a DictCursor connection to the public Ensembl MySQL server.

    Connects to the homo_sapiens_core_47_36i schema on ensembldb.ensembl.org
    using the anonymous account.

    Returns:
        A MySQLdb DictCursor for the Ensembl homo_sapiens_core_47_36i database.
    """
    ensemblHost = "ensembldb.ensembl.org"
    ensemblUser = "anonymous"
    ensembldbname = "homo_sapiens_core_47_36i"
    ensembldb = MySQLdb.connect(host=ensemblHost,user=ensemblUser,db=ensembldbname)
    return ensembldb.cursor(MySQLdb.cursors.DictCursor)

####################
#
#Operations on UCSC genome browser data
#
####################
def fetchRefSeq(genome = 'hg18',lookupval = 'name'):
    """Returns a dictionary of RefSeq genes (by chromosome and strand with 'name' parameter as key) from UCSC genome browser (equivalent to RefSeq ID)"""
    cursor=gbdbConnect(gbdbname=genome)
    select="SELECT * FROM refGene"
    cursor.execute(select)
    rows=cursor.fetchall()
    output={}
    for chr in genomelib.chr_names:
        output[chr]={}
        output[chr]['+']={}
        output[chr]['-']={}
    for row in rows:
        if row['chrom'] in genomelib.chr_names:
            output[row['chrom']][row['strand']][row[lookupval]]=row
    return output 

def fetchRefSeqIntervals(genome = 'hg18'):
    """Returns a dictionary of RefSeq SplicedInterval objects keyed by transcript name.

    Queries the refGene table of the UCSC Genome Browser database and
    constructs an intervallib.SplicedInterval for each transcript.

    Args:
        genome: UCSC genome database name (default: 'hg18').

    Returns:
        A dictionary mapping RefSeq transcript names to SplicedInterval
        objects.
    """
    cursor = gbdbConnect(gbdbname=genome)
    select = "SELECT * from refGene"
    cursor.execute(select)
    rows = cursor.fetchall()
    output = {}
    for row in rows:
        exonStarts = map(int,row['exonStarts'].rstrip().split(","))
        exonEnds = map(int,row['exonEnds'].rstrip().split(","))
        start = int(row['txStart'])
        exonOffsets = [x-start for x in exonStarts]
        exonLengths = []
        for i in len(exonStarts):
            exonLengths.append(exonEnds-exonStarts+1)
        output[row['name']] = intervallib.SplicedInterval(row['chrom'],row['txStart'],row['txEnd'],row['strand'],",".join([str(x) for x in exonLengths]),",".join([str(x) for x in exonOffsets]),name=row['name2'])
    return output

def fetchRefSeqIntervalsIndexed(genome='hg18',proteinCodingOnly=False,verbose=False):
    """
    Returns a dictionary of RefSeq SplicedIntervals (by chromosome and strand) from UCSC table browser.
    Indexed lists are sorted prior to return for easy search
    Same as fetchRefSeqIntervals but indexed by chrom and strand
    """
    cursor=gbdbConnect(gbdbname=genome)
    select="SELECT * FROM refGene"
    if verbose:
        sys.stderr.write("Fetching RefSeq Sequences...\n")
    cursor.execute(select)
    rows=cursor.fetchall()
    output={}
    for chr in genomelib.chr_names:
        output[chr]={}
        output[chr]['+']=[]
        output[chr]['-']=[]
    if verbose:
        sys.stderr.write("Creating index by chr and strand...\n")
    
    for row in rows:
        if proteinCodingOnly and not row['name'].startswith('NM'):
            continue
        try:
            exonStarts = map(int,row['exonStarts'].rstrip().split(",")[:-1])
            exonEnds = map(int,row['exonEnds'].rstrip().split(",")[:-1])
        except:
            print("\t".join(["%s:%s" % (k,v) for k,v in row.iteritems()]))
        start = int(row['txStart'])
        exonOffsets = [x-start for x in exonStarts]
        exonLengths = []
        for i in xrange(len(exonStarts)):
            exonLengths.append(exonEnds[i]-exonStarts[i]+1)
        if row['chrom'] in genomelib.chr_names:
            output[row['chrom']][row['strand']].append(intervallib.SplicedInterval(row['chrom'],row['txStart'],row['txEnd'],row['strand'],",".join([str(x) for x in exonLengths]),",".join([str(x) for x in exonOffsets]),name=row['name2']))
    
    #Sort 
    if verbose:
        sys.stderr.write("Sorting:\n")
    tstart = time.time()
    for key in output.keys():
        if verbose:
            sys.stderr.write("\t%s\t" % key)
        output[key]['+'].sort()
        output[key]['-'].sort()
        tend = time.time()
        if verbose:
            sys.stderr.write('%0.2f sec\n' % (tend-tstart))
        tstart = time.time()
    return output

def getIntervalFromRefSeq(lookupval,genome='hg18',lookupkey= 'name2',verbose=False):
    """Returns SplicedInterval objects for RefSeq transcripts matching a lookup value.

    Queries the UCSC refGene table for rows where lookupkey equals lookupval
    and constructs an intervallib.SplicedInterval for each matching transcript.

    Args:
        lookupval: The value to search for (e.g. a gene symbol or transcript
            ID).
        genome: UCSC genome database name (default: 'hg18').
        lookupkey: refGene column to search against (default: 'name2', which
            corresponds to the gene symbol).
        verbose: If True, print the SQL query and row count to stderr
            (default: False).

    Returns:
        A list of SplicedInterval objects for the matching transcripts.
    """
    cursor = gbdbConnect(gbdbname=genome)
    select = """SELECT * FROM refGene WHERE %s = '%s'""" % (lookupkey,lookupval)
    if verbose:
        sys.stderr.write("Query: "+select+"\nFetching RefSeq Record(s)\n")
    cursor.execute(select)
    rows=cursor.fetchall()
    if verbose:
        sys.stderr.write("%d Rows returned...\n" % len(rows))
    output = []
    for row in rows:
        try: 
            exonStarts = map(int,row['exonStarts'].rstrip().split(",")[:-1])
            exonEnds = map(int,row['exonEnds'].rstrip().split(",")[:-1])
        except:
            print("\t".join(["%s:%s" % (k,v) for k,v in row.iteritems()]))
        start = int(row['txStart'])
        exonOffsets = [x-start for x in exonStarts]
        exonLengths = []
        for i in xrange(len(exonStarts)):
            exonLengths.append(exonEnds[i]-exonStarts[i]+1)
        output.append(intervallib.SplicedInterval(row['chrom'],row['txStart'],row['txEnd'],row['strand'],",".join([str(x) for x in exonLengths]),",".join([str(x) for x in exonOffsets]),name=row['name2']))
    return output

def getIntervalFromAll_mRNA(lookupval,genome='hg18',lookupkey='qName',verbose=False):
    """Returns SplicedInterval objects from the UCSC all_mrna alignment table.

    Queries the all_mrna table for mRNA BLAT alignments matching lookupval
    in the specified column, and constructs a SplicedInterval for each row.

    Args:
        lookupval: The value to search for (e.g. a GenBank accession).
        genome: UCSC genome database name (default: 'hg18').
        lookupkey: all_mrna column to search (default: 'qName', the query
            sequence name).
        verbose: If True, print the SQL query and row count to stderr
            (default: False).

    Returns:
        A list of SplicedInterval objects for the matching alignments.
    """
    cursor = gbdbConnect(gbdbname=genome)
    select = """SELECT * FROM all_mrna WHERE %s = '%s'""" % (lookupkey,lookupval)
    if verbose:
        sys.stderr.write("Query: "+select+"\nFetching all_mrna Record(s)\n")
    cursor.execute(select)
    rows=cursor.fetchall()
    if verbose:
        sys.stderr.write("%d Rows returned...\n" % len(rows))
    output = []
    for row in rows:
        try:
            exonStarts = map(int,row['tStarts'].rstrip().split(",")[:-1])
            blockSizes = map(int,row['blockSizes'].rstrip().split(",")[:-1])
            exonEnds = [exonStarts[i]+blockSizes[i] for i in xrange(len(exonStarts))]
        except:
            print("\t".join(["%s:%s" % (k,v) for k,v in row.iteritems()]))
        start = int(row['tStart'])
        exonOffsets = [x-start for x in exonStarts]
        exonLengths = [exonEnds[i]-exonStarts[i]+1 for i in xrange(len(exonStarts))]
        output.append(intervallib.SplicedInterval(row['tName'],start,int(row['tEnd']),row['strand'],",".join([str(x) for x in exonLengths]),",".join([str(x) for x in exonOffsets]),name=row['qName']))
    return output

def refseqTSS():
    """Uses fetchRefSeq to retrieve current RefSeq Sequences and then returns a sorted list of tuples (as value of chr.strand dictionaries) containing ('refSeqID','chr','tss','orientation')"""
    refSeqs=fetchRefSeq()
    output={}
    for chr in genomelib.chr_names:
        output[chr]=[]
        for strand in ['+','-']:
            for k in refSeqs[chr][strand]:
                v=refSeqs[chr][strand][k]
                if v['strand'] == "+":
                    tss=v['txStart']
                elif v['strand'] == "-":
                    tss=v['txEnd']
                tssInfo=(v['name'],v['chrom'],int(tss),v['strand'])
                output[chr].append(tssInfo)
            output[chr].sort(lambda x,y:cmp(x[2],y[2]))
    return output

def fetchwgRNA():
    """Returns all wgRNA entries from the UCSC Genome Browser indexed by chromosome, strand, and name.

    Queries the wgRna table of the default genome (hg18) and organises
    results into a nested dictionary structure.

    Returns:
        A dictionary of the form output[chr][strand][name] = row_dict for
        each wgRNA entry on a standard chromosome.
    """
    cursor=gbdbConnect()
    select="SELECT * FROM wgRna"
    cursor.execute(select)
    rows=cursor.fetchall()
    output={}
    for chr in genomelib.chr_names:
        output[chr]={}
        output[chr]['+']={}
        output[chr]['-']={}
    for row in rows:
        if row['chrom'] in genomelib.chr_names:
            output[row['chrom']][row['strand']][row['name']]=row
    return output


#Tests for known annotation
def hostRefSeq(chr,start,end,strand):
    """
    Checks to see if interval is within a host RefSeq gene (does not test strand!!).  If no, returns False.  
    If yes, returns a list of dictionaries for each host RefSeq gene.  Keys are consistent with field names 
    from UCSC table refGene.
    """
    cursor=gbdbConnect()
    selSQL="SELECT * from refGene WHERE chrom='%s' AND txStart<='%d' AND txEnd>='%d'" % (chr,int(start),int(end))
    cursor.execute(selSQL)
    rows=cursor.fetchall()
    results=[]
    if cursor.rowcount==0:
        return False
    else:
        for row in rows:
            results.append(row)
        return results

def testCpG(chr,start,end):
    """Tests whether a genomic interval overlaps a CpG island in the UCSC database.

    Queries the cpgIslandExt table for CpG islands that overlap the given
    coordinates.

    Args:
        chr: Chromosome name (e.g. 'chr1').
        start: Start coordinate (0-based).
        end: End coordinate.

    Returns:
        The first matching row as a dictionary, or False if no CpG island
        overlaps the interval.
    """
    cursor=gbdbConnect()
    selSQL="SELECT * from cpgIslandExt WHERE chrom='%s' AND chromStart<='%d' AND chromEnd>='%d'" % (chr,int(start),int(end))
    cursor.execute(selSQL)
    if cursor.rowcount==0:
        return False
    else:
        return cursor.fetchone()

def testwgRNA(chr,start,end,strand):
    """
    Checks to see if interval is entirely within a known wgRNA gene (including miRNA). Does consider strand!!!
    If no flanking host wgRNA, returns False. If yes, returns a list of dictionaries for each host wgRNA gene.
    Keys are consistent with field names from UCSC table wgRNA.
    """
    cursor=gbdbConnect()
    selSQL="SELECT * from wgRna WHERE chrom='%s' AND strand='%s' AND chromStart<='%d' AND chromEnd>='%d'" % (chr,strand,int(start),int(end))
    cursor.execute(selSQL)
    rows=cursor.fetchall()
    results=[]
    if cursor.rowcount==0:
        return False
    else:
        for row in rows:
            results.append(row)
        return results

def hostmRNA(chr,start,end,strand):
    """Returns mRNA alignments that span a given genomic interval from the UCSC database.

    Queries a chromosome-specific mRNA table (named <chr>_mrna) for
    alignments that contain the interval [start, end].

    Args:
        chr: Chromosome name (e.g. 'chr1').
        start: Start coordinate of the query interval.
        end: End coordinate of the query interval.
        strand: Strand orientation (not currently used in the SQL query).

    Returns:
        A list of row dictionaries for overlapping mRNA alignments, or False
        if none are found.
    """
    cursor=gbdbConnect()
    selSQL="SELECT * from %s_mrna WHERE tName='%s' AND tStart<='%d' AND tEnd>='%d'" % (chr,chr,int(start),int(end))
    cursor.execute(selSQL)
    rows=cursor.fetchall()
    results=[]
    if cursor.rowcount==0:
        return False
    else:
        for row in rows:
            results.append(row)
        return results

def fetchLincRNA(fname="/seq/compbio/lgoff/lincRNAs/hg18_lincRNA_Guttman.bed"):
    """Reads a lincRNA BED file and returns intervals indexed by chromosome.

    Parses a three-column BED file (chr, start, end) and organises the
    resulting intervals into a dictionary keyed by chromosome name.

    Args:
        fname: Path to a BED file of lincRNA intervals (default: hg18
            Guttman et al. lincRNA catalogue).

    Returns:
        A dictionary mapping chromosome names to lists of interval
        dictionaries, each with keys 'chr', 'start' (int), and 'end' (int).
    """
    handle=open(fname,'r')
    lincs={}
    for chr in genomelib.chr_names:
        lincs[chr]=[]
    for line in handle:
        if line.startswith("#"):continue
        fields=['chr','start','end']
        vals=line.rstrip().split("\t")
        d=dict(zip(fields,vals))
        d['start'],d['end']=int(d['start']),int(d['end'])
        lincs[d['chr']].append(d)
    return lincs

def fetchmiRNASeeds(fname="/seq/compbio/lgoff/smallRNAs/genomes/human/microRNA/mature.fa",species = 'hsa'):
    """Reads a miRBase FASTA file and returns a dictionary mapping seed sequences to miRNA names.

    Extracts the 7-nt seed sequence (positions 2-8 of the mature miRNA) for
    each entry matching the given species prefix.

    Args:
        fname: Path to a miRBase mature miRNA FASTA file (default: internal
            Broad Institute path).
        species: Two- or three-letter miRBase species prefix to filter by
            (default: 'hsa' for Homo sapiens).

    Returns:
        A dictionary mapping 7-nt seed sequences (str) to the first token
        of the miRNA name (str).
    """
    handle = open(fname,'r')
    seeds = {}
    iter = sequencelib.FastaIterator(handle)
    for i in iter:
        if i.name.startswith(species):
            seeds[i.sequence[1:8]] = i.name.split()[0]
    return seeds

#############
#Added for lincRNA pipeline (only works on valor)
############

def findRepeatOverlap(interval,cursor=None):
    """Returns RepeatMasker annotations that overlap a given genomic interval.

    Queries the rmsk table of the local UCSC mirror for repeat elements that
    partially or fully overlap the interval.

    Args:
        interval: An intervallib interval object with chr, start, end, and
            genome attributes.
        cursor: An optional pre-existing MySQLdb DictCursor.  If None, a new
            connection to the local valor UCSC mirror is opened.

    Returns:
        A list of row dictionaries for overlapping repeats, or False if none
        are found.
    """
    if cursor == None:
        cursor = valorGbdbConnect(interval.genome)
    selSQL = "SELECT * from rmsk WHERE genoName = '%s' AND (genoStart >= '%d' OR genoEnd >= '%d') AND (genoStart <= '%d' OR genoEnd <= '%d')" % (interval.chr,interval.start,interval.start,interval.end,interval.end)
    cursor.execute(selSQL)
    rows = cursor.fetchall()
    results=[]
    if cursor.rowcount==0:
        return False
    else:
        for row in rows:
            results.append(row)
        return results
    
def findUCSCOverlap(interval,cursor=None):
    """Returns UCSC knownGene entries (with RefSeq mapping) that overlap a given interval.

    Queries the knownGene table joined to knownToRefSeq on the local UCSC
    mirror for known genes that partially or fully overlap the interval.

    Args:
        interval: An intervallib interval object with chr, start, end, and
            genome attributes.
        cursor: An optional pre-existing MySQLdb DictCursor.  If None, a new
            connection to the local valor UCSC mirror is opened.

    Returns:
        A list of row dictionaries for overlapping known genes, or False if
        none are found.
    """
    if cursor == None:
        cursor = valorGbdbConnect(interval.genome)
    selSQL = "SELECT * from knownGene kg LEFT JOIN knownToRefSeq krs ON kg.name = krs.name WHERE kg.chrom = '%s' AND (kg.txStart >= '%d' OR kg.txEnd >= '%d') AND (kg.txStart <= '%d' OR kg.txEnd <= '%d')" % (interval.chr,interval.start,interval.start,interval.end,interval.end)
    cursor.execute(selSQL)
    rows = cursor.fetchall()
    results = []
    if cursor.rowcount == 0:
        return False
    else:
        for row in rows:
            results.append(row)
        return results
