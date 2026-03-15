"""Utilities for processing lincRNA (long intergenic non-coding RNA) transcript models.

Processes BED-format lincRNA annotations to fetch spliced sequences, insert
records into a MySQL database, generate transcript model PNG plots, and export
sequences to FASTA format.  Requires a local MySQL instance at the Broad
Institute and the intervallib package.
"""
import os
import sys

import intervallib

#from seqtools import dbConn
import MySQLdb


def main(bedFile,lincLotID):
    """Processes a BED file of lincRNA models and inserts them into the database.

    For each transcript in the BED file, fetches its spliced sequence,
    creates a PNG transcript model plot, and bulk-inserts all records into the
    lgoff_nextgen MySQL database using mysqlimport.

    Args:
        bedFile: Path to a BED-format file of lincRNA transcript models.
        lincLotID: Integer identifier for the lincRNA lot/batch being
            processed; used as a foreign key in the database insert.
    """

    #Setup environment
    if not os.path.exists('transcriptModels'):
        os.mkdir('transcriptModels')
    
    host="mysql.broadinstitute.org"
    user="lgoff"
    password=""
    db="lgoff_nextgen"
    
    tmpFname = 'transcripts.tab'
    tmpHandle = open(tmpFname,'w')
    
    #Make Database connection
    #db = getDb()
    
    #Make generator
    iter = intervallib.parseBed(bedFile)
    
    #Main loop
    for i in iter:
        #Fetch Sequence
        i.fetchSplicedSequence()
    
        #Make master tab-delim for insert
        print("\t".join(['NULL',i.name,i.chr,str(i.start),str(i.end),i.strand,",".join([str(x) for x in i.exonLengths]),",".join([str(x) for x in i.exonOffsets]),i.splicedSequence,str(lincLotID)]), file=tmpHandle)
        #insertRecord(i,lincLotID,db=db)
        
        #Make plots
        drawModelPNG(i,outDir='transcriptModels',verbose=True)
        
    
    
    #Close tmp file
    tmpHandle.close()
    
    #Do large insert into database
    os.system("mysqlimport -h %s -u %s -p%s %s %s") % (host,user,password,db,tmpFname)
    
    
    return

def drawModelPNG(bedRecord,outDir=os.getcwd(),verbose=False):
    """Generates a PNG transcript model image for a single BED record.

    Delegates to the BED record's makePNG method and optionally prints
    progress information to stdout.

    Args:
        bedRecord: An intervallib BED interval object that exposes a
            makePNG(outDir) method and a name attribute.
        outDir: Directory path where the PNG file will be written
            (default: current working directory).
        verbose: If True, print status messages indicating which transcript
            model is being drawn (default: False).
    """
    if verbose:
        print("Making transcript model plot...")
    bedRecord.makePNG(outDir)
    if verbose:
        print("\t"+bedRecord.name)
    return

def insertRecord(lincRNA,lincLotID):
    """Inserts a single lincRNA transcript record into the database.

    Constructs and executes an INSERT SQL statement for the transcripts table.
    The function references a module-level db cursor variable which must be
    set before calling.  Note: this function is known to be non-functional;
    use the bulk mysqlimport approach in main() instead.

    Args:
        lincRNA: An intervallib interval object with attributes: name, chr,
            start, end, strand, exonLengths, exonOffsets, and splicedSequence.
        lincLotID: Integer lot identifier to associate with the transcript
            record in the database.
    """
    
    cursor = db.cursor()
    insert="INSERT INTO transcripts VALUES (NULL,'%s','%s','%d','%d','%s','%s','%s','%s','%d');" % (lincRNA.name,lincRNA.chr,lincRNA.start,lincRNA.end,lincRNA.strand,",".join([str(x) for x in lincRNA.exonLengths]),",".join([str(x) for x in lincRNA.exonOffsets]),lincRNA.splicedSequence,int(lincLotID))
    cursor.execute(insert)
    try:
        db.commit()
        print(insert)
    except:
        db.rollback()
    return

def getDb():
    """Opens and returns a connection to the Broad Institute MySQL database.

    Connects to the lgoff_nextgen database on mysql.broadinstitute.org with
    a hard-coded user and empty password.

    Returns:
        A MySQLdb connection object for the lgoff_nextgen database.
    """
    host="mysql.broadinstitute.org"
    user="lgoff"
    password=""
    db="lgoff_nextgen"
    broadDb=MySQLdb.connect(host=host,user=user,db=db,passwd=password)
    return broadDb

def bed2Fa(fname):
    """Takes a .bed file input and makes a .fa file to be used for creating a reference set of sequences"""
    outHandle = open(fname.rstrip(".bed")+".fa",'w')
    iter = intervallib.parseBed(fname)
    
    for i in iter:
        i.fetchSplicedSequence()
        print(i.toFasta(), file=outHandle)
        sys.stderr.write(i.name+"\n")
    return    

##########################
#Setup Main
##########################

if __name__=="__main__":
    bedFile = sys.argv[1]
    lincLotID = sys.argv[2]
    main(bedFile,lincLotID)
