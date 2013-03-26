'''
Created on Jun 3, 2010

@author: lgoff
'''
import intervallib
import os,sys
#from seqtools import dbConn
import MySQLdb

def main(bedFile,lincLotID):
    
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
        print >>tmpHandle, "\t".join(['NULL',i.name,i.chr,str(i.start),str(i.end),i.strand,",".join([str(x) for x in i.exonLengths]),",".join([str(x) for x in i.exonOffsets]),i.splicedSequence,str(lincLotID)])
        #insertRecord(i,lincLotID,db=db)
        
        #Make plots
        drawModelPNG(i,outDir='transcriptModels',verbose=True)
        
    
    
    #Close tmp file
    tmpHandle.close()
    
    #Do large insert into database
    os.system("mysqlimport -h %s -u %s -p%s %s %s") % (host,user,password,db,tmpFname)
    
    
    return

def drawModelPNG(bedRecord,outDir=os.getcwd(),verbose=False):
    if verbose:
        print "Making transcript model plot..."
    bedRecord.makePNG(outDir)
    if verbose:
        print "\t"+bedRecord.name
    return

def insertRecord(lincRNA,lincLotID):
    """Does not work for some reason..."""
    
    cursor = db.cursor()
    insert="INSERT INTO transcripts VALUES (NULL,'%s','%s','%d','%d','%s','%s','%s','%s','%d');" % (lincRNA.name,lincRNA.chr,lincRNA.start,lincRNA.end,lincRNA.strand,",".join([str(x) for x in lincRNA.exonLengths]),",".join([str(x) for x in lincRNA.exonOffsets]),lincRNA.splicedSequence,int(lincLotID))
    cursor.execute(insert)
    try:
        db.commit()
        print insert
    except:
        db.rollback()
    return

def getDb():
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
        print >>outHandle, i.toFasta()
        sys.stderr.write(i.name+"\n")
    return    

##########################
#Setup Main
##########################

if __name__=="__main__":
    bedFile = sys.argv[1]
    lincLotID = sys.argv[2]
    main(bedFile,lincLotID)
