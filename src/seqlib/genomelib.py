'''
Created on Aug 28, 2010

This is a port of the genome.py module from seqtools (it is a work in progress)

@author: lgoff
'''
############
#Imports
############
import sequencelib
import random
from pygr import seqdb, sqlgraph, annotation, worldbase, cnestedlist
import sys
#######
#Constants
#######

purines=['A','G']
pyrimidines=['C','T','U']

chr_names = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
             'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19',
             'chr20','chr21','chr22','chrX','chrY']

genome_length = 3080419480

chr_lengths = {'chr1':247249719,
               'chr2':242951149,
               'chr3':199501827,
               'chr4':191273063,
               'chr5':180857866,
               'chr6':170899992,
               'chr7':158821424,
               'chr8':146274826,
               'chr9':140273252,
               'chr10':135374737,
               'chr11':134452384,
               'chr12':132349534,
               'chr13':114142980,
               'chr14':106368585,
               'chr15':100338915,
               'chr16':88827254,
               'chr17':78774742,
               'chr18':76117153,
               'chr19':63811651,
               'chr20':62435964,
               'chr21':46944323,
               'chr22':49691432,
               'chrX':154913754,
               'chrY':57772954
               }

genbases = {'A': 843953565, 'C': 584268578, 'T': 845168978, 'G': 584621685, 'N': 222406671}
genfreqs = {'A': 0.27397358394837834, 'C': 0.18967175795161509, 'T': 0.27436814482162669, 'G': 0.18978638746954035, 'N': 0.072200124834946186}

###############
#BROAD SETTINGS
###############
#genome_build = 'hg18'
#genome_dir = '/seq/compbio-hp/lgoff/genomes/'+genome_build
#genome_file = genome_build+".fa"
#hg19_genome_file = '/fg/compbio-t/lgoff/magda/references/human/genome/hg19/hg19.fa'
#hg18_genome_file = '/fg/compbio-t/lgoff/magda/references/human/genome/hg18/hg18.fa'
#mm9_genome_file = '/fg/compbio-t/lgoff/magda/references/mouse/genome/mm9/mm9.fa'
#rmgenome_dir = "/seq/compbio-hp/lgoff/smallRNAs/genomes/human_repeatmasked/"
#
#mammals_alignments_dir = '/ahg/scr3/mammals/ucsc/multiz44way/'

################
#Valor Settings
################
genome_build = 'hg18'
genome_dir = '/n/rinn_data1/indexes/human/'+genome_build
genome_file = genome_build+".fa"
hg19_genome_file = '/n/rinn_data1/indexes/human/hg19/hg19.fa'
hg18_genome_file = '/n/rinn_data1/indexes/human/hg18/hg18.fa'
mm9_genome_file = '/n/rinn_data1/indexes/igenomes/Mus_musculus/UCSC/mm9/Sequence/Chromosomes/mm9.fa'
#rmgenome_dir = "/seq/compbio-hp/lgoff/smallRNAs/genomes/human_repeatmasked/"

#mammals_alignments_dir = '/ahg/scr3/mammals/ucsc/multiz44way/'


bed_fields = ['chr','start','end','label','score','strand']
#######
#Functions
#######
def fetch_genbases(genhandle,genbases={}):
    bases = ['A','T','G','C','N']
    geniter = sequencelib.FastaIterator(genhandle)
    for genseq in geniter:
        print genseq['name']
        seq = genseq['sequence'].upper()
        for b in bases:
            genbases[b] = seq.count(b) + genbases.get(b,0)
    return genbases

def fetch_genome_freqs():
    """Specifically returns a dictionary containing frequencies of every 7mer in hg18"""
    freqfile = '/seq/compbio-hp/lgoff/smallRNAs/genomes/human/hg18/hg18_7mer_frequencies.txt'
    freqhandle = open(freqfile,'r')
    freqs = {}
    for line in freqhandle:
        vals = line.rstrip().split()
        freqs[vals[0]] = float(vals[1])
    return freqs


def random_region(n,m=1):
    '''Generate a random region of max length "n" and min length "m" (default m=1).'''
    c = random.choice(chr_names)
    strand= random.choice(["+","-"])
    start = random.randint(1,chr_lengths[c])
    end = start+random.randint(m,n)
    return c, start, end, strand

def isMasked(s):
    maskedChars='actgnN'
    for c in s:
        if c in maskedChars:
            return True
    return False


#######################
#pygr specific
#######################
#SeqPath = pygr.Data.Bio.Seq.Genome.HUMAN.hg18

def pygrConnect(genome="hg18",useWorldbase = False):
    if useWorldbase:
        if genome == "hg18":
            res=worldbase.Bio.Seq.Genome.HUMAN.hg18()
        elif genome == "hg19":
            res=worldbase.Bio.Seq.Genome.HUMAN.hg19()
        elif genome == "mm9":
            res=worldbase.Bio.Seq.Genome.MOUSE.mm9()
        elif genome == "mm8":
            res=worldbase.Bio.Seq.Genome.MOUSE.mm8()
        else:
            raise AssertionError ("No genome by that name in worldbase. (that I'm currently aware of...)")
    else:
        if genome == "hg18":
            res = seqdb.SequenceFileDB(hg18_genome_file)
        elif genome == "hg19":
            res = seqdb.SequenceFileDB(hg19_genome_file)
        elif genome == "mm9":
            res = seqdb.SequenceFileDB(mm9_genome_file)
        else:
            raise AssertionError ("I'm not sure how to handle that genome build yet...sorry. Please create a seqquenceFileDB for this genome.")
    return res

#pygr annotation layers
#This is very closely tied to valor
class UCSCStrandDescr(object):
    def __get__(self, obj, objtype):
        if obj.strand == '+':
            return 1
        else:
            return -1

class UCSCSeqIntervalRow(sqlgraph.TupleO):
    orientation = UCSCStrandDescr()

serverInfo = sqlgraph.DBServerInfo(host='localhost',user='root',passwd='')

def build_rmsk_nlmsa(genome="hg19"):
    #This is horse shit...
    
    seqDB = pygrConnect(genome)
    rmsk = sqlgraph.SQLTable('hg19.rmsk',serverInfo=serverInfo,itemClass=UCSCSeqIntervalRow,primaryKey="lookupName")
    annodb = annotation.AnnotationDB(rmsk,
                                     seqDB,
                                     sliceAttrDict=dict(id='genoName',
                                                        start='genoStart',
                                                        stop='genoEnd',
                                                        orientation='orientation'
                                                        ),
                                     annotationType='repeat:')
    al = cnestedlist.NLMSA('/n/rinn_data1/indexes/human/'+genome+'/repeat_'+genome,'w',pairwiseMode=True)
    for k in annodb:
        al.addAnnotation(annodb[k])
    al.build()

def refGene_nlmsa(genome="hg19"):
    #Needed to add primary key 'lookupName' to hg19.refGene for this to work (pygr requires unique ids for an annotation)
    #This is really CRAP....I don't know how or why anyone will every be able to use this....
    
    try:
        al = cnestedlist.NLMSA('/n/rinn_data1/indexes/human/'+genome+'/refGene/refGene_'+genome,'r')
    except:
        sys.stderr.write("Could not find NLMSA index, attempting to build one...\n")
        seqDB = pygrConnect(genome)
        sys.stderr.write("Found genome...\n")
        refGene = sqlgraph.SQLTable('hg19.refGene',serverInfo=serverInfo,itemClass=UCSCSeqIntervalRow,primaryKey="lookupName")
        sys.stderr.write("Got table from Valor UCSC...\n")
        annodb = annotation.AnnotationDB(refGene,
                                         seqDB,
                                         sliceAttrDict=dict(id='chrom',
                                                            start='txStart',
                                                            stop='txEnd',
                                                            orientation='orientation'
                                                            ),
                                         annotationType='refGene:')
        sys.stderr.write("annodb created...\n")
        sys.stderr.write('Creating NLMSA object at /n/rinn_data1/indexes/human/'+genome+'/refGene/refGene_'+genome+'...\n')
        al = cnestedlist.NLMSA('/n/rinn_data1/indexes/human/'+genome+'/refGene/refGene_'+genome,'w',pairwiseMode=True)
        for k in annodb:
            al.addAnnotation(annodb[k])
        al.build(saveSeqDict=True)
        sys.stderr.write("Done!\n")
    return al

################
#MISC
################
def fetchSequence(chrom,start,end,strand,genome="hg18"):
    connection=pygrConnect(genome)
    start,end=int(start),int(end)
    seq=connection[chrom][start:end]
    if strand == "-":
        seq=-seq
    return seq
