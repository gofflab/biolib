'''
Created on Jun 23, 2011

@author: lgoff
'''
from pygr import annotation, mapping
from pygr import worldbase

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
		self.id = cols[0] # sequence ID
		self.start = int(cols[3]) - 1 # correct for 1-based coords
		self.stop = int(cols[4])
		if cols[6] == '+': # convert to Pygr convention
			self.orientation = 1
		elif cols[6] == '-':
			self.orientation = -1
		else:
			raise ValueError('Bad strand: %s' % cols[6])
		for s in cols[8].split(';')[:-1]: # parse attributes
			attr, val = s.strip().split(' ')
			#print '%s: %s' % (attr,val)
			if ',' in val:
				setattr(self, attr, val.split(','))
			else:
				setattr(self, attr, val)

def read_gff3(filename, genome):
	d = {} # for different types of sliceDBs
	ifile = file(filename)
	for line in ifile: # parse all the GFF3 lines
		if line.startswith('#'): # ignore this line
			continue
		row = GFF3Row(line)
		try:
			d.setdefault(row.type, {})[row.gene_id] = row
		except AttributeError:
			pass # no type or ID so ignore...
	ifile.close()
	annotations = {}
	for atype,sliceDB in d.items(): # create annotation DBs
		adb = annotation.AnnotationDB(sliceDB, genome)
		annotations[atype] = adb
	return annotations


#from pygr import cnestedlist,seqdb
#import glob
#
#mafdir = "/n/rinn_data1/indexes/human/hg19/alignments/hg19_ucsc_multiz46way/maf/unzipped"
#
#mafFiles = glob.glob(mafdir+"/*.maf")
#
#genomes = {'hg19':seqdb.SequenceFileDB('/n/rinn_data1/indexes/human/hg19/hg19.fa'),
#            'mm9':seqdb.SequenceFileDB('/n/rinn_data1/indexes/igenomes/Mus_musculus/UCSC/mm9/Sequence/Chromosomes/mm9.fa')
#}
#
#genomeUnion=seqdb.PrefixUnionDict(genomes)
#al = cnestedlist.NLMSA('hg19_vs_mm9','w',genomeUnion,mafFiles = mafFiles)

from pygr import cnestedlist

msa = cnestedlist.NLMSA('hg19_vs_mm9','r')

ival = msa.seqDict['hg19.chr7'][27180996:27183287] #HOXA5 in human

for x in msa[ival]:
	print repr(x)
#
#  OR
#  
for x,y,e in msa[ival].edges():
	print "%s\t%s\t%s\n%s\t%s\t%s\n" % (x,(~(msa.seqDict))[x],repr(x),y,(~(msa.seqDict))[y],repr(y))


