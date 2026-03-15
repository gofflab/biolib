"""Plotting utilities for genomic and epigenomic data visualisation.

Provides helper functions for generating publication-quality plots of
chromatin mark occupancy and other aggregate genomic features using R via
Rscript.
"""
import os


def chromatinAggPlots(basename):
    """Generates chromatin aggregate plots as a multi-panel PDF using R.

    Writes an R script that reads three data files produced by an upstream
    pipeline step, then calls Rscript to execute it and produce a PDF of
    aggregate chromatin mark occupancy profiles centred on smRNA predictions.

    Required input files (all derived from basename):
        - basename.vec: Tab-delimited matrix of signal values.
        - basename.row: Tab-delimited BED-like annotation of rows.
        - basename.col: Tab-delimited column name file.

    Output:
        - basename.pdf: Multi-panel PDF with one line plot per chromatin mark.
        - basename.q: The R script used to generate the plot (retained).

    Args:
        basename: Base path/name shared by all input files and used for the
            output PDF and R script.

    Returns:
        The return code of the Rscript invocation (0 on success).
    """
    myScript = """
colNames<-read.table("%s.col",colClasses="character",header=F,sep="\\t")
myCols<-unlist(colNames)

myNames<-read.table("%s.row",sep="\\t",header=F)
colnames(myNames)<-c("chr","start","end","name","score","strand")

data<-read.table("%s.vec",sep="\\t",header=F)
colnames(data)<-myCols

myColSplit<-unlist(strsplit(myCols,"_"))
myIndex<-cbind(myColSplit[seq(1,length(myColSplit),2)],myColSplit[seq(2,length(myColSplit),2)])

myIndex[,2]<-gsub("\\\\.","-",myIndex[,2])
myIndex<-as.data.frame(myIndex)
myIndex[,2]<-as.integer(format(myIndex[,2]))

marks<-unique(myIndex[,1])
pos<-unique(myIndex[,2])

#Set colors
cols<-c("green","lightgreen","blue","red","darkred","orange")

################
#Print chromatin aggregation plot data vector
################

pdf("%s.pdf",width=11,height=8.5)
par(mfrow=c(2,3))

for (i in 1:length(marks)) {
    plot(pos,data[myIndex[,1]==marks[i]]/min(data[myIndex[,1]==marks[i]]),col=cols[i],type="l",lwd="1.5",ylab="Fold change in mark occupancy relative to min",xlab="Position relative to start of smRNA prediction",ylim=c(1,2.4),main=paste("H1 flanking density of ",marks[i]))  
    abline(v=0,lty="dashed",col="lightgrey")
   }
dev.off()
    
    """ % (basename,basename,basename,basename)
    handle=open(basename+".q",'w')
    handle.write(myScript)
    handle.close()
    myCommand = """Rscript --vanilla %s.q""" % basename
    res = os.system(myCommand)
    return res
