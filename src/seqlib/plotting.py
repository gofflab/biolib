'''
Created on Jul 13, 2010

@author: lgoff
'''
import os

def chromatinAggPlots(basename):
    """
    Makes chromatin aggregate plots
    
    requires:
        basename.vec
        basename.row
        basename.col
        
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