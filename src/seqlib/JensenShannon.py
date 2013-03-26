#!/usr/bin/env python

"""
JensenShannon.py

Created by Loyal Goff on Nov 10, 2010.
Copyright (c) 2010
"""
from scipy import *
from numpy import *
import time
from scipy.stats.distributions import entropy
import rpy2.robjects as r
import rpy2.robjects.numpy2ri
#efficnent js_div
def js_div_matrix(a):
    a=array(a)
    W=zeros((a.shape[0],a.shape[0]))
    e=-entropy(a.transpose())
    for i in range(a.shape[0]):
        val_range=range(i+1,a.shape[0])
        sumAB=tile(a[i,:],(a.shape[0]-i-1,1))+a[val_range,:]
        result=0.5*(e[i]+e[val_range,:]-sum((sumAB)*nan_to_num(log(((sumAB)/2))),1))
        W[val_range,i]=result
        W[i,val_range]=result
    return W

def make_probs(a):
    sums = sum(a,1)
    res = zeros(a.shape)
    for i in xrange(a.shape[0]):
        res[i,:]=a[i,:]/sums[i]
    return res

def js_div(A,B):
    half=(A+B)/2
    return 0.5*kl_div(A,half)+0.5*kl_div(B,half)

def kl_div(A,B):
    return sum(multiply(A,log(A/B)))

def main():
    pass

if __name__ == "__main__":
    fname = '/seq/rinnscratch/cole/linc_db/differential_expr/lincs_all_data/genes.fpkm_tracking'
    colLabs = ['HUVEC','NHLF','HepG2','H1-hESC','NHEK','HSMM','K562','MCF7','GM12878','HeLa','Lung','Adipose','Prostate','Adrenal','Testes','Thyroid','Colon','Kidney','Breast','Lymph','Brain','Liver','Ovary','Skel_Musc','Heart','Leukocytes','Placenta']
    myCols = range(6,87,3)
    a = loadtxt(fname,dtype="float",skiprows=1,usecols=myCols)
    sums = sum(a,1)
    a = a[sums>0,:]
    #a = random.random((200,30))
    #a[0,10] = 0.0
    #a[2,15] = 0.0
    #a[178,22] = 0.0
    #a[178,2] = 0.0
    #a[178,11] = 0.0
    #a = a[:2000,:]
    
#    r.r.pdf('isoform_row_JS.pdf')
    #Rows
#    rowMat = make_probs(a)
#    rowJS = js_div_matrix(rowMat)
#    rowJS_dist = sqrt(rowJS)
#    rowDist = r.r['as.dist'](rowJS_dist)
#    rowHclust = r.r.hclust(rowDist)
#    rowDendro = r.r['as.dendrogram'](rowHclust)
#    r.r.plot(rowHclust,main='',xlab='',ylab='JS-distance')
#   r.r['dev.off']()
    
    
    r.r.pdf('isoform_column_JS.pdf')
    #Columns
    #colMat = log(a[sum(a,1)>0,]+1).transpose()
    colMat = a[sum(a,1)>0,].transpose()
    #colMat = a.transpose()
    colMat = make_probs(colMat)
    print colMat[1:5,1:5]
    colJS = js_div_matrix(colMat)
    print colJS
    colJS_dist = sqrt(colJS)
    
    colDist = r.r['as.dist'](colJS_dist)
    colHclust = r.r.hclust(colDist)
    colHclust[3] = colLabs
    colDendro = r.r['as.dendrogram'](colHclust)
    r.r.plot(colHclust,main="JS Distance",xlab="",sub="",ylab="JS-distance on FPKM")

#    
#    #colMat = a[sum(a,1)>0,].transpose()
#    #coldist = r.r.dist(r.r.log2(colMat+0.001))
#    coldist = r.r.dist(colMat)
#    colHclust = r.r.hclust(coldist)
#    colHclust[3] = colLabs
#    colDendro = r.r['as.dendrogram'](colHclust)
#
#    r.r.plot(colHclust,main="Euclidean",sub="",xlab="",ylab="Euclidean-distance on log2 FPKM")
#    
    
    colcor = r.r.cor(colMat.transpose())
    #print colcor
    colcor = 1-(array(colcor)**2)
    #print colcor
    coldist = r.r['as.dist'](colcor)
    colHclust = r.r.hclust(coldist)
    colHclust[3] = colLabs
    colDendro = r.r['as.dendrogram'](colHclust)
    #print '%s took %0.3f ms' % (js_div_matrix.func_name, (t2-t1)*1000.0)
    r.r.plot(colHclust,main="Pearson",sub="",xlab="",ylab="Pearson-distance on FPKM")
    #heatmap
    
    r.r['dev.off']()