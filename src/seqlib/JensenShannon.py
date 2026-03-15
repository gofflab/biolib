#!/usr/bin/env python

"""Jensen-Shannon divergence utilities for comparing probability distributions.

Provides functions to compute the Jensen-Shannon (JS) divergence between pairs
of discrete probability distributions and to construct pairwise JS divergence
matrices from a collection of distributions.  The JS divergence is a
symmetrised, smoothed version of the Kullback-Leibler (KL) divergence and is
defined as::

    JS(A || B) = 0.5 * KL(A || M) + 0.5 * KL(B || M),  where M = (A + B) / 2

Because it is bounded in [0, ln 2] (or [0, 1] in bits), its square root is a
proper metric known as the Jensen-Shannon distance.

Originally created by Loyal Goff on Nov 10, 2010.
"""

import rpy2.robjects as r
from numpy import *
from scipy import *
from scipy.stats.distributions import entropy


#efficnent js_div
def js_div_matrix(a):
    """Compute a pairwise Jensen-Shannon divergence matrix efficiently.

    For each pair of rows ``i`` and ``j`` in ``a``, computes::

        JS(i, j) = 0.5 * (H(M) - 0.5*(H(i) + H(j)))

    where ``M = (a[i] + a[j]) / 2`` and ``H`` denotes Shannon entropy.
    The implementation avoids an O(n^2) loop over pairs by vectorising
    the inner computation row-by-row, giving O(n) outer iterations.

    Args:
        a: A 2-D array-like of shape ``(n, d)`` where each row is a
            probability distribution over ``d`` categories (rows should
            sum to 1 for the result to be a true JS divergence).

    Returns:
        A symmetric ``(n, n)`` NumPy array ``W`` where ``W[i, j]`` is the
        Jensen-Shannon divergence between rows ``i`` and ``j`` of ``a``.
        Diagonal entries are 0.
    """
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
    """Normalise each row of a 2-D array to sum to 1.

    Divides each row of ``a`` by its sum, converting raw counts or
    unnormalised weights into proper probability distributions.

    Args:
        a: A 2-D NumPy array of shape ``(n, d)`` with non-negative
            entries.  Each row must have a positive sum.

    Returns:
        A 2-D NumPy array of the same shape as ``a`` where each row
        sums to 1.0.
    """
    sums = sum(a,1)
    res = zeros(a.shape)
    for i in range(a.shape[0]):
        res[i,:]=a[i,:]/sums[i]
    return res

def js_div(A,B):
    """Compute the Jensen-Shannon divergence between two distributions.

    The JS divergence is defined as::

        JS(A || B) = 0.5 * KL(A || M) + 0.5 * KL(B || M)

    where ``M = (A + B) / 2`` is the mixture distribution and
    ``KL(P || Q) = sum(P * log(P / Q))``.  The result is symmetric and
    always non-negative.

    Args:
        A: A 1-D array-like representing the first probability
            distribution.  All entries should be positive for a
            well-defined result.
        B: A 1-D array-like representing the second probability
            distribution of the same length as ``A``.

    Returns:
        A scalar float equal to the Jensen-Shannon divergence between
        ``A`` and ``B``.
    """
    half=(A+B)/2
    return 0.5*kl_div(A,half)+0.5*kl_div(B,half)

def kl_div(A,B):
    """Compute the Kullback-Leibler divergence of distribution A from B.

    Calculates the KL divergence using the formula::

        KL(A || B) = sum(A * log(A / B))

    where the sum is taken element-wise.  The result is non-negative and
    equals zero only when ``A`` and ``B`` are identical.  Note that the
    KL divergence is not symmetric: ``kl_div(A, B) != kl_div(B, A)``
    in general.

    Args:
        A: A 1-D array-like representing the first (reference)
            probability distribution.  All entries should be positive.
        B: A 1-D array-like representing the second probability
            distribution of the same length as ``A``.  All entries
            should be positive to avoid division by zero.

    Returns:
        A scalar float equal to the KL divergence KL(A || B).
    """
    return sum(multiply(A,log(A/B)))

def main():
    """Entry point placeholder; no operation is performed."""
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
    print(colMat[1:5,1:5])
    colJS = js_div_matrix(colMat)
    print(colJS)
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
