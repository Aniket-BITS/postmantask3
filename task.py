#!/usr/bin/env python
# coding: utf-8

# In[2]:


import sys
import os
def read_sequences(file1, file2):
    with open(file1, 'r') as f1, open(file2, 'r') as f2:
        seq1 = f1.readline().strip()
        seq2 = f2.readline().strip()
    return seq1, seq2
def global_alignment(seq1, seq2, match_score=1, mismatch_score=-1, gap_penalty=1):
    #seq1,seq2=read_sequences('sample1','sample2')
    #seq1,seq2='CGTGAATTCAT','GACTTAC'
    am=[[0]]
    for i in seq2:
        am.append([0])
    #filling up MATRIX
    #introducing scoring system
    g=gap_penalty
    mis=mismatch_score
    mch=match_score
    def filler(h,v,d,sim=0):
        h=h+g
        v=v+g
        if sim==1:
            d=d+mch
        else:
            d=d+mis
        return max(h,v,d)
    for j in range(1,len(seq1)+1):
        am[0].append(g*j)
    for i in range(1,len(seq2)+1):
        am[i][0]=(g*i)
    for i in range(1,len(seq2)+1):
        for j in range(1,len(seq1)+1):
            if seq1[j-1]==seq2[i-1]:
                sim=1 #matching case
            else:
                sim=0
            am[i].append(filler(am[i][j-1],am[i-1][j],am[i-1][j-1],sim))
    #tracing back
    j=len(seq1)
    i=len(seq2)
    alignedseq1=""
    alignedseq2=""
    alignmentscore=0
    while i > 0 or j > 0:
        if j > 0 and am[i][j]==am[i][j-1]+g:
            alignedseq1= seq1[j - 1] + alignedseq1
            alignedseq2= "-" + alignedseq2
            j -= 1
        elif i > 0 and am[i][j] ==am[i-1][j] + g:
            alignedseq1 = "-" + alignedseq1
            alignedseq2 = seq2[i - 1] + alignedseq2
            i -= 1
        else:
            alignedseq1 = seq1[j - 1] + alignedseq1
            alignedseq2 = seq2[i - 1] + alignedseq2
            i -= 1
            j -= 1
    return alignedseq1,alignedseq2
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python global_alignment.py <sequence_file1> <sequence_file2>")
        sys.exit(1)

    file1 = sys.argv[1]
    file2 = sys.argv[2]

    seq1, seq2 = read_sequences(file1, file2)

    match_score = 1
    mismatch_score = -1
    gap_penalty = -1

    alignedseq1, alignedseq2 = global_alignment(seq1, seq2, match_score, mismatch_score, gap_penalty)

    print("Aligned Sequence 1:", alignedseq1)
    print("Aligned Sequence 2:", alignedseq2)


# In[ ]:




