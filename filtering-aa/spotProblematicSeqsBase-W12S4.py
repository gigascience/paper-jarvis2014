#! /usr/bin/env python
import sys
window=12
step=1
maxGapContent=0.6
alnNamesFile=sys.argv[1]
outprefix="satePrankOrthologs.W12S1G6_%s"%(alnNamesFile.split(".")[1])
aaMatrix="blosum62.txt"
ignoreSeqs=["HUMAN", "ANOCA","CHEMY","ALLIG"]

#fileLocation CHANGE LINE 29 

########################################################################
from numpy import *
from spotProblematicSeqsModules import *
from Bio import AlignIO
import re
import os

alnNames=open(alnNamesFile)
alnNames=alnNames.readlines()

matrixAA=readMatrix(aaMatrix)

out1=open("%s_WrongA.txt"%(outprefix), "w")
out2=open("%s_WrongC.txt"%(outprefix), "w")

for filename in alnNames:

    wrong={} #store all windows for this alignment
    
    filename=filename.rstrip()
    fileLocation="8000orthologs/%s/sate.prank.pep.aligned"%(filename)
    
    #remove human and lizard seqs if any
    extractSeqs(fileLocation,outprefix,ignoreSeqs) # writes clean.temp in current dir without ignoreSeqs

    align = AlignIO.read("%s_clean.temp"%(outprefix), "fasta")
    
    results={}
    j=0
    while j<len(align[1].seq)-window:
        
        start=j+1
        end=j+window

        windowAlign=align[:, j:j+window] 
        scoresProb=getAveProbPerWindowPerSeq(windowAlign,maxGapContent)
#        scoresMatrix=getPairwiseScore(windowAlign, matrixAA)
        scoresMatrixPerCol=getColumnDistScoreToClosest(windowAlign, matrixAA) 
        
        for record in align:
            name=record.id 
            species=name.split("_")[-1]
            
            seq=str(record.seq)
            startSeq=len(re.sub("-", "", seq[:j]))+1
            endSeq=len(re.sub("-", "", seq[:j+window]))
        
            id=name.split("_")[-1]
            
            P=scoresProb[name][0]
            Z=scoresProb[name][1]
            
            D=scoresMatrixPerCol[name][0]
            Z62=scoresMatrixPerCol[name][1]
            
            if P!="NA" and float(P)<0.7: #ignore windows with too many gaps or those that are very conserved
            
                #model for ALL in W15S3 but actually works better also with W12
                result=1+2.21*P+0.36*Z62+1.32*Z+0.06*D+(-1.5)*D*P+(-0.70)*Z*P+(-0.10)*Z*Z62                
                
                if result > 0: # it's WRONG or at least very divergent according to the alignment
                    
                    if not wrong.has_key(name):
                        wrong[name]=[]
                
                    if len(wrong[name])>0:

                        if wrong[name][-1][1]>=startSeq: #check for coordinates overlap

                            wrong[name][-1][1]=endSeq
                            wrong[name][-1][3]=end
                            wrong[name][-1][-1]+=1 #n.o window overlapped
                            wrong[name][-1][-2]+=result
                        else:
                            wrong[name].append([startSeq, endSeq, start,end,result,1])
                    else:
                        wrong[name].append([startSeq, endSeq, start,end,result,1])
        
        j+=step
        
    for name in wrong:
        info=wrong[name]
        
        i=0
        while i <len(info):
            infoA=info[i]
            ave=infoA[4]
            nWindows=infoA[5]
            aveResults=round(ave/nWindows, 2)
            infoA[-2]=aveResults
            
            if aveResults > 3:
                
                infoA[-1]="C"
                a="\t".join(["%s" % el for el in infoA])      
                out1.write("%s\t%s\t%s\n"%(filename, name, a))
            
            else:
            
                infoA[-1]="A"
                a="\t".join(["%s" % el for el in infoA])      
                out2.write("%s\t%s\t%s\n"%(filename, name, a))
            
            i+=1

out1.close()
out2.close()

