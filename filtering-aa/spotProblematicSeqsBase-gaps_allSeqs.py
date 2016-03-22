#! /usr/bin/env python
import sys

alnNamesFile=sys.argv[1]
outprefix="satePrankOrthologs%s_allSeqs"%(alnNamesFile.split(".")[1])
aaMatrix="blosum62.txt"
#ignoreSeqs=["HUMAN", "ANOCA","CHEMY","ALLIG"]
ignoreSeqs=[]
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

out1=open("%s_inserts_allSeqs.txt"%(outprefix), "w")

for filename in alnNames:

    wrong={} #store all windows for this alignment
    
    filename=filename.rstrip()
    #fileLocation="8000orthologs/%s/sate.prank.pep.aligned"%(filename)
    
    #remove human and lizard seqs if any
    #extractSeqs(fileLocation,outprefix,ignoreSeqs) # writes clean.temp in current dir without ignoreSeqs

#    align = AlignIO.read("%s_clean.temp"%(outprefix), "fasta")
    align = AlignIO.read("8000orthologs/%s/sate.prank.pep.aligned"%(filename), "fasta")

    oneSeqInsertion={}
    twoSeqInsertion={}
    threeSeqInsertion={}
    
    for entry in align:
        name=entry.id
        oneSeqInsertion[name]=[]
        twoSeqInsertion[name]=[]
        threeSeqInsertion[name]=[]

    j=0
    while j<len(align[1].seq):
        
        position=j+1

        col=align[:,j] 

        cleancol=col.replace("-", "")

        if len(cleancol)==1:

            for name in align:
                species=name.id
                currentSite=name.seq[j]
                
                seq=str(name.seq)
                siteSeq=len(re.sub("-", "", seq[:j]))+1               
            
                if currentSite!="-":
                    if len(oneSeqInsertion[species])>0:
#                        if oneSeqInsertion[species][-1][1]+1==position:
#                            oneSeqInsertion[species][-1][1]=position
#                            oneSeqInsertion[species][-1][3]=siteSeq
                        if oneSeqInsertion[species][-1][1]+1==siteSeq:
                            oneSeqInsertion[species][-1][3]=position
                            oneSeqInsertion[species][-1][1]=siteSeq    
                        else:
                            oneSeqInsertion[species].append([siteSeq,siteSeq, position, position])
                    else:
                        oneSeqInsertion[species].append([siteSeq,siteSeq, position, position])
                            
        if len(cleancol)==2:

            for name in align:
                species=name.id
                currentSite=name.seq[j]
                
                seq=str(name.seq)
                siteSeq=len(re.sub("-", "", seq[:j]))+1      
                
                if currentSite!="-":
                    if len(twoSeqInsertion[species])>0:
                        if twoSeqInsertion[species][-1][1]+1==siteSeq:
                            twoSeqInsertion[species][-1][1]=siteSeq
                            twoSeqInsertion[species][-1][3]=position
                        else:
                            twoSeqInsertion[species].append([siteSeq,siteSeq, position, position])   
                    else:
                        twoSeqInsertion[species].append([siteSeq,siteSeq, position, position])     

        if len(cleancol)==3:

            for name in align:
                species=name.id
                currentSite=name.seq[j]
                
                seq=str(name.seq)
                siteSeq=len(re.sub("-", "", seq[:j]))+1      
                
                if currentSite!="-":
                    if len(threeSeqInsertion[species])>0:
                        if threeSeqInsertion[species][-1][1]+1==siteSeq:
                            threeSeqInsertion[species][-1][1]=siteSeq
                            threeSeqInsertion[species][-1][3]=position
                        else:
                            threeSeqInsertion[species].append([siteSeq,siteSeq, position, position])   
                    else:
                        threeSeqInsertion[species].append([siteSeq,siteSeq, position, position])    
        j+=1
            
    for name in oneSeqInsertion:
        if len(oneSeqInsertion[name])>0:
            
            info=oneSeqInsertion[name]
            
            i=0
            while i <len(info):
                infoA=info[i]

                a="\t".join(["%s" % el for el in infoA])      
                out1.write("%s\t%s\t%s\tI1\n"%(filename, name, a))
                
                i+=1

    for name in twoSeqInsertion:
        if len(twoSeqInsertion[name])>0:
            
            info=twoSeqInsertion[name]
            
            i=0
            while i <len(info):
                infoA=info[i]

                a="\t".join(["%s" % el for el in infoA])      
                out1.write("%s\t%s\t%s\tI2\n"%(filename, name, a))
                
                i+=1
                
    for name in threeSeqInsertion:
        if len(threeSeqInsertion[name])>0:
            
            info=threeSeqInsertion[name]
            
            i=0
            while i <len(info):
                infoA=info[i]

                a="\t".join(["%s" % el for el in infoA])      
                out1.write("%s\t%s\t%s\tI3\n"%(filename, name, a))
                
                i+=1
                
out1.close()

