from numpy import *

def getAveProbPerWindowPerSeq(align, pGaps):
    
    scorePerSpeciesPerWindow={} # [species]=[Z-score w1, Z-score w2...]
    ####################################################################################
    # amino acid frequencies, gaps and X will be given a zero prob
    aminoacids=["X","A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "-"]

    #start up the aa freqs dicts
    aaFreqsDict={}          
    for element in aminoacids:

        aaFreqsDict[element]=[] 
    

    currentWindowSpeciesProbs={} # [species]=[probsite1, probsite2...]
    aveWindowSpeciesProbs={} # [species]=windowprob
    window=len(align[1].seq)
    
    i=0
    while i<window:
        
        col=align[:, i]
        cleancol=col.replace("-", "").replace("X", "")
        
        # populate the amino acids frequencies per site per column 
        for key in aaFreqsDict:

            if key!="-" and key!="X" and len(cleancol)>0 :
                aaFreqsDict[key]=round(cleancol.count(key)*1.00/len(cleancol), 10)

            else:
                aaFreqsDict[key]="NA"
                
        # save the prob of a species site aa
        for name in align:
            species=name.id
            currentAA=name.seq[i]
            siteProb=aaFreqsDict[currentAA]
            if not currentWindowSpeciesProbs.has_key(species):
                currentWindowSpeciesProbs[species]=[]        
                
#            if siteProb!="NA":#
            currentWindowSpeciesProbs[species].append(siteProb)
#            else:
#                currentWindowSpeciesProbs[species].append(siteProb)#####
                
        i+=1
        
    #calculate ave prob per window  per seq
    allAveWindow=[]
    
    for species in currentWindowSpeciesProbs:
        probs=currentWindowSpeciesProbs[species]
        probsNumeric=[x for x in probs if x!="NA" ]
        if len(probsNumeric)>pGaps*window: #if not mostly-gap seq
        
            aveProb=array(probsNumeric).mean()
            aveWindowSpeciesProbs[species]=round(aveProb, 2)
            allAveWindow.append(aveProb)
            
        else:
            aveWindowSpeciesProbs[species]="NA"
            
    # calculate Z-score
    if allAveWindow>0: #not a gap-only window
        aveWindow=array(allAveWindow).mean()
        stdevWindow=array(allAveWindow).std()
        
    for species in aveWindowSpeciesProbs:
        probSpecies=aveWindowSpeciesProbs[species]
        Z=0
        if not probSpecies=="NA": #not a gap-only window
            if not stdevWindow<0.001: #very conserved window
                dev=abs(aveWindow-probSpecies)
                Z=round(dev/stdevWindow, 2)
            else:
                Z=0
        else:
            Z="NA"
            
        if not scorePerSpeciesPerWindow.has_key(species):
            scorePerSpeciesPerWindow[species]=[]
            
        score=[probSpecies,Z]
        scorePerSpeciesPerWindow[species]=score
      
    return scorePerSpeciesPerWindow

def readMatrix(nameMatrix):
    #takes an half similarity matrix and turns into a dict

    input=open(nameMatrix)
    lines=input.readlines()

    names=lines[0].rstrip().split()

    aaDict={} #aaDict[(aa1,aa2)]=value

    i=1
    while i<len(lines):
        info=lines[i].rstrip().split()
        aa1=info[0]

        j=1
        while j<len(info):
            aa2=names[j-1]
            value=int(info[j])
            pair=(aa1,aa2)
            aaDict[pair]=value
            pairI=(aa2,aa1)
            aaDict[pairI]=value
            j+=1
            
        i+=1
    return aaDict

def scorePairAlignedSeqsWithMatrix(seq1,seq2,matrixAA):
    #takes 2 aa seqs as strings with the same size, 
    #and calculates the score of aligned aas ignoring gaps

    totalScore=[]
    gaps=0
    i=0
    while i<len(seq1):
        if seq1[i]!="-" and seq2[i]!="-" and seq1[i]!="X" and seq2[i]!="X":
            totalScore.append(matrixAA[(seq1[i],seq2[i])])
        else:
            gaps+=1
        i+=1
        
    if gaps==len(seq1):
        
        total="NA"
    
    else:
        total=array(totalScore).mean()
    
    return total

def getPairwiseScore(windowAlign, matrixAA):
    
    scoreSlice={} # [species]=lowestScore
    allSeqs={}
    forAverage=[]
    
    for entry in windowAlign:
        name=entry.id
        allSeqs[name]=str(entry.seq)

    for name in allSeqs:
        allScores=[]
        others=allSeqs.keys()
        others.remove(name)
        
        for alt in others:
            seq1=allSeqs[name]
            seq2=allSeqs[alt]
            score=scorePairAlignedSeqsWithMatrix(seq1, seq2, matrixAA)
            
            if score!="NA":
                allScores.append(score)
            
        if len(allScores)>0:
            score=round(allScores[-1], 2)
            forAverage.append(score)
            scoreSlice[name]=[score]
        
    aveWindow=array(forAverage).mean()
    stdevWindow=array(forAverage).std()
        
    for species in allSeqs:
        if scoreSlice.has_key(species):
            score=round(scoreSlice[species][0], 2)
            
            if not stdevWindow<0.001: #very conserved window
                dev=abs(aveWindow-score)
                Z=round(dev/stdevWindow, 2)
                
            else:
                Z=0
                
            scoreSlice[species].append(Z)
        else:
            scoreSlice[species]=["NA", "NA"]
        
    return scoreSlice

def extractSeqs(fileLocation, alnName, listSpecies):
    
    import re
    seqs={}
    ids={}

    align = open("%s"%(fileLocation))
    info=align.readline()
    while info:
        if re.match(">", info):
            key=info.rstrip()[1:]
            seqs[key]=""
            
            species=key.split("_")[-1]
            ids[species]=key
            
        else:
            seqs[key]+=info.rstrip()
            
        info=align.readline()
        
    for name in listSpecies:
        if ids.has_key(name):
            seqs.pop(ids[name])
            
    out1=open("%s_clean.temp"%(alnName), "w")
    for key in seqs:
        out1.write(">%s\n%s\n"%(key, seqs[key]))
    out1.close()

def getColumnDistScoreToClosest(windowAlign, matrixAA):
    
    scoreSlice={} # [species]=lowestScore
    allSeqs={}
    forAverage=[]
    allScores=[]
    
    for entry in windowAlign:
        name=entry.id
        allSeqs[name]=str(entry.seq)

    for name in allSeqs:
        seq1=allSeqs[name]
        others=allSeqs.keys()
        others.remove(name)
        windowScore=[]
        
        i=0
        while i<len(seq1):
            scoreCol=[]
            #score each position of the alignment 
            for alt in others:
            
                seq2=allSeqs[alt][i]
                score=scorePairAlignedSeqsWithMatrix(seq1[i], seq2, matrixAA)
                if score!="NA":
                    scoreCol.append(score)
            
            if len(scoreCol)>0: 
                scoreCol.sort()
                windowScore.append(scoreCol[-1]) #get the highest score for this column  
                
            i+=1
        
        if len(windowScore)>0:
            scoreAVEwindow=array(windowScore).mean()
            scoreAVEwindow=round(scoreAVEwindow, 2)
            
            allScores.append(scoreAVEwindow)
            scoreSlice[name]=[scoreAVEwindow]   
    
    aveWindow=array(allScores).mean()
    stdevWindow=array(allScores).std()
    
    for species in allSeqs:
        if scoreSlice.has_key(species):
            score=scoreSlice[species][0]
            
            if not stdevWindow<0.001: #very conserved window
                dev=abs(aveWindow-score)
                Z=round(dev/stdevWindow, 2)
            else:
                Z=0
                
            scoreSlice[species].append(Z)
        else:
            scoreSlice[species]=["NA", "NA"]
        
    return scoreSlice
