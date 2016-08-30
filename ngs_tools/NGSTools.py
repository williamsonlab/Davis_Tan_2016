# -*- coding: utf-8 -*-
"""
Created on Sat Jun  7 14:28:28 2014

@author: jhdavis
"""

import pylab
import HTSeq
import itertools
import qMS
from string import maketrans
import operator

rRNASeq_correct   = 'GGTTAAGCGACTAAGCGTACACGGTGGATGCCCTGGCAGTCAGAGGCGATGAAGGACGTGCTAATCTGCGATAAGCGTCGGTAAGGTGATATGAACCGTTATAACCGGCGATTTCCGAATGGGGAAACCCAGTGTGTTTCGACACACTATCATTAACTGAATCCATAGGTTAATGAGGCGAACCGGGGGAACTGAAACATCTAAGTACCCCGAGGAAAAGAAATCAACCGAGATTCCCCCAGTAGCGGCGAGCGAACGGGGAGCAGCCCAGAGCCTGAATCAGTGTGTGTGTTAGTGGAAGCGTCTGGAAAGGCGTGCGATACAGGGTGACAGCCCCGTACACAAAAATGCACATGCTGTGAGCTCGATGAGTAGGGCGGGACACGTGGTATCCTGTCTGAATATGGGGGGACCATCCTCCAAGGCTAAATACTCCTGACTGACCGATAGTGAACCAGTACCGTGAGGGAAAGGCGAAAAGAACCCCGGCGAGGGGAGTGAAAAAGAACCTGAAACCGTGTACGTACAAGCAGTGGGAGCACGCTTAGGCGTGTGACTGCGTACCTTTTGTATAATGGGTCAGCGACTTATATTCTGTAGCAAGGTTAACCGAATAGGGGAGCCGAAGGGAAACCGAGTCTTAACTGGGCGTTAAGTTGCAGGGTATAGACCCGAAACCCGGTGATCTAGCCATGGGCAGGTTGAAGGTTGGGTAACACTAACTGGAGGACCGAACCGACTAATGTTGAAAAATTAGCGGATGACTTGTGGCTGGGGGTGAAAGGCCAATCAAACCGGGAGATAGCTGGTTCTCCCCGAAAGCTATTTAGGTAGCGCCTCGTGAATTCATCTCCGGGGGTAGAGCACTGTTTCGGCAAGGGGGTCATCCCGACTTACCAACCCGATGCAAACTGCGAATACCGGAGAATGTTATCACGGGAGACACACGGCGGGTGCTAACGTCCGTCGTGAAGAGGGAAACAACCCAGACCGCCAGCTAAGGTCCCAAAGTCATGGTTAAGTGGGAAACGATGTGGGAAGGCCCAGACAGCCAGGATGTTGGCTTAGAAGCAGCCATCATTTAAAGAAAGCGTAATAGCTCACTGGTCGAGTCGGCCTGCGCGGAAGATGTAACGGGGCTAAACCATGCACCGAAGCTGCGGCAGCGACACTATGTGTTGTTGGGTAGGGGAGCGTTCTGTAAGCCTGTGAAGGTGTGCTGTGAGGCATGCTGGAGGTATCAGAAGTGCGAATGCTGACATAAGTAACGATAAAGCGGGTGAAAAGCCCGCTCGCCGGAAGACCAAGGGTTCCTGTCCAACGTTAATCGGGGCAGGGTGAGTCGACCCCTAAGGCGAGGCCGAAAGGCGTAGTCGATGGGAAACAGGTTAATATTCCTGTACTTGGTGTTACTGCGAAGGGGGGACGGAGAAGGCTATGTTGGCCGGGCGACGGTTGTCCCGGTTTAAGCGTGTAGGCTGGTTTTCCAGGCAAATCCGGAAAATCAAGGCTGAGGCGTGATGACGAGGCACTACGGTGCTGAAGCAACAAATGCCCTGCTTCCAGGAAAAGCCTCTAAGCATCAGGTAACATCAAATCGTACCCCAAACCGACACAGGTGGTCAGGTAGAGAATACCAAGGCGCTTGAGAGAACTCGGGTGAAGGAACTAGGCAAAATGGTGCCGTAACTTCGGGAGAAGGCACGCTGATATGTAGGTGAAGCGACTTGCTCGTGGAGCTGAAATCAGTCGAAGATACCAGCTGGCTGCAACTGTTTATTAAAAACACAGCACTGTGCAAACACGAAAGTGGACGTATACGGTGTGACGCCTGCCCGGTGCCGGAAGGTTAATTGATGGGGTTAGCCGCAAGGCGAAGCTCTTGATCGAAGCCCCGGTAAACGGCGGCCGTAACTATAACGGTCCTAAGGTAGCGAAATTCCTTGTCGGGTAAGTTCCGACCTGCACGAATGGCGTAATGATGGCCAGGCTGTCTCCACCCGAGACTCAGTGAAATTGAACTCGCTGTGAAGATGCAGTGTACCCGCGGCAAGACGGAAAGACCCCGTGAACCTTTACTATAGCTTGACACTGAACATTGAGCCTTGATGTGTAGGATAGGTGGGAGGCTTTGAAGTGTGGACGCCAGTCTGCATGGAGCCGACCTTGAAATACCACCCTTTAATGTTTGATGTTCTAACGTTGACCCGTAATCCGGGTTGCGGACAGTGTCTGGTGGGTAGTTTGACTGGGGCGGTCTCCTCCTAAAGAGTAACGGAGGAGCACGAAGGTTGGCTAATCCTGGTCGGACATCAGGAGGTTAGTGCAATGGCATAAGCCAGCTTGACTGCGAGCGTGACGGCGCGAGCAGGTGCGAAAGCAGGTCATAGTGATCCGGTGGTTCTGAATGGAAGGGCCATCGCTCAACGGATAAAAGGTACTCCGGGGATAACAGGCTGATACCGCCCAAGAGTTCATATCGACGGCGGTGTTTGGCACCTCGATGTCGGCTCATCACATCCTGGGGCTGAAGTAGGTCCCAAGGGTATGGCTGTTCGCCATTTAAAGTGGTACGCGAGCTGGGTTTAGAACGTCGTGAGACAGTTCGGTCCCTATCTGCCGTGGGCGCTGGAGAACTGAGGGGGGCTGCTCCTAGTACGAGAGGACCGGAGTGGACGCATCACTGGTGTTCGGGTTGTCATGCCAATGGCACTGCCCGGTAGCTAAATGCGGAAGAGATAAGTGCTGAAAGCATCTAAGCACGAAACTTGCCCCGAGATGAGTTCTCCCTGACTCCTTGAGAGTCCTGAAGGAACGTTGAAGACGACGACGTTGATAGGCCGGGTGTGTAAGCGCAGCGATGCGTTGAGCTAACCGGTACTAATGAACCGTGAGGCTTAACCTT'
rRNASeq_weeks     = 'GGTTAAGCGACTAAGCGTACACGGTGGATGCCCTGGCAGTCAGAGGCGATGAAGGACGTGCTAATCTGCGATAAGCGTCGGTAAGGTGATATGAACCGTTATAACCGGCGATTTCCGAATGGGGAAACCCAGTGTGTTTCGACACACTATCATTAACTGAATCCATAGGTTAATGAGGCGAACCGGGGGAACTGAAACATCTAAGTACCCCGAGGAAAAGAAATCAACCGAGATTCCCCCAGTAGCGGCGAGCGAACGGGGAGCAGCCCAGAGCCTGAATCAGTGTGTGTGTTAGTGGAAGCGTCTGGAAAGGCGTGCGATACAGGGTGACAGCCCCGTACACAAAAATGCACATGCTGTGAGCTCGATGAGTAGGGCGGGACACGTGGTATCCTGTCTGAATATGGGGGGACCATCCTCCAAGGCTAAATACTCCTGACTGACCGATAGTGAACCAGTACCGTGAGGGAAAGGCGAAAAGAACCCCGGCGAGGGGAGTGAAAAAGAACCTGAAACCGTGTACGTACAAGCAGTGGGAGCACGCTTAGGCGTGTGACTGCGTACCTTTTGTATAATGGGTCAGCGACTTATATTCTGTAGCAAGGTTAACCGAATAGGGGAGCCGAAGGGAAACCGAGTCTTAACTGGGCGTTAAGTTGCAGGGTATAGACCCGAAACCCGGTGATCTAGCCATGGGCAGGTTGAAGGTTGGGTAACACTAACTGGAGGACCGAACCGAGTAATGTTGAAAAATTAGCGGATGACTTGTGGCTGGGGGTGAAAGGCCAATCAAACCGGGAGATAGCTGGTTCTCCCCGAAAGCTATTTAGGTAGCGCCTCGTGAATTCATCTCCGGGGGTAGAGCACTGTTTCGGCAAGGGGGTCATCCCGACTTACCAACCCGATGCAAACTGCGAATACCGGAGAATGTTATCACGGGAGACACACGGCGGGTGCTAACGTCCGTCGTGAAGAGGGAAACAACCCAGACCGCCAGCTAAGGTCCCAAAGTCATGGTTAAGTGGGAAACGATGTGGGAAGGCCCAGACAGCCAGGATGTTGGCTTAGAAGCAGCCATCATTTAAAGAAAGCGTAATAGCTCACTGGTCGAGTCGGCCTGCGCGGAAGATGTAACGGGGCTAAACCATGCACCGAAGCTGCGGCAGCGACACTATGTGTTGTTGGGTAGGGGAGCGTTCTGTAAGCCTGTGAAGGTGTGCTGTGAGGCATGCTGGAGGTATCAGAAGTGCGAATGCTGACATAAGTAACGATAAAGCGGGTGAAAAGCCCGCTCGCCGGAAGACCAAGGGTTCCTGTCCAACGTTAATCGGGGCAGGGTGAGTCGACCCCTAAGGCGAGGCCGAAAGGCGTAGTCGATGGGAAACAGGTTAATATTCCTGTACTTGGTGTTACTGCGAAGGGGGGACGGAGAAGGCTATGTTGGCCGGGCGACGGTTGTCCCGGTTTAAGCGTGTAGGCTGGTTTTCCAGGCAAATCCGGAAAATCAAGGCTGAGGCGTGATGACGAGGCACTACGGTGCTGAAGCAACAAATGCCCTGCTTCCAGGAAAAGCCTCTAAGCATCAGGTAACATCAAATCGTACCCCAAACCGACACAGGTGGTCAGGTAGAGAATACCAAGGCGCTTGAGAGAACTCGGGTGAAGGAACTAGGCAAAATGGTGCCGTAACTTCGGGAGAAGGCACGCTGATATGTAGGTGAAGCGACTTGCTCGTGGAGCTGAAATCAGTCGAAGATACCAGCTGGCTGCAACTGTTTATTAAAAACACAGCACTGTGCAAACACGAAAGTGGACGTATACGGTGTGACGCCTGCCCGGTGCCGGAAGGTTAATTGATGGGGTTAGCCGCAAGGCGAAGCTCTTGATCGAAGCCCCGGTAAACGGCGGCCGTAACTATAACGGTCCTAAGGTAGCGAAATTCCTTGTCGGGTAAGTTCCGACCTGCACGAATGGCGTAATGATGGCCAGGCTGTCTCCACCCGAGACTCAGTGAAATTGAACTCGCTGTGAAGATGCAGTGTACCCGCGGCAAGACGGAAAGACCCCGTGAACCTTTACTATAGCTTGACACTGAACATTGAGCCTTGATGTGTAGGATAGGTGGGAGGCTTTGAAGTGTGGACGCCAGTCTGCATGGAGCCGACCTTGAAATACCACCCTTTAATGTTTGATGTTCTAACGTTGACCCGTAATCCGGGTTGCGGACAGTGTCTGGTGGGTAGTTTGACTGGGGCGGTCTCCTCCTAAAGAGTAACGGAGGAGCACGAAGGTTGGCTAATCCTGGTCGGACATCAGGAGGTTAGTGCAATGGCATAAGCCAGCTTGACTGCGAGCGTGACGGCGCGAGCAGGTGCGAAAGCAGGTCATAGTGATCCGGTGGTTCTGAATGGAAGGGCCATCGCTCAACGGATAAAAGGTACTCCGGGGATAACAGGCTGATACCGCCCAAGAGTTCATATCGACGGCGGTGTTTGGCACCTCGATGTCGGCTCATCACATCCTGGGGCTGAAGTAGGTCCCAAGGGTATGGCTGTTCGCCATTTAAAGTGGTACGCGAGCTGGGTTTAGAACGTCGTGAGACAGTTCGGTCCCTATCTGCCGTGGGCGCTGGAGAACTGAGGGGGGCTGCTCCTAGTACGAGAGGACCGGAGTGGACGCATCACTGGTGTTCGGGTTGTCATGCCAATGGCACTGCCCGGTAGCTAAATGCGGAAGAGATAAGTGCTGAAAGCATCTAAGCACGAAACTTGCCCCGAGATGAGTTCTCCCTGACTCCTTGAGAGTCCTGAAGGAACGTTGAAGACGACGACGTTGATAGGCCGGGTGTGTAAGCGCAGCGATGCGTTGAGCTAACCGGTACTAATGAACCGTGAGGCTTAACCTT'

mapRY = {'A' : 'R', 'G' : 'R',
         'C' : 'Y', 'T' : 'Y',
         'N' : 'N'}

def calcMatches(primer, rRNASeq=rRNASeq_correct):
    posList = [0]*(len(rRNASeq))
    for i in range(len(rRNASeq)-len(primer)):
        off = hammingDist(primer, rRNASeq[i:i+len(primer)])
        if off > 20:
            print off
            print i
        posList[i] = len(primer) - off
    return posList

def assignReadsByPos(revReadFile, totalReads, lengthRead, primerList, cogRNASeq=None):
    if (cogRNASeq is None):
        primer = revReadFile.split('_')[-1].split('.')[0]
        cogRNAName = str(primerList[primerList['rtPrimer'].str.slice(0,16).str.upper() == primer]['name'].values[0])
        cogRNAFile = open('/home/jhdavis/data/2014/2014_06_07-ShapeSeqAnalysis/23SSeqs/'+cogRNAName+'_rRNASeq.txt', 'r')
        cogRNAFile.readline()
        cogRNASeq = cogRNAFile.readline()
    else:
        cogRNAName = 'noMatch'
        cogRNASeq = rRNASeq_correct
    fastq_fileR2 = itertools.islice(HTSeq.FastqReader(revReadFile), totalReads)

    posHash = {}
    counts = 0
    while counts < totalReads:
        try:
            r2 = fastq_fileR2.next()
        except StopIteration:
            break
        s = r2.seq
        s = s[:lengthRead]
        try:
            posHash[str(cogRNASeq.index(s[:lengthRead]))] += 1
        except KeyError:
            posHash[str(cogRNASeq.index(s[:lengthRead]))] = 1
        except ValueError:
            try:
                posHash['-1'] +=1
            except KeyError:
                posHash['-1'] = 1
        counts = counts + 1
    
    return [posHash, cogRNAName]


def findIndexCode (fastq_file, numReads):
    codeHash = {}
    for read in itertools.islice(fastq_file, numReads):
        code = ''.join([mapRY[str(b)] for b in read.seq[0:4]])
        try:
            codeHash[code] = codeHash[code] + 1
        except KeyError:
            codeHash[code] = 1
    return codeHash
    
def plotCodeHash(fastq_file, numReads):
    codeHash = findIndexCode(fastq_file, numReads)
    f = pylab.figure(figsize = (10, 5))
    a = f.add_subplot(111)
    a.bar(range(len(codeHash.keys())), [codeHash[i] for i in qMS.sort_nicely(codeHash.keys())], align='center')
    a.set_xticks(range(len(codeHash.keys())))
    a.set_xticklabels(qMS.sort_nicely(codeHash.keys()), rotation=90)
    a.set_xlim(-0.5, len(codeHash.keys()))

def sortReadsUsingCode(read1, read2, readsToUse=1000000000, trunc=25):
    fastq_fileR1 = HTSeq.FastqReader(read1)
    fastq_fileR2 = HTSeq.FastqReader(read2)

    counts = 0
    readPairs = itertools.izip(fastq_fileR1, fastq_fileR2)
    plusReadsR2 = []
    minusReadsR2 = []

    plusReadsR1 = []
    minusReadsR1 = []

    while counts < readsToUse:
        try:
            r1, r2 = readPairs.next()
            r1Trunc = r1[0:trunc]
            r2Trunc = r2[0:trunc]
            
            code = ''.join([mapRY[str(b)] for b in r1.seq[0:4]])
            hdRRRY = hammingDist(code, 'RRRY')
            hdYYYR = hammingDist(code, 'YYYR')
            
            if hdRRRY < 2:
                plusReadsR1.append(r1Trunc.seq)
                plusReadsR2.append(r2Trunc.seq)
            if hdYYYR < 2:
                minusReadsR1.append(r1Trunc.seq)
                minusReadsR2.append(r2Trunc.seq)
            counts = counts + 1 
        except StopIteration:
            break
    return [plusReadsR1, plusReadsR2, minusReadsR1, minusReadsR2]

def assignR2ByPos(p1, p2, m1, m2, rRNASeq=rRNASeq_correct):
    plusSeqsHash = {}
    for s in p2:
        try:
            plusSeqsHash[s] = plusSeqsHash[s]+1
        except KeyError:
            plusSeqsHash[s] = 1

    minusSeqsHash = {}
    for s in m2:
        try:
            minusSeqsHash[s] = minusSeqsHash[s]+1
        except KeyError:
            minusSeqsHash[s] = 1
    
    plusHitsIndexHash = {}
    plusMissesIndexHash = {}
    for s in plusSeqsHash.keys():
        try:
            plusHitsIndexHash[str(rRNASeq.index(s))] += plusSeqsHash[s]
        except ValueError:
            plusMissesIndexHash[s] = plusSeqsHash[s]
        except KeyError:
            plusHitsIndexHash[str(rRNASeq.index(s))] = plusSeqsHash[s]
        
    minusHitsIndexHash = {}
    minusMissesIndexHash = {}
    for s in minusSeqsHash.keys():
        try:
            minusHitsIndexHash[str(rRNASeq.index(s))] += minusSeqsHash[s]
        except ValueError:
            minusMissesIndexHash[s] = minusSeqsHash[s]
        except KeyError:
            minusHitsIndexHash[str(rRNASeq.index(s))] = minusSeqsHash[s]
    
    return [plusHitsIndexHash, plusMissesIndexHash, minusHitsIndexHash, minusMissesIndexHash]


def truncateReads(read1, read2, truncateLength=35, readsToUse = 1000000000, outputPath=None):
    
    fastq_fileR1 = HTSeq.FastqReader(read1)
    fastq_fileR2 = HTSeq.FastqReader(read2)
    if outputPath is None:
        outputPath1 = '.'.join(read1.split('.')[:-1])
        outputPath2 = '.'.join(read2.split('.')[:-1])
        print outputPath1
        print outputPath2
    outFileR1 = open(outputPath1+'_'+'truncated_'+str(truncateLength)+'.fastq', 'w')
    outFileR2 = open(outputPath2+'_'+'truncated_'+str(truncateLength)+'.fastq', 'w')

    readPairs = itertools.izip(fastq_fileR1, fastq_fileR2)
    count = 0
    while count < readsToUse:
        try:
            r1, r2 = readPairs.next()
            r1Trunc = r1[0:truncateLength]
            r2Trunc = r2[0:truncateLength]
            r1Trunc.write_to_fastq_file(outFileR1)
            r2Trunc.write_to_fastq_file(outFileR2)
            count = count+1
        except StopIteration:
            count = readsToUse+1
            break
    outFileR1.close()
    outFileR2.close()
    
def rc(seq):
    complements = maketrans('acgtACGT', 'tgcaTGCA')
    rcseq = seq.translate(complements)[::-1]
    return rcseq

def hammingDist(s1, s2):
    assert len(s1) == len(s2)
    return sum(map(operator.ne, s1, s2))

def buildPrimerFileListHash(RTs, outputPath, read1, read2):
    outFileHash = {}
    
    if outputPath is None:
        outputPath1 = '.'.join(read1.split('.')[:-1])
        outputPath2 = '.'.join(read2.split('.')[:-1])
    
    for i in RTs:
        outFileHash[i] = [open(outputPath1+'_sub_'+i+'.fastq', 'w'),
                            open(outputPath2+'_sub_'+i+'.fastq', 'w'), 0, [0, 0], [0, 0], [0, 0], [0, 0], [0, 0]]
    outFileHash['noMatch'] = [open(outputPath1+'_sub_noMatch.fastq', 'w'),
                              open(outputPath2+'_sub_noMatch.fastq', 'w'), 0, [0, 0], [0, 0], [0, 0], [0, 0], [0, 0]]
    return outFileHash

def assignReadsToPrimers(primerList, read1, read2, maxHamDist=1, rtMatch=13,
                         adapterLength=30, outputPath=None, readsToUse=1000000000, trunc=1000, writeFile=False):
    if not 'rtPrimer' in primerList.columns:
        primerList['rtPrimer'] = primerList['seq'].str.slice(adapterLength+4, stop=None).str.upper()
    RTs = [i.upper()[:rtMatch] for i in list(primerList['rtPrimer'].values)]
    
    
    outFileHash = buildPrimerFileListHash(RTs, outputPath, read1, read2)
    
    fastq_fileR1 = HTSeq.FastqReader(read1)
    fastq_fileR2 = HTSeq.FastqReader(read2)
    readPairs = itertools.izip(fastq_fileR1, fastq_fileR2)
    count = 0
    while count < readsToUse:
        try:
            r1, r2 = readPairs.next()
            r1Trunc = r1[0:trunc]
            r2Trunc = r2[0:trunc]

            found = False
            for primerSeq in RTs:
                hd = hammingDist(primerSeq, r1Trunc.seq[4:rtMatch+4])
                if hd <= maxHamDist:
                    found = True

                    outFileHash[primerSeq][2] += 1
                    code = ''.join([mapRY[str(b)] for b in r1Trunc.seq[0:4]])
                    hdRRRY = hammingDist(code, 'RRRY')
                    hdYYYR = hammingDist(code, 'YYYR')
                    outFileHash[primerSeq][hdRRRY+3][0] += 1
                    outFileHash[primerSeq][hdYYYR+3][1] += 1
                    if writeFile:
                        r1Trunc.write_to_fastq_file(outFileHash[primerSeq][0])
                        r2Trunc.write_to_fastq_file(outFileHash[primerSeq][1])
                    break
            if not found:
                outFileHash['noMatch'][2] += 1
                code = ''.join([mapRY[str(b)] for b in r1Trunc.seq[0:4]])
                hdRRRY = hammingDist(code, 'RRRY')
                hdYYYR = hammingDist(code, 'YYYR')
                outFileHash['noMatch'][hdRRRY+3][0] += 1
                outFileHash['noMatch'][hdYYYR+3][1] += 1

                if writeFile:
                    r1Trunc.write_to_fastq_file(outFileHash['noMatch'][0])
                    r2Trunc.write_to_fastq_file(outFileHash['noMatch'][1])
            count+=1
        except StopIteration:
                break
    for k in outFileHash.keys():
        outFileHash[k][0].close()
        outFileHash[k][1].close()
    return outFileHash

def generateTruncSeqs(rna, primerList, processItems = None, path=None, adapterLength=30, intBarLength=3, printOutput=True):
    if not 'rtPrimer' in primerList.columns:
        primerList['rtPrimer'] = primerList['seq'].str.slice(adapterLength+4, stop=None).str.upper()
    if processItems is None:
        processItems = range(len(primerList))
    for i in processItems:
        primerSeq = rc(primerList.loc[i]['rtPrimer']).upper()
        primerName = primerList.loc[i]['name']
        bindLocation = rna.find(primerSeq[:-intBarLength])
        if bindLocation == -1:
            print "Error in primer "+ primerName + ", cannot find binding location"
        else:
            s = rna[0:bindLocation+len(primerSeq)-intBarLength] + primerSeq[-intBarLength:]
            if printOutput:
                print primerName + ':' + 'OK'
                print s
            if not path==None:
                outFile = open(path+primerName+'_'+rc(primerSeq)+'_rRNASeq.fa', 'w')
                outFile.write('>'+primerName+'\n')
                outFile.write(s)
                outFile.close()
    return primerList

