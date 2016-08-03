"""
functions to get the closest target genes and distance to closest TSS/TES for motif sites
---------------------------------------------------------
Usage: python getDtss.py refseq-file motifChrom out-path
e.g. python getDtss.py /path/refseq_file chr17 /outpath/
"""
import re, os, csv, sys, glob, gzip, time, bz2
import numpy as np
from array import array
from bisect import bisect
from collections import defaultdict
import countWig
from scipy.stats import binom
import pickle, random

def getRefSeqDict(refSeqReader):
    """store transcripts loci info in dictionaries"""
    tssDict = defaultdict(list)
    tesDict = defaultdict(list)
    geneNameDict = defaultdict(tuple)
    sizeDict = defaultdict(int)
    geneRangeDict = defaultdict(lambda : defaultdict(tuple))
    for row in refSeqReader:
        try:
            chrom, strand, txStart, txEnd, geneName, transcriptId = row[2], row[3], int(row[4]), int(row[5]), row[12], row[1]
            size = txEnd - txStart
            if strand == '+':
                if not geneName in sizeDict:
                    tssDict[chrom].append(txStart)
                    tesDict[chrom].append(txEnd)
                    geneNameDict[(chrom,txStart)] = (geneName, transcriptId)
                    sizeDict[geneName] = size
                    #geneRangeDict[chrom][Interval(txStart,txEnd)] = (geneName, transcriptId)
                    geneRangeDict[chrom][(txStart,txEnd)] = (geneName, transcriptId)
                else:
                    if size > sizeDict[geneName]:
                        tssDict[chrom].append(txStart)
                        tesDict[chrom].append(txEnd)
                        geneNameDict[(chrom,txStart)] = (geneName, transcriptId)
                        sizeDict[geneName] = size
                        geneRangeDict[chrom][(txStart,txEnd)] = (geneName, transcriptId)
            else:
                if not geneName in sizeDict:
                    tssDict[chrom].append(txEnd)
                    tesDict[chrom].append(txStart)
                    geneNameDict[(chrom,txEnd)] = (geneName, transcriptId)
                    sizeDict[geneName] = size
                    geneRangeDict[chrom][(txStart,txEnd)] = (geneName, transcriptId)
                else:
                    if size > sizeDict[geneName]:
                        tssDict[chrom].append(txEnd)
                        tesDict[chrom].append(txStart)
                        geneNameDict[(chrom,txEnd)] = (geneName, transcriptId)
                        sizeDict[geneName] = size
                        geneRangeDict[chrom][(txStart,txEnd)] = (geneName, transcriptId)
        except ValueError:
            pass
    return tssDict, tesDict, geneNameDict, geneRangeDict

def getClosestGene(tssDict,geneNameDict,motifChrom,motifStart,motifEnd):
    """find closest gene feature to motif entries"""
    for chrom in tssDict:
        tssDict[chrom].sort()##don't sort beforehand, change the original dictionary
    motifMid = motifStart+(motifEnd-motifStart)/2
    try:
        i = bisect(tssDict[motifChrom],motifMid)
        if i < len(tssDict[motifChrom]) and i > 0:
            if abs(motifMid - tssDict[motifChrom][i-1]) >= abs(tssDict[motifChrom][i] - motifMid):
                (geneName, transcriptId) = geneNameDict[(motifChrom,tssDict[motifChrom][i])]
                dist = abs(tssDict[motifChrom][i] - motifMid)
            else:
                (geneName, transcriptId) = geneNameDict[(motifChrom,tssDict[motifChrom][i-1])]
                dist = abs(motifMid - tssDict[motifChrom][i-1])
        elif i == len(tssDict[motifChrom]):
            if not i == 0:
                (geneName, transcriptId) = geneNameDict[(motifChrom,tssDict[motifChrom][i-1])]
                dist = abs(motifMid - tssDict[motifChrom][i-1])
            else:
                (geneName, transcriptId) = (None, None)
                dist = None
        elif i == 0:
            (geneName, transcriptId) = geneNameDict[(motifChrom,tssDict[motifChrom][i])]
            dist = abs(motifMid - tssDict[motifChrom][i])
    except KeyError:
        geneName, dist, transcriptId = None, None, None
    return geneName, dist, transcriptId

def getDist(tsDict,motifChrom,motifStart,motifEnd):
    """find distance to tss/tes for motif entries"""
    for chrom in tsDict:
        tsDict[chrom].sort()##don't sort beforehand, change the original dictionary
    motifMid = motifStart+(motifEnd-motifStart)/2
    try:
        i = bisect(tsDict[motifChrom],motifMid)
        if i < len(tsDict[motifChrom]) and i > 0:
            if abs(motifMid - tsDict[motifChrom][i-1]) >= abs(tsDict[motifChrom][i] - motifMid):
                dist = abs(tsDict[motifChrom][i] - motifMid)
            else:
                dist = abs(motifMid - tsDict[motifChrom][i-1])
        elif i == len(tsDict[motifChrom]):
            if not i == 0:
                dist = abs(motifMid - tsDict[motifChrom][i-1])
            else:
                dist = None
        elif i == 0:
            dist = abs(motifMid - tsDict[motifChrom][i])
    except KeyError:
        dist  = None
    return dist

def getTargetGene(geneRangeDict, intervalStartDict, intervalEndDict, motifChrom, motifStart, motifEnd, txWindow):
    """push motif target genes into a list"""
    targetList = []
    transcriptList = []
    #itl = IntervalList( intv for intv in geneRangeDict[motifChrom])
    #itl.sort()
    try:
        startList = [start-1-txWindow for (start,end) in intervalStartDict[motifChrom]]
        ##-1 changes 1-based encoding into 0-based
        endList = [end+txWindow for (start,end) in intervalEndDict[motifChrom]]
        intervalStartList = [(start,end) for (start,end) in intervalStartDict[motifChrom]]
        intervalEndList = [(start,end) for (start,end) in intervalEndDict[motifChrom]]
        iStart = bisect(startList, motifEnd)
        iEnd = bisect(endList, motifStart)
        s = intervalStartList[:iStart]##bed intervals are 0-based, so overlapping with bed start is not accounted for
        e = intervalEndList[iEnd:]##motifStart is 1-based, so overlapping with bed end is acccounted for
        overlapping = list(set(s)&set(e))
        #motifInterval = Interval(motifStart, motifEnd)
        #overlapping = [ x for x in intervals if (x[1]>motifInterval[0] and x[0]<motifInterval[0]) \
        #or (x[1]>motifInterval[1] and x[0]<motifInterval[1]) ]
        #print 'overlapping',overlapping

        if not len(overlapping) == 0:
            for x in overlapping:
                #print x
                (geneName, transcriptId) = geneRangeDict[motifChrom][x]
                if not geneName in targetList:
                    targetList.append(geneName)
                if not (geneName,transcriptId) in transcriptList:
                    transcriptList.append((geneName,transcriptId))
        #if len(targetList) > 1:
            #print targetList, motifChrom, motifInterval
    except KeyError:
        pass
    return targetList, transcriptList

def main(argv):
    if len(argv) < 4:
        sys.stderr.write("Usage: %s refseq-file motifChrom out-path\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: refseq-files %r was not found!\n' % argv[1])
        return 1
    if not os.path.isdir(argv[3]):
        sys.stderr.write('Error: out-path %r was not found!\n' % argv[3])
        return 1

    infile = sys.argv[1]
    motifChrom = str(sys.argv[2])
    outpath = sys.argv[3]

    refSeqFile = open(infile,'rt')
    refSeqReader = csv.reader(refSeqFile, delimiter='\t')
    tssDict, tesDict, geneNameDict, geneRangeDict = getRefSeqDict(refSeqReader)

    #intervalStartDict = countBed.sortStart(geneRangeDict)
    #intervalEndDict = countBed.sortEnd(geneRangeDict)

    motifdir = os.path.join(outpath,"bedMotifs")
    targetdir = os.path.join(outpath,"targetMotifs")
    if not os.path.isdir(motifdir):
        print "Error: path-to-motif-bed-files invalid, please specify a valid outpath to store all calculated scores."
        sys.exit()
    if not os.path.isdir(targetdir):
        os.mkdir(targetdir)

    for infile in glob.glob(os.path.join(motifdir,"*"+motifChrom+".bed.gz")):
        (filepath,filename) = os.path.split(infile)
        tfname = filename.split(motifChrom)[0]
        gcoordsfile = gzip.open(infile)
        gcoords = csv.reader(gcoordsfile, delimiter='\t')
        targetfile = gzip.open(os.path.join(targetdir,tfname+motifChrom+'target.txt.gz'),'w')
        writer = csv.writer(targetfile, delimiter='\t')
	for test in gcoords:
		motifStart, motifEnd = int(test[1]), int(test[2])
		#tssdist = getDist(tssDict,motifChrom,motifStart,motifEnd)
		closest = getClosestGene(tssDict,geneNameDict,motifChrom,motifStart,motifEnd)
		row = ["__".join([closest[0],closest[2]]),closest[1]]#target_gene, dist
		writer.writerows([row])
	targetfile.close()

if __name__=='__main__':
    sys.exit(main(sys.argv))
