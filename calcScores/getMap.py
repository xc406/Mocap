"""
calculating mapability scores for motif regions 
-------------------------------------------------
Usage: python getMap.py path-to-map-files map-file-name out-path motifChrom seqlen window 
e.g. python getMap.py /path/ wgEncodeCrgMapabilityAlign100merchr17.bedGraph.bz2 /outpath/ chr17 50 100
"""
import re, os, csv, sys, glob, gzip, time, bz2
import numpy as np
from array import array
from collections import defaultdict
import countWig, countBed
import pickle, random

def updateMap(inpath,infile,outpath,motifChrom,seqlen,window=100):
    """update mapability scores from bedGraph files"""
    motifdir = os.path.join(outpath,"bedMotifs")
    mapdir = os.path.join(outpath,"mapMotifs")
    if not os.path.isdir(motifdir):
        print "Error: path-to-motif-bed-files invalid, please specify a valid outpath to store all calculated scores."
        sys.exit()
    if not os.path.isdir(mapdir):
        os.mkdir(mapdir) 
    mapinfile = os.path.join(inpath, infile)##better in per chromosome bedGraph files
    (mappath,mapfilename) = os.path.split(mapinfile)
    expName = mapfilename.split(motifChrom)[0]

    print 'updating', expName, motifChrom
    with bz2.BZ2File(mapinfile,'r') as bedGraphFile:
        bbFile = os.path.join(inpath,expName+motifChrom+'.bb')
        if not os.path.isfile(bbFile):
            countBed.compressBed4(bedGraphFile, expName, bbFile)
        coordDict, valuesDict = countBed.getBinBedCoord(bbFile, expName)
        arrayDict=defaultdict(list)
        for bedfile in glob.glob(os.path.join(motifdir,"*"+motifChrom+".bed.gz")):
            (filepath,filename) = os.path.split(bedfile)
            tfname = filename.split(motifChrom)[0]
            gcoordsfile = gzip.open(bedfile,'r')
            gcoords = csv.reader(gcoordsfile, delimiter='\t')
            mapfile = gzip.open(os.path.join(mapdir,tfname+motifChrom+'map'+str(seqlen)+'.txt.gz'),'w')
            writer = csv.writer(mapfile, delimiter='\t')
            for test in gcoords:
                motifStart, motifEnd = int(test[1]), int(test[2])
                if not motifChrom in arrayDict:
                    arrayDict[motifChrom] = countBed.buildBedHist(motifChrom, coordDict, valuesDict, expName)
                xs, xvals, sums = arrayDict[motifChrom]
                avg = countBed.queryHist(xs, xvals, sums, motifStart-window, motifEnd+window)[0]
                row = [avg]
                writer.writerows([row])
            mapfile.close()
    return 0

def updateOvlpRegions(inpath,outpath,motifChrom='chr17',window=0, stranded=False):
    """exclude motif sites overlapping ENCODE blacklisted regions"""
    for infile in glob.glob(os.path.join(outpath,"*"+motifChrom+".bed")):
        (filepath,filename) = os.path.split(infile)
        tfname = filename.split(motifChrom)[0]
        gcoordsfile = open(infile,'rt')
        gcoords = csv.reader(gcoordsfile, delimiter='\t')

        for exclinfile in glob.glob(os.path.join(inpath,"*Excludable.bed.gz")):
            #wgEncodeDacMapabilityConsensusExcludable.bed.gz
            (regpath,regfilename) = os.path.split(exclinfile)
            expName = regfilename.split('.')[0]
            #print 'updating', expName
            with gzip.open(exclinfile,'rt') as bedFile:
                bed = csv.reader(bedFile, delimiter='\t')
                annoIntvlDict = countBed.getBed6Anno(bed,expName)
                intervalStartDict = countBed.sortStart(annoIntvlDict)
                intervalEndDict = countBed.sortEnd(annoIntvlDict)
                exclfile = open(os.path.join(outpath+"/"+expName,tfname+motifChrom+'.bed'),'wt')
                writer = csv.writer(exclfile, delimiter='\t')
                if stranded:
                    for test in gcoords:
                        exonList = []
                        motifStart, motifEnd, motifStrand = int(test[1])-window, int(test[2])+window, test[-1]
                        regionList, valueList = countBed.getMotifAnno(annoIntvlDict,
                                intervalStartDict,intervalEndDict,motifChrom,
                                motifStart,motifEnd,window)
                        if valueList != []:##[(exon_name,strand),]
                            for i in xrange(len(valueList)):
                                if valueList[i][1] == motifStrand:
                                    exonList.append(regionList[i])
                        #if exonList != []:
                            #row = ["_".join(exonList)]
                            #writer.writerows([row])
                        if exonList == []:
                            #row = ['ok']
                            writer.writerows([test])
                else:
                    for test in gcoords:
                        motifStart, motifEnd = int(test[1])-window, int(test[2])+window
                        regionList, valueList = countBed.getMotifAnno(annoIntvlDict,
                                intervalStartDict,intervalEndDict,motifChrom,
                                motifStart,motifEnd,window)
                        #if regionList != []:
                            #row = ["_".join(regionList)]
                            #writer.writerows([row])
                        if regionList == []:
                            #row = ['ok']
                            writer.writerows([test])
                exclfile.close()
    return 0

def main(argv):
    if len(argv) < 7:
        sys.stderr.write("Usage: %s path-to-map-files map-file-name out-path motifChrom seqlen window\n" % argv[0])
        return 1
    if not os.path.exists(argv[1]):
        sys.stderr.write('Error: path-to-map-files %r was not found!\n' % argv[1])
        return 1
    if not os.path.isfile(os.path.join(argv[1],argv[2])):
        sys.stderr.write('Error: map-file-name %r was not found!\n' % argv[2])
        return 1
    if not os.path.exists(argv[3]):
        sys.stderr.write('Error: out-path %r was not found!\n' % argv[3])
        return 1
    

    inpath = sys.argv[1]#"/home/xc406/Downloads/hg19_wgEncodeMapability"
    infile = sys.argv[2]#"wgEncodeCrgMapabilityAlign50mer"+motifChrom+".bedGraph.bz2"
    outpath = sys.argv[3]
    motifChrom = str(sys.argv[4])
    seqlen = int(sys.argv[5])#options are 24,36,40,50,75,100
    window = int(sys.argv[6])

    updateMap(inpath=inpath,infile=infile,outpath=outpath,motifChrom=motifChrom,seqlen=seqlen,window=window)

if __name__=='__main__':
    sys.exit(main(sys.argv))
