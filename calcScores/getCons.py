"""
calculating conservation scores from phastCons and phyloP conservation tracks in fixed step wiggle format
e.g. chr15.phastCons46way.placental.wigFix.gz 
----------------------------------------------------------------------
Usage: python getCons.py path-to-cons-files out-path motifChrom tfName
e.g. python getCons.py /path/ /outpath/ chr17 CTCF
"""
import re, os, csv, sys, glob, gzip, time, bz2
import numpy as np
from array import array
from collections import defaultdict
import countWig
import pickle, random

def updateCons(inpath, tfname, motifChrom, outpath):
    """update conservation scores in gzipped fixedStep wiggle format"""

    #check if directories exists
    motifdir = os.path.join(outpath,"bedMotifs")
    consdir = os.path.join(outpath,"consMotifs")
    if not os.path.isdir(motifdir):
	print "Error: path-to-motif-bed-files invalid, please specify a valid outpath to store all calculated scores."
        sys.exit()
    	#os.mkdir(motifdir)
    if not os.path.isdir(consdir):
    	os.mkdir(consdir)
    gcoordsfile = gzip.open(os.path.join(motifdir,tfname+motifChrom+".bed.gz"))
    consfile = gzip.open(os.path.join(consdir,tfname+motifChrom+'cons.txt.gz'),'w')
    writer = csv.writer(consfile, delimiter='\t')
    l = []
    consTypes = ["phastCons100way","phastCons46way","phastCons46way.placental","phastCons46way.primates",
                "phyloP100way","phyloP46way","phyloP46way.placental","phyloP46way.primate"]

    for consType in consTypes:
        infile = os.path.join(inpath, motifChrom+"."+consType+".wigFix.gz")
        (wigpath,wigfilename) = os.path.split(infile)
        chrom = wigfilename.split('.')[0]
        consName = '_'.join(wigfilename.split('.')[1:-2])

        print 'updating', consName

        gcoordsfile.seek(0)
        gcoords = csv.reader(gcoordsfile, delimiter='\t')

        with gzip.open(infile) as wigFile:
            bwFile = os.path.join(wigpath,motifChrom+"."+consName+'.bw')
            if not os.path.isfile(bwFile):
                countWig.compressFixWig(wigFile, consName, bwFile)
            stepDict, startDict, valuesDict = countWig.getBinFixStart(bwFile,consName)
            start = startDict[consName][chrom]
            arrayDict = countWig.buildFixHist(chrom,stepDict,startDict,valuesDict,consName)
            r = []
	
	    for test in gcoords:
                motifStart, motifEnd = int(test[1]), int(test[2])
                #print motifStart, motifEnd
                avg = 0
                startlist = [start[i] for i in xrange(len(start)-1) if (motifStart >= start[i] and motifStart < start[i+1]) or (motifEnd >= start[i] and motifEnd < start[i+1])]
                if motifEnd > start[-1]:
                    startlist.append(start[-1])
                for i in xrange(len(startlist)):
                    #if avg != 0:
                        #if motifEnd >= startlist[i]:##cases of partial overlap need to renormalize over two fragments
                    ss = startlist[i]
                    xs, xvals, sums, ll = arrayDict[ss]
                    if motifStart < ss <= motifEnd <= ss+ll-1:##left out, right in
                        if avg == 'NA' and i == len(startlist)-1:
                            avg = 0
                        avg += countWig.queryHist(xs,xvals, sums, ss, motifEnd)[0] *(motifEnd - ss + 1) /(motifEnd - motifStart + 1)
                    elif ss <= motifStart < motifEnd <= ss+ll-1:##in array
                        avg = countWig.queryHist(xs,xvals, sums, motifStart, motifEnd)[0]
                    elif motifStart < ss and ss+ll-1 < motifEnd:##motif > array
                        if avg == 'NA':
                            avg = 0
                        avg += countWig.queryHist(xs,xvals, sums, ss, ss+ll-1)[0] * ll /(motifEnd - motifStart + 1)
                    elif ss <= motifStart <= ss+ll-1 < motifEnd:##left in, right out
                        if avg == 'NA' and i == len(startlist)-1:
                            avg = 0
                        avg += countWig.queryHist(xs,xvals, sums, motifStart, ss+ll-1)[0] *(ss + ll - motifStart) /(motifEnd - motifStart + 1)
                    elif ss+ll-1 < motifStart:
                        if avg == 0:
                            #print '...', motifStart, motifEnd, ss, ll
                            avg = 'NA'
                    elif motifEnd < ss:
                        print "Error: motifStart < motifEnd < ss "
                        if avg == 0:
                            avg = 'NA'

                r.append(avg)
            l.append(tuple(r))
    wl = zip(*l)

    for i in wl:
        writer.writerows([list(i)])

    consfile.close()

    return 0

def main(argv):
    if len(argv) < 4:
        sys.stderr.write("Usage: %s path-to-cons-files out-path motifChrom tfName\n" % argv[0])
        return 1
    if not os.path.exists(argv[1]):
        sys.stderr.write('Error: path-to-cons-files %r was not found!\n' % argv[1])
        return 1
    if not os.path.exists(argv[2]):
        sys.stderr.write('Error: out-path %r was not found!\n' % argv[2])
        return 1

    inpath = sys.argv[1]
    motifChrom = str(sys.argv[3])
    tfName = str(sys.argv[4])
    outpath = sys.argv[2]

    updateCons(inpath,tfName,motifChrom,outpath)

if __name__=='__main__':
    sys.exit(main(sys.argv))
