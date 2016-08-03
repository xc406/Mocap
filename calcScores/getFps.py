"""
calculating footprint scores from DNase-Seq/ATAC-Seq wiggle files
-------------------------------------------------------------------
Usage: python getFps.py wig-file out-path motifChrom ctName cutoff tfName expName
e.g. python getFps.py /path/file.wig /outpath/ chr17 Gm12878 39 CTCF atacseq
"""
import re, os, csv, sys, glob, gzip, time, bz2
import numpy as np
from array import array
from bisect import bisect
from collections import defaultdict
import countWig
from scipy.stats import binom
import pickle, random

def updateFPS(infile, outpath, tfname= "CTCF", motifChrom="chr15", ctName="Gm12878", dgfCutoff=36, expName="atacseq"):#, stranded=False):
        """calculate footprint scores from discontinuous variableStep wiggle files"""
	
	#check if directories exists
	motifdir = os.path.join(outpath,"bedMotifs")
	fpsdir = os.path.join(outpath,"fpsMotifs")
	if not os.path.isdir(motifdir):
		print "Error: path-to-motif-bed-files invalid, please specify a valid outpath to store all calculated scores."
		sys.exit()
		#os.mkdir(motifdir)
	if not os.path.isdir(fpsdir):
		os.mkdir(fpsdir)
	
	print 'updating fps for ', ctName, motifChrom

	wigfilename = re.split(".wig",infile)[0]#os.path.join(inpath,"SRR8912"+"68"+"sort_cut")
	gcoordsfile = gzip.open(os.path.join(motifdir,tfname+motifChrom+".bed.gz"),'r')
	gcoords = csv.reader(gcoordsfile, delimiter='\t')
	fpsfile = gzip.open(os.path.join(fpsdir,tfname+ctName+motifChrom+'fps.txt.gz'),'w')
	writer = csv.writer(fpsfile, delimiter='\t')

	##non-strand-specific fps for atac-seq
	if expName=="atacseq":
		wigFile = open(infile,'rt')
		bwFile = wigfilename+motifChrom+'.bw'
		if not os.path.isfile(bwFile):
	    		countWig.compressVarWig(wigFile, ctName, wigfilename)
		coordDict, valuesDict = countWig.getBinVarCoord(bwFile,ctName)
		arrayDict = defaultdict(list)
		for test in gcoords:
	    		if not motifChrom in arrayDict:
				arrayDict[motifChrom] = countWig.buildVarHist(motifChrom,coordDict,valuesDict,ctName)
	    		xs, xvals, sums = arrayDict[motifChrom]
	    		motifStart, motifEnd = int(test[1]), int(test[2])

	    		flankWin = round((motifEnd - motifStart + 1)*1.75)#35
	    		flankL = max(0, int(motifStart - flankWin))
	    		flankR = int(motifEnd + flankWin)
	    		countTotL = countWig.queryHist(xs, xvals, sums, flankL, motifEnd, varWindow=True)[2]
	    		countTotR = countWig.queryHist(xs, xvals, sums, motifStart, flankR, varWindow=True)[2]
	    		countCent = countWig.queryHist(xs, xvals, sums, motifStart, motifEnd)[2]
	    		count = countWig.queryHist(xs, xvals, sums, motifStart-100, motifEnd+100, varWindow=True)[2]
			count = count-countCent
			if count >= dgfCutoff:
        			acces = 1.0
        			fragP = array("d")
				fragN = array("d")
				for i in xrange(flankL,motifEnd):
                			c = countWig.queryHist(xs, xvals, sums, i, i+1, varWindow=True)[2]
                			fragP.append(c)

        			#Centp = motifStart-flankL
        			countTotL = sum(fragP)#countWig.queryHist(xs, xvals, sums, flankL, motifEnd, varWindow=True)[2]
        			for i in xrange(motifStart,flankR):
                			c = countWig.queryHist(xs, xvals, sums, i, i+1, varWindow=True)[2]
                			fragN.append(c)

        			#Centn = motifEnd+1-motifStart
        			countTotR = sum(fragN)#countWig.queryHist(xs, xvals, sums, motifStart, flankR, varWindow=True)[2]
				try:
                			pp = binom.cdf(countCent,countTotL,float(motifEnd+1-motifStart)/(motifEnd+1-flankL))
                			pn = binom.cdf(countCent,countTotR,float(motifEnd+1-motifStart)/(flankR+1-motifStart))
                			fos = pp*pn
        			except ZeroDivisionError:
                			fos = 1.0

        			##fdr correction
        			fosArray = array("d")
        			for s in xrange(500):
                			random.shuffle(fragP)
					random.shuffle(fragN)
                			try:
                        			pp = binom.cdf(sum(fragP[(motifStart-flankL):]),sum(fragP),float(motifEnd+1-motifStart)/(motifEnd+1-flankL))
                        			pn = binom.cdf(sum(fragN[:(motifEnd+1-motifStart)]),sum(fragN),float(motifEnd+1-motifStart)/(flankR+1-motifStart))
                        			fosArray.append(pp*pn)
                			except ZeroDivisionError:
                        			fosArray.append(1.0)
        			fosCutoff = np.sort(fosArray)[4]
        			#round(sum(1 for s in fosArray if s <= a)/500.0,2) <= 0.01:
        			if fos <= fosCutoff:
                			fps = 1.0##profile
        			else:
                			fps = 0.0##no profile
        			#print tf_name+'\t'+motifChrom+'\t'+str(motifStart)+'\t'+str(motifEnd)+'\t'+str(count)+'\t'+str(fos)+'\t'+str(fdr)
			else:
                        	fos = 1.0+(1/(count+1.0))
                        	fps = -1.0
                        	fosCutoff = -1.0
                        	acces = 0.0
                	row = [count,acces,fos,fosCutoff,fps]
                	writer.writerows([row])

	#run strand-specific calls for DNase-Seq or DGF
    	else:
		shortname = os.path.join(os.path.split(wigfilename)[0],re.split("_",os.path.split(wigfilename)[1])[0])
		infilep = shortname + "_p.wig"#os.path.join(path,"wgEncodeUwDgf"+ctName+"Aln_p_cut.wig")
		infilen = shortname + "_n.wig"#os.path.join(path,"wgEncodeUwDgf"+ctName+"Aln_n_cut.wig")
		wigFilep = open(infilep,'rt')
		wigFilen = open(infilen,'rt')
		bwFilep = shortname+'_p'+motifChrom+'.bw'
		bwFilen = shortname+'_n'+motifChrom+'.bw'
		if not os.path.isfile(bwFilep):
			countWig.compressVarWig(wigFilep,expName,shortname+'_p')
		if not os.path.isfile(bwFilen):
			countWig.compressVarWig(wigFilen,expName,shortname+'_n')
		coordDictp, valueDictp = countWig.getBinVarCoord(bwFilep,ctName)
		coordDictn, valueDictn = countWig.getBinVarCoord(bwFilen,ctName)
		arrayDictp = defaultdict(list)
		arrayDictn = defaultdict(list)

		for test in gcoords:
			if not motifChrom in arrayDictp:
				arrayDictp[motifChrom] = countWig.buildVarHist(motifChrom,coordDictp,valueDictp,ctName)
			if not motifChrom in arrayDictn:
				arrayDictn[motifChrom] = countWig.buildVarHist(motifChrom,coordDictn,valueDictn,ctName)

			motifStart,motifEnd = int(test[1]),int(test[2])
			flankWin = round((motifEnd - motifStart + 1)*2)#1.75)
			flankL= max(0, int(motifStart-flankWin))
			flankR = int(motifEnd + flankWin)

			xsp, xvalsp, sumsp = arrayDictp[motifChrom]
			xsn, xvalsn, sumsn = arrayDictn[motifChrom]

			countCentp = countWig.queryHist(xsp, xvalsp, sumsp, motifStart, motifEnd)[2]
			countCentn = countWig.queryHist(xsn, xvalsn, sumsn, motifStart, motifEnd)[2]
			countp = countWig.queryHist(xsp, xvalsp, sumsp, motifStart-100, motifEnd+100, varWindow=True)[2]
			countn = countWig.queryHist(xsn, xvalsn, sumsn, motifStart-100, motifEnd+100, varWindow=True)[2]
			
			count = countp+countn-countCentp-countCentn
			if count >= dgfCutoff:
				acces = 1.0
				fragP = array("d")
				fragN = array("d")
				for i in xrange(flankL,motifEnd):
					c = countWig.queryHist(xsp, xvalsp, sumsp, i, i+1, varWindow=True)[2]
					fragP.append(c)

				Centp = motifStart-flankL
				countTotLp = sum(fragP)#countWig.queryHist(xs, xvals, sums, flankL, motifEnd, varWindow=True)[2]
				#countCentp = sum(fragP[Centp:])#countWig.queryHist(xs, xvals, sums, motifStart, motifEnd)[2]
				#countp = countWig.queryHist(xs, xvals, sums, motifStart-100, motifEnd+100, varWindow=True)[2]

				for i in xrange(motifStart,flankR):
                                	c = countWig.queryHist(xsn, xvalsn, sumsn, i, i+1, varWindow=True)[2]
                                	fragN.append(c) 

				Centn = motifEnd+1-motifStart
				countTotRn = sum(fragN)#countWig.queryHist(xs, xvals, sums, motifStart, flankR, varWindow=True)[2]
				#countCentn = sum(fragN[:Centn])#countWig.queryHist(xs, xvals, sums, motifStart, motifEnd)[2]
				#countn = countWig.queryHist(xs, xvals, sums, motifStart-100, motifEnd+100, varWindow=True)[2]

				try:
					pp = binom.cdf(countCentp,countTotLp,float(motifEnd+1-motifStart)/(motifEnd+1-flankL))
					pn = binom.cdf(countCentn,countTotRn,float(motifEnd+1-motifStart)/(flankR+1-motifStart))
					fos = pp*pn
				except ZeroDivisionError:
					fos = 1.0

				##fdr correction
				fosArray = array("d")
				for s in xrange(500):
					random.shuffle(fragP)
					random.shuffle(fragN)
					try:
						pp = binom.cdf(sum(fragP[Centp:]),sum(fragP),float(motifEnd+1-motifStart)/(motifEnd+1-flankL))
						pn = binom.cdf(sum(fragN[:Centn]),sum(fragN),float(motifEnd+1-motifStart)/(flankR+1-motifStart))
						fosArray.append(pp*pn)
					except ZeroDivisionError:
						fosArray.append(1.0)
				fosCutoff = np.sort(fosArray)[4]
				#round(sum(1 for s in fosArray if s <= a)/500.0,2) <= 0.01:
				if fos <= fosCutoff:
					fps = 1.0##profile
				else:
					fps = 0.0##no profile
				#print tf_name+'\t'+motifChrom+'\t'+str(motifStart)+'\t'+str(motifEnd)+'\t'+str(count)+'\t'+str(fos)+'\t'+str(fdr)
			else:
				fos = 1.0+(1/(count+1.0))
				fps = -1.0
				fosCutoff = -1.0
				acces = 0.0
			row = [count,acces,fos,fosCutoff,fps]
			writer.writerows([row])
	fpsfile.close()	
        return 0 

def main(argv):
    if len(argv) < 8:
        sys.stderr.write("Usage: %s wig-file out-path motifChrom ctName cutoff tfName expName\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: wig-files %r was not found!\n' % argv[1])
        return 1
    if not os.path.isdir(argv[2]):
        sys.stderr.write('Error: out-path %r was not found!\n' % argv[2])
        return 1

    infile = sys.argv[1]
    outpath = sys.argv[2]
    motifChrom = str(sys.argv[3])#chr12
    ctName = str(sys.argv[4])
    dgfCutoff = float(sys.argv[5])#res.mocha$cutoff from run.Mocha()
    tfname = str(sys.argv[6])
    expname = str(sys.argv[7])#"atacseq" or "dnaseseq"

    #startTime = time.time()

    updateFPS(infile, outpath, tfname=tfname, motifChrom=motifChrom, ctName=ctName, dgfCutoff=dgfCutoff, expName=expname)

    #print 'total time', time.time() - startTime

if __name__=='__main__':
    sys.exit(main(sys.argv))
