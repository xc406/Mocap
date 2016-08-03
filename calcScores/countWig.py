"""
functions and scripts to get count data from wig files
for genomic regions defined in gff/bed files
returns a 1-based "bed-4" format count file with columns
chr\tstart\tend\tcount(from[start-100,end+100])
-----------------------------------------------------------
Usage: python countWig.py gff-or-bed-file wig-file out-path ctName
e.g. python countWig.py /outpath/bedMotifs/CTCFchr17.bed.gz file.wig /outpath/ K562
"""
import os, sys, csv, re, gzip
from array import array
from copy import copy
from collections import defaultdict
import time
import struct
import itertools

def getRange(gffFile):
	"""get genomic ranges defined in gff files"""
	##gff_reader shifted the base automatically for 0-based bam output
	featurelist = []
	for feature in gffFile:
		chrom, start, end = feature[0],int(feature[3]),int(feature[4])
		#chrom = feature.iv.chrom
		#start = feature.iv.start + 1
		#end = int(feature.iv.end)
		#strand = feature.iv.strand
		#print start, end
		featurelist.append((chrom,start,end))#,strand))
	return featurelist	

def grouper(n, iterable):
	"""split iterable list by n
		grouper(3, 'ABCDEF') --> ABC DEF"""
    	it = iter(iterable)
	chunklist = []
	chunk = (0,0)
	while chunk:
		if chunk != (0,0) and chunk != ('',):##not storing first or last vals
			chunklist.append(chunk)
		chunk = tuple(itertools.islice(it,n))
	return chunklist

def compressVarWig(wigFile, ctName, bwFileName):
	"""encode wiggle files in binary format as header(cci)+numString(if)"""
	wigFile.seek(0)
	dataString = wigFile.read()
	dataList = dataString.split("variableStep")
	#print 'open for compression', len(dataList)
	for d in dataList[1:]:
		chrom = d.split('\n',1)[0].split('=')[-1]
		with open(bwFileName+chrom+'.bw', 'wb') as f:
			if len(chrom) == 4:
				chrom = '0'+chrom[-1]
			elif len(chrom) == 5:
				chrom = chrom[3:5]
			numString = d.split('\n',1)[1]
			header = struct.Struct('cci')
			ls = list(chrom)
			numLine = len(numString.split('\n'))-1
			ls.append(numLine)
			l = header.pack(*ls)
			f.write(l)
			numst = struct.Struct('if'*numLine)
			lociVal = re.split('\n|\t', numString)
			del lociVal[-1]
			chunklist = grouper(2,lociVal)
			num = []
			for x in chunklist:
				num.append(int(x[0]))
				num.append(float(x[1]))
			l = numst.pack(*num)
			f.write(l)
	return 0

def compressVarWigAll(wigFile, ctName, bwFile):
        """encode wiggle files in binary format as header(cci)+numString(if)"""
        with open(bwFile, 'wb') as f:
                wigFile.seek(0)
                dataString = wigFile.read()
                dataList = dataString.split("variableStep")
                #print 'open for compression', len(dataList)
                for d in dataList[1:]:
                        chrom = d.split('\n',1)[0].split('=')[-1]
                        if len(chrom) == 4:
                                chrom = '0'+chrom[-1]
                        elif len(chrom) == 5:
                                chrom = chrom[3:5]
                        numString = d.split('\n',1)[1]
                        header = struct.Struct('cci')
                        ls = list(chrom)
                        numLine = len(numString.split('\n'))-1
                        ls.append(numLine)
                        l = header.pack(*ls)
                        f.write(l)
                        numst = struct.Struct('if'*numLine)
                        lociVal = re.split('\n|\t', numString)
                        del lociVal[-1]
                        chunklist = grouper(2,lociVal)
                        num = []
                        for x in chunklist:
                                num.append(int(x[0]))
                                num.append(float(x[1]))
                        l = numst.pack(*num)
                        f.write(l)
        return 0

def compressFixWig(wigFile, ctName, bwFile):
	"""fixedStep chrom=chr19 start=3025580 step=1
	encode wiggle files in binary format as header(cciii)+numString(f)"""
	with open(bwFile, 'wb') as f:
		wigFile.seek(0)
		dataString = wigFile.read()
		dataList = dataString.split('fixedStep')
		for d in dataList[1:]:
			chrom = d.split(' ')[1].split('=')[-1]
			start = int(d.split(' ')[2].split("=")[-1])
			step = int(d.split('\n')[0].split(' ')[3].split("=")[-1])
			if len(chrom) == 4:
				chrom = '0'+chrom[-1]
			elif len(chrom) == 5:
				chrom = chrom[3:5]
			numString = d.split('\n',1)[-1]	
			header = struct.Struct('cciii')
			ls = list(chrom)
			numLine = len(numString.split('\n'))-1
			ls.append(start)
			ls.append(step)
			ls.append(numLine)
			l = header.pack(*ls)
			f.write(l)
			numst = struct.Struct('f'*numLine)
			vals = numString.split('\n')
			del vals[-1]
			num = [ float(x) for x in vals ]
			l = numst.pack(*num)
			f.write(l)
	return 0	 

def getBinVarCoord(bwFile, ctName):
	"""get chromosomal coordinates and values stored in dictionaries 
		unpack from bw files """
	coordDict = defaultdict(lambda: defaultdict(list))
	valuesDict = defaultdict(lambda: defaultdict(list))
	numst = struct.Struct('if')
	chunkSize = struct.calcsize('if')
	headerst = struct.Struct('cci')
	headerSize = struct.calcsize('cci')
	with open(bwFile,'rb') as f:
		flag = 1
		header = f.read(headerSize)
		h1,h2,s = headerst.unpack(header)
		if h1 == '0':
			chrom = 'chr'+h2
		else:
			chrom = 'chr'+h1+h2
		while flag:
			for i in xrange(s):
				chunk = f.read(chunkSize)
				data = numst.unpack(chunk)
				coordDict[ctName][chrom].append(int(data[0]))
				valuesDict[ctName][chrom].append(data[1])
			header = f.read(headerSize)
			l = len(header)
			if l < headerSize:
				flag = 0
			else:
				h1, h2, s = headerst.unpack(header)
				if h1 == '0': 
					chrom = 'chr'+h2
				else:
					chrom = 'chr'+h1+h2
	return coordDict, valuesDict

def getBinFixStart(bwFile, dataType):
	"""get chromosomal coordinates and values stored in dictionaries
		unpack from bw files"""
	stepDict = defaultdict(lambda: defaultdict(list))
	startDict = defaultdict(lambda: defaultdict(list))
	valuesDict = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
	numst = struct.Struct('f')
	chunkSize = struct.calcsize('f')
	headerst = struct.Struct('cciii')
	headerSize = struct.calcsize('cciii')
	with open(bwFile,'rb') as f:
		flag = 1
		header = f.read(headerSize)
		h1,h2,start,step,s = headerst.unpack(header)
		if h1 == '0':
			chrom = 'chr'+h2
		else:
			chrom = 'chr'+h1+h2
		startDict[dataType][chrom].append(start)
		stepDict[dataType][chrom].append(step)
		while flag:
			for i in xrange(s):
				chunk = f.read(chunkSize)
				data = numst.unpack(chunk)
				valuesDict[dataType][chrom][start].append(data[0])
			header = f.read(headerSize)
			l = len(header)
			if l < headerSize:
				flag = 0
			else:
				h1,h2,start,step,s = headerst.unpack(header)
				if h1 == '0':
					chrom = 'chr'+h2
				else:
					chrom = 'chr'+h1+h2
			startDict[dataType][chrom].append(start)
			stepDict[dataType][chrom].append(step)
	return stepDict, startDict, valuesDict

def getVarCoord(wigFile, ctName):
        """get chromosomal coordinates and values stored in dictionaries"""
        coordDict = defaultdict(lambda: defaultdict(list))
        valuesDict = defaultdict(lambda: defaultdict(list))
	wigFile.seek(0)
	wigCsv = csv.reader(wigFile, delimiter = '\t')
        for row in wigCsv:
                #print row
                if 'variableStep' in row[0]:
                        chrom = row[0].split('=')[-1]
			#if len(chrom) == 4:
			#	chrom = '0'+chrom[-1]
			#elif len(chrom) == 5:
			#	chrom = chrom[3:5]
                else:
                        coordDict[ctName][chrom].append(int(row[0]))
                        valuesDict[ctName][chrom].append(float(row[1]))
        return coordDict, valuesDict

def getFixStart(wigFile, dataType):
	"""get chromasomal start and values for fixedStepped wig and store in dictionaries"""
	stepDict = defaultdict(lambda: defaultdict(list))
	startDict = defaultdict(lambda: defaultdict(list))
	valuesDict = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
	for row in wigFile:
		#print row
		if "fixedStep" in row[0]:
			chrom = row[0].split(" ")[1].split("=")[-1]
			start = row[0].split(" ")[2].split("=")[-1]
			#print start, chrom
			startDict[dataType][chrom].append(int(start))
			step = row[0].split(" ")[3].split("=")[-1]
			#print step, start, chrom
			stepDict[dataType][chrom].append(int(step))
		else:
			valuesDict[dataType][chrom][int(start)].append(float(row[0]))
	return stepDict, startDict, valuesDict

def buildFixHist(chrom, stepDict, startDict, valuesDict, dataType):
	"""build histogram from wig input with discontinuous values
	e.g.
	fixedStep chrom=chr12 start=60085 step=1
	0.121
	0.102
	0.184
	0.184
	0.155
	0.184
	0.184
	0.184
	-1.194
	0.184
	"""
	arrayDict = defaultdict(list)
	start = startDict[dataType][chrom]
	for s in xrange(len(start)):
		values = valuesDict[dataType][chrom][start[s]]
		#print values
		step =  stepDict[dataType][chrom][s]
		n = len(values)
		x, lastval = start[s], -100000000
		xs, xvals = array("I"), array("d")	# clear and reset histogram, index is unsigned integer, value is double
		xs.append(start[s]-1)
		xvals.append(0)
		for i in xrange(n):		# get the histogram
			x = start[s] + i*step
			val = values[i]
			if val != lastval:
		    		xs.append(x)
		    		xvals.append(val)
			lastval = val
		xs.append(start[s] + n*step)
		xvals.append(0)	# end the histogram explicitly
		#print xs
		nlen = len(xs)
		sums = [0] * nlen
 		sums[0] = xvals[0]*(xs[1]-xs[0])
		for i in xrange(1,nlen-1):
			sums[i] = sums[i-1] + (xs[i+1]-xs[i])*xvals[i]
		arrayDict[start[s]] = [xs, xvals, sums, n]
	
	return arrayDict

def buildVarHist(chrom,coordDict,valuesDict,ctName):#coord,values):
        """build histogram from wig input with discontinuous values
                e.g.
                variableStep chrom=chr13
                19021446        1.00
                19022345        1.00
                19022949        3.00
                19022956        1.00
                19025166        5.00
                19025399        1.00
                19025986        1.00
                19026391        1.00
                19026727        1.00
                """
        coord = [-1] + coordDict[ctName][chrom] # prepend left boundary
        values = [0] + valuesDict[ctName][chrom]
        n = len(values)
        x, lastx, lastval = coord[0], coord[0]-1, -100000000
        xs, xvals = array("i"), array("d")      # clear and reset histogram, index is integer, value is double
        for i in xrange(n):             # get the histogram
                x = coord[i]
                val = values[i]
                if x != lastx+1:
                        xs.append(lastx+1)
                        xvals.append(0)
                        #lastx = lastx+1
                        lastval = 0
                if val != lastval:
                        xs.append(x)
                        xvals.append(val)
                lastx = x
                lastval = val
        xs.append(coord[n-1]+1) # append right boundary
        xvals.append(0)
        nlen = len(xs)
        sums = [0] * nlen
        sums[0] = xvals[0]*(xs[1]-xs[0])
        for i in xrange(1,nlen-1):
                sums[i] = sums[i-1] + (xs[i+1]-xs[i])*xvals[i]
        sums[nlen-1] = sums[nlen-2]
        #print "xs",xs
        #print "xvals",xvals
        #print "sums",sums
        return xs, xvals, sums

def queryHist(xs, xvals, sums, start, end, varWindow=False):
	"""query histogram to get average value of a defined genomic region"""
	if varWindow:
		start = max(min(start,xs[-1]),xs[2])
		end = min(max(end,xs[2]),xs[-1])
	#if start < xs[0]:
	#	print 'start out of range'
	#elif end > xs[-1]:
	#	print 'end out of range'
	n = len(xs) # last ending elemented is appended
	ll, rr = 0, n-1
	while ll <= rr:
		m = (ll+rr)/2
		if xs[m] > start:
			rr = m - 1
		else:
			ll = m + 1
	li = rr#ll	# get the left side index
	ll, rr = 0, n-1
	while ll <= rr:
		m = (ll+rr)/2
		if xs[m] > end:
			rr = m - 1
		else:
			ll = m + 1
	ri = rr # get the right side index
	if ri < li:
		return 0.0	# nothing in between, maybe wrong [start, end]?
	sum = sums[ri-1] - sums[li-1]
	sum -= (start - xs[li]) * xvals[li]		# remove the left extra area
	sum += (end - xs[ri] + 1) * xvals[ri]	# remove the right extra area

	##brute force
	check = False
	if check:
		ans = 0
		for i in range(0,n):
			if start <= xs[i] and (i==n-1 or end >= xs[i+1]):
				if i!=n-1:
					ans += (xs[i+1]-xs[i]) * xvals[i]
			elif xs[i] <= start and (i==n-1 or start < xs[i+1]):
				if i!=n-1:
					ans += (xs[i+1] - start) * xvals[i]
			elif xs[i] <= end and (i==n-1 or end < xs[i+1]):
				ans += (end-xs[i]+1) * xvals[i]
		if ans != sum:
			print "No - Wrong Answer:", ans, sum
			return -1
		else:
			print "Yes", ans

	size = end - start + 1
	avg = sum/size

	return avg, size, sum
	
def main(argv):
	if len(argv) < 5:
		sys.stderr.write("Usage: %s gff-or-bed-file wig-file out-path ctName\n" % argv[0])
		return 1
	if not os.path.isfile(argv[1]):
		sys.stderr.write('Error: gff-or-bed-file %r was not found!\n' % argv[1])
		return 1
	if not os.path.isfile(argv[2]):
		sys.stderr.write('Error: wig-file %r was not found!\n' % argv[2])
		return 1
	if not os.path.isdir(argv[3]):
                sys.stderr.write('Error: out-path %r was not found!\n' % argv[3])
                return 1
	
	# example usage
	wigFile = open(sys.argv[2],'rt')
	(path,fname) = os.path.split(sys.argv[2])
    	(shortname, extension) = os.path.splitext(fname)

	outpath = sys.argv[3]
	bwFile = os.path.join(outpath,shortname+'_compress.bw')
	ctName = sys.argv[4]	

	#compress the input wig file to bw file
	if not os.path.isfile(bwFile):
		compressVarWigAll(wigFile,ctName,bwFile)

	coordDict, valuesDict = getBinVarCoord(bwFile, ctName)

	gff = False
	arrayDict = defaultdict(list)
	if gff:
		tfchr = re.split(".gff",os.path.split(sys.argv[1])[1])[0]
		ofile = gzip.open(os.path.join(outpath,tfchr+ctName+'count.txt.gz'),'w')
		writer = csv.writer(ofile,delimiter="\t")
		gff = open(sys.argv[1],'rt')
		##gffFile = HTSeq.GFF_Reader(gff)
		gffFile = csv.reader(gff, delimiter = '\t')
		features = getRange(gffFile)
		intvlen = len(features)
		##build the arrays (3 needed)
		for i in xrange(intvlen):
			chrom,start,end = features[i][0],features[i][1],features[i][2]
			if not chrom in arrayDict:
				arrayDict[chrom] = buildVarHist(chrom,coordDict,valuesDict,ctName)
			xs, xvals, sums = arrayDict[chrom]
			sum = queryHist(xs, xvals, sums, start, end)[2]	# pass the three arrays to the query
			writer.writerows([[chrom,start,end,sum]])
			#print avg
	else:#bed-4 format[start,end]
		bed = gzip.open(sys.argv[1])
		tfchr = re.split(".bed",os.path.split(sys.argv[1])[1])[0]
		bedFile = csv.reader(bed,delimiter="\t")
		ofile = gzip.open(os.path.join(outpath,tfchr+ctName+'count.txt.gz'),'w')
		writer = csv.writer(ofile,delimiter="\t")
		for row in bedFile:
			chrom,start,end=str(row[0]),int(row[1])-100,int(row[2])+100
			if not chrom in arrayDict:
				arrayDict[chrom] = buildVarHist(chrom,coordDict,valuesDict,ctName)
			xs, xvals, sums = arrayDict[chrom]
			sum = queryHist(xs,xvals,sums,start,end)[2]
			#print avg
			writer.writerows([[chrom,start,end,sum]])
	ofile.close()

if __name__=='__main__':
    sys.exit(main(sys.argv))
