"""
calculating sequence features for motif regions
from fasta files e.g. 'hg19_chrom%02d' % chrom_num
---------------------------------------------------
Usage: python getGC.py path-to-fasta-files out-path motifChrom
e.g. python getGC.py /path/ /outpath/ chr17
"""
import re, os, csv, sys, glob, gzip, time, bz2
import numpy as np
from array import array
from bisect import bisect
from collections import defaultdict
import countWig, countBed
from Bio import SeqIO
from scipy.stats import binom
import pickle, random

def updateGC(inpath, motifChrom, outpath):
    """update GC content, CpG count and repeat sequence for all genomic regions
        e.g. updateGC(path-to-fasta-files chr17 chrom17 out-path)"""
    motifdir = os.path.join(outpath,"bedMotifs")
    gcdir = os.path.join(outpath,"gcMotifs")
    if not os.path.isdir(motifdir):
        print "Error: path-to-motif-bed-files invalid, please specify a valid outpath to store all calculated scores."
        sys.exit()
    if not os.path.isdir(gcdir):
        os.mkdir(gcdir)

    print 'updating GC/CgG for', motifChrom
    winl=[0,100,400]
    for infile in glob.glob(os.path.join(motifdir,"*"+motifChrom+".bed.gz")):
        (filepath,filename) = os.path.split(infile)
        tfname = filename.split(motifChrom)[0]
        gcoordsfile = gzip.open(infile)
        gcoords = csv.reader(gcoordsfile, delimiter='\t')
        gcfile = gzip.open(os.path.join(gcdir,tfname+motifChrom+'gc.txt.gz'),'w')
        writer = csv.writer(gcfile, delimiter='\t')

	if not motifChrom[3:]=="X":
        	chrnum = int(motifChrom[3:])
	else:
		chrnum = 23
        fname = 'hg19_chrom%02d' % chrnum

        with open(os.path.join(inpath,fname)) as f:
            seq_record = SeqIO.read(f,"fasta")
            sl = len(seq_record.seq)
            for test in gcoords:
		l = []
                for win in winl:
                    motifStart, motifEnd = max(int(test[1])-win,1), min(int(test[2])+win,sl)
                    seq = seq_record.seq[motifStart-1:motifEnd]##(]
                    size = float(motifEnd-motifStart+1)
		    cpg = float(seq.count('CG') + seq.count('Cg') + seq.count('cG'))
                    c = float(seq.count('c') + seq.count('C'))
                    g = float(seq.count('g') + seq.count('G'))
                    try:
                    	r = (cpg * size)/((c+g)/2)**2
                    except ZeroDivisionError:
                        r = -1
                    gc = (c+g)/size
                    #print "nl", gc, "\t", motifStart,"\t", motifEnd, "\t", test[1], "\t", test[2], "\t", str(seq)
                    if size == float(seq.count('N')):
                        #print 'N repeat'#, seq
                        row = ['NA', 'NA', 'NA', 5]
                        #writer.writerows([row])
                        l.extend(row)
                    elif size == float(seq.count('c')+seq.count('g')+seq.count('a')+seq.count('t')):
                        #print 'repeat found'#, gc, seq
                        row = [cpg, r, gc, 3]#, str(seq)]
                        #writer.writerows([row])
                        l.extend(row)
                    elif size > float(seq.count('N')) and float(seq.count('N')) > 0:
                        #print 'partial N repeat'#, seq
                        row = ['NA', 'NA', 'NA', 4]#, str(seq)]
                        #writer.writerows([row])
                        l.extend(row)
                    elif float(seq.count('N')) == 0 and size > float(seq.count('c')+seq.count('g')+seq.count('a')+seq.count('t')) and size > float(seq.count('C')+seq.count('G')+seq.count('A')+seq.count('T')):
                        #print 'partial repeat'#, gc, seq
                        row = [cpg, r, gc, 2]#, str(seq)]
                        #writer.writerows([row])
                        l.extend(row)
                    elif size == float(seq.count('C')+seq.count('G')+seq.count('A')+seq.count('T')):
                        #print 'non-repeat', gc, seq
                        row = [cpg, r, gc, 1]#, str(seq)]
                        #writer.writerows([row])
                        l.extend(row)
                writer.writerows([l])

    	    gcfile.close()
    return 0 

def main(argv):
    if len(argv) < 4:
        sys.stderr.write("Usage: %s path-to-fasta-files out-path motifChrom\n" % argv[0])
        return 1
    if not os.path.exists(argv[1]):
        sys.stderr.write('Error: path-to-fasta-files %r was not found!\n' % argv[1])
        return 1
    if not os.path.exists(argv[2]):
        sys.stderr.write('Error: out-path %r was not found!\n' % argv[2])
        return 1

    #chromNum = str(sys.argv[4])
    #if chromNum[-2] == '0':
    #    motifChrom = 'chr'+chromNum[-1]
    #elif not chromNum[5:] == '23':
    #    motifChrom = 'chr'+chromNum[5:]
    #else:
    #    motifChrom = 'chrX'

    inpath = sys.argv[1]#"/home/xc406/data/hg19_chromFa"#hg19_chrom15
    outpath = sys.argv[2]
    motifChrom = sys.argv[3]

    updateGC(inpath=inpath,motifChrom=motifChrom,outpath=outpath)
    #updateCpG(inpath=inpath,bedpath=bedpath,motifChrom=motifChrom,chromNum=chromNum,outpath=outpath)

if __name__=='__main__':
    sys.exit(main(sys.argv))
