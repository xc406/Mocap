"""this script parses atac-seq/dhs bam alignment files 
to extract 5' cut the transposon insertion/DNase I cut sites
requires pysam version 0.7~0.8.4
-------------------------------------------------------------
Usage: python extractCutSite.py bam_file chrom_size output_dir stranded paired
e.g. python extractCutSite.py example.bam human.hg19.chrom.size /out_dir/ 1 paired
"""
import sys
import os
import csv
import pysam
from os import getenv
from collections import defaultdict
import subprocess

def cutBam2(bamfile,hi,cl,stranded=False):
	"""cut pair-end sequencing bam files"""
	dupCounter = defaultdict(int)
	rlen, cap = 100, 50
	if stranded:
                ofilename = (shortname+'_p',shortname+'_n')
                ofilep = pysam.Samfile(os.path.join(temppath,ofilename[0]+'.bam'),'wb',
                        referencenames=[c for c, l in hi],
                        referencelengths=[int(l) for c, l in hi])
                ofilen = pysam.Samfile(os.path.join(temppath,ofilename[1]+'.bam'),'wb',
                        referencenames=[c for c, l in hi],
                        referencelengths=[int(l) for c, l in hi])
	else:
                ofilename = (shortname+'_b',)
                ofile = pysam.Samfile(os.path.join(temppath,ofilename[0]+'.bam'),'wb',
                        referencenames=[c for c, l in hi],
                        referencelengths=[int(l) for c, l in hi])
	for s in bamfile:
		a = pysam.AlignedRead()
		if s.mapq >= 30 and (s.rname in cl) and s.is_duplicate == False:
		##discard low q reads and reads mapped to chrY chrM chrRandom
			if rlen > s.rlen:
                                rlen = s.rlen
                        a.rname = cl.index(s.rname)
			a.mapq = s.mapq
			a.qname = s.qname
			a.flag = s.flag
			a.tags = s.tags
			a.rnext = s.rnext
			a.is_paired = s.is_paired
			a.cigar = [(0,1)]

			if a.flag == 83 or a.flag == 147:##minus strand -5
				a.pos = int(s.pos)+int(s.rlen)-5
				a.seq = s.seq[s.rlen-5]
				a.qual = s.qual[s.rlen-5]
				a.pnext = int(s.pnext)+4
				a.tlen = int(a.pnext) - int(a.pos)
				if a.tlen > 0:
					a.tlen = -a.tlen
				dupCounter[(a.rname,a.pos,a.flag)] +=1
				if dupCounter[(a.rname,a.pos,a.flag)] <= cap:
					if stranded:
						ofilen.write(a)
					else:
						ofile.write(a)
					if dupCounter[(a.rname,a.pos,a.flag)] == cap:
						a.is_duplicate ==True

			elif a.flag == 163 or a.flag == 99:##plus strand +4
				a.pos = int(s.pos)+4
				a.seq = s.seq[4]
				a.qual = s.qual[4]
				a.pnext = int(s.pnext)+int(s.rlen)-5
				a.tlen = int(a.pnext) - int(a.pos)
				if a.tlen < 0:
					a.tlen = -a.tlen
				dupCounter[(a.rname,a.pos,a.flag)] +=1
				if dupCounter[(a.rname,a.pos,a.flag)] <= cap:
					if dupCounter[(a.rname,a.pos,a.flag)] == cap:	
						a.is_duplicate ==True
					if stranded:
						ofilep.write(a)
					else:
						ofile.write(a)
	bamfile.close()
	print "read length: ", rlen
	if stranded:
		ofilep.close()
		ofilen.close()
	else:
		ofile.close()
	return ofilename

def cutBam1(bamfile,hi,cl,stranded=False):
	"""cut single-end sequencing bam files"""
        dupCounter = defaultdict(int)
	rlen, cap = 100, 50
	if stranded:
		ofilename = (shortname+'_p',shortname+'_n')
                ofilep = pysam.Samfile(os.path.join(temppath,ofilename[0]+'.bam'),'wb', 
			referencenames=[c for c, l in hi],
			referencelengths=[int(l) for c, l in hi])
                ofilen = pysam.Samfile(os.path.join(temppath,ofilename[1]+'.bam'),'wb', 
			referencenames=[c for c, l in hi],
			referencelengths=[int(l) for c, l in hi])
        else:
		ofilename = (shortname+"_b",)

		ofile = pysam.Samfile(os.path.join(temppath,ofilename[0]+'.bam'),'wb', 
			referencenames=[c for c, l in hi], 
			referencelengths=[int(l) for c, l in hi])
        for s in bamfile:

		if s.mapq >= 25 and (s.rname in cl) and s.is_duplicate == False:
                ##discard low q reads and reads mapped to chrY chrM chrRandom
			a = pysam.AlignedRead()
			if rlen > s.rlen:
				rlen = s.rlen
			a.rname = cl.index(s.rname)
                        a.mapq = s.mapq
                        a.qname = s.qname
                        a.flag = s.flag
                        a.tags = s.tags
                        a.rnext = s.rnext
                        a.is_paired = s.is_paired
                        a.cigar = [(0,1)]#s.cigar
                        if a.flag == 16:##minus strand -5
                                a.pos = int(s.pos)+int(s.rlen)-1
                                a.seq = s.seq[s.rlen-1]
                                a.qual = s.qual[s.rlen-1]
                                dupCounter[(a.rname,a.pos,a.flag)] +=1
                                if dupCounter[(a.rname,a.pos,a.flag)] <= cap:
					if stranded:
                                        	ofilen.write(a)
					else:
						ofile.write(a)
                                	if dupCounter[(a.rname,a.pos,a.flag)] == cap:
                                        	a.is_duplicate ==True

                        elif a.flag == 0:##plus strand +4
                                a.pos = int(s.pos)
                                a.seq = s.seq[0]
                                a.qual = s.qual[0]
                                dupCounter[(a.rname,a.pos,a.flag)] +=1
                                if dupCounter[(a.rname,a.pos,a.flag)] <= 50:
					if stranded:
                                        	ofilep.write(a)
					else:
						ofile.write(a)
                                	if dupCounter[(a.rname,a.pos,a.flag)] == 50:
                                        	a.is_duplicate ==True

        bamfile.close()
	print "read length: ", rlen
	if stranded:
        	ofilep.close()
		ofilen.close()
	else:
		ofile.close()
	return ofilename

def countDuplicate2(bamfile):
    """deprecated"""
    duplicateCounter = defaultdict(int)
    for s in bamfile:
	if s.mapq >= 30 and s.rname < 20 and s.is_duplicate == False:
		if s.flag == 83 or s.flag == 147:
			pos = int(s.pos) + int(s.rlen) -5
        		duplicateCounter[(s.rname,pos)] += 1
			#if duplicateCounter[(s.rname,pos)] < 31:
				#ofile.write(s)
		elif s.flag == 163 or s.flag == 99:
			pos = int(s.pos)+4
			duplicateCounter[(s.rname,pos)] += 1
			
    return duplicateCounter

def countDuplicate1(bamfile):
    """deprecated"""
    duplicateCounter = defaultdict(int)
    for s in bamfile:
	if s.mapq >= 20 and s.rname < 20 and s.is_duplicate == False:
		if s.flag == 16:
			pos = int(s.pos) + int(s.rlen) -5
			duplicateCounter[(s.rname,pos)] += 1
		elif s.flag == 0:
			pos = int(s.pos)+4
			duplicateCounter[(s.rname,pos)] += 1

    return duplicateCounter

def shiftBam(bamfile,ofile):
	for s in bamfile:
		a = pysam.AlignedRead()
		a.qname = s.qname
		a.flag = s.flag
		a.rname = s.rname
		if s.mapq > 30:
			a.mapq = s.mapq
			a.tags = s.tags
			a.mrnm = s.mrnm
			a.mpos = s.mpos
			a.cigar = s.cigar
			if a.flag == 83 or a.flag == 147:
				a.pos = int(s.pos) + int(s.rlen) -5
				a.seq = s.seq
				a.qual = s.qual
				a.pnext = int(s.pnext) + 4
			elif a.flag == 163 or a.flag == 99:
				a.pos = int(s.pos) + 4
				a.seq = s.seq
				a.qual = s.qual
				a.pnext = int(s.pnext) + int(s.rlen) - 5
			a.tlen = int(a.pnext) - int(a.pos)
			ofile.write(a)
	ofile.close()

def main(argv):
    if len(argv) < 6:
        sys.stderr.write("Usage: %s bam_file chrom_size output_dir stranded paired\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: bam_file %r was not found!\n' % argv[1])
        return 1
    if not os.path.isfile(argv[2]):
        sys.stderr.write('Error: chrom_size %r was not found!\n' % argv[2])
        return 1
    if not os.path.exists(argv[3]):
        sys.stderr.write('Error: output_dir %r was not found!\n' % argv[3])
        return 1

    global temppath, opath, shortname

    ##check pysam version 0.7~0.8.4 required
    if pysam.__version__>='0.9.0' or pysam.__version__<'0.7':
                print "pysam version 0.7 or 0.8 required. run pip install pysam==0.8.4"
                sys.exit()

    opath = sys.argv[3]

    if getenv('MYTMP'):
        temppath = getenv('MYTMP')
    else:
	temppath = opath

    infile = sys.argv[1]

    (path,fname) = os.path.split(infile)
    (shortname, extension) = os.path.splitext(fname)

    chromsize = sys.argv[2]
    hi = [l[:-1].split() for l in open(chromsize,"rt")]
    referencenames=[c for c, l in hi]    
    bamfile = pysam.Samfile(infile,'rb')
    bh = bamfile.header['SQ']
    cl = [x for x,y in enumerate(bh) if y['SN'] in referencenames]

    stranded = bool(int(sys.argv[4]))
    paired = sys.argv[5]

    if paired=="paired":    
    	ofilename = cutBam2(bamfile, hi, cl, stranded)
    else:
	ofilename = cutBam1(bamfile, hi, cl, stranded)

    for i in xrange(len(ofilename)):
	#if pysam.__version__>='0.9.0' or pysam.__version__<'0.7':
		#print "pysam version 0.8 required. run pip install pysam==0.8.4"
		#break
		#pysam.sort("-o",os.path.join(opath, ofilename[i]+'_cut'),os.path.join(temppath, ofilename[i]+'.bam'))
	if '0.9.0'>pysam.__version__>='0.7.0':
    		pysam.sort(os.path.join(temppath, ofilename[i]+'.bam'), os.path.join(opath, ofilename[i]+'_cut'))
    		pysam.index(os.path.join(opath, ofilename[i]+'_cut.bam'))
 

if __name__=='__main__':
    sys.exit(main(sys.argv))



