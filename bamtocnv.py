#import pysam
import random
import os
import os.path
from os import path
import tempfile
import re
import sys

'''
bam = "/mnt/NL200_1/prod/workspace/IFA20170803002/OncoH/output/170005964BD/normal/5_recal_bam/170005964BD_normal_sort_markdup_realign_recal.bam"
bam_path = "/mnt/X500/farmers/zhengyd/cnvTest_2/OncoH_test/bam/"
'''

def graBam(bam, outb_path, ratio, chrom, start, end):
    sys.stdout.write("0.0000%\r")
    sys.stdout.flush()
    random.seed(1234)
    b1, b2 = splitBam(bam, outb_path)
    for i in range(len(ratio)):
        newb = creatBam(b1, b2, ratio[i], chrom, start, end)
    sys.stdout.write("99.9999%\r")
    sys.stdout.flush()
    newb = re.search(r'[^/]+$', newb).group(0)
    dirsList = []
    dirsList = os.listdir(outb_path)
    for f in dirsList:
        if f != newb:
            filepath = os.path.join(outb_path,f)
            os.remove(filepath)
    sys.stdout.write("100%!finish!\n")

def sort_index(bam):
    d, f = os.path.split(bam); 
    f = f.split('.')[0]
    path_sort = d + "/" + f + '_sort.bam'
    #print path_sort
    pysam.sort('-o', path_sort, bam)
    pysam.index(path_sort)
    return path_sort
    
def splitBam(bam, outb_path):
    
    
    b1 = outb_path + path.split(bam)[1].split('_')[0]+ "_1.bam"
    b2 = outb_path + path.split(bam)[1].split('_')[0]+ "_2.bam"
    
    infile = pysam.AlignmentFile(bam, "rb")

    out1 = pysam.AlignmentFile(b1, "wb", template=infile)
    out2 = pysam.AlignmentFile(b2, "wb", template=infile)
    for s in infile:
        if s.flag != 0x400:
            u = random.random()
            if u < 0.5:
                out1.write(s)
            elif u >= 0.5:
                out2.write(s)
    infile.close()
    out1.close()
    out2.close()
    
    b1 = sort_index(b1)
    b2 = sort_index(b2)
    
    sys.stdout.write("50.0000%\r")
    sys.stdout.flush()
    return b1,b2

def creatBam(b1, b2, ratio, chrom, start, end):
    
    # init
    inf1 = pysam.AlignmentFile(b1, "rb")
    inf2 = pysam.AlignmentFile(b2, "rb")
    
    query1 = set()
    query2 = set()
    
    for i in range(len(chrom)):
        iter = inf1.fetch(str(chrom[i]), start[i], end[i])
        for x in iter:
            query1.add(x.query_name)
        iter = inf2.fetch(str(chrom[i]), start[i], end[i])
        for x in iter:
            query2.add(x.query_name)
    
    inf1.close()
    inf2.close()
    
    inf1 = pysam.AlignmentFile(b1, "rb")
    inf2 = pysam.AlignmentFile(b2, "rb")
    
    
    d, f = os.path.split(b1); 
    f = f.split('.')[0]
    path_out = d + "/" + f + "_%s.bam"% int(ratio*10)
    outf = pysam.AlignmentFile(path_out, "wb", template = inf1)
    
    if ratio < 1:
        for s in inf1:
            u = random.random()
            if (s.query_name in query1) and u > ratio:
                pass
            else:
                outf.write(s)

    else:  
        for s in inf1:
            outf.write(s)
        for s in inf2:
            u = random.random()
            if s.query_name in query2:
                if u <= ratio - 1:
                    outf.write(s)     
                    
    outf.close()
    newb = sort_index(path_out)
    
    inf1.close()
    inf2.close()
    
    return newb
