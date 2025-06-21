#-*- coding: UTF-8 -*
import numpy as np
import argparse
import re
import os

np.seterr(invalid='ignore')

def converse(x,mu,sigma):
    return (1./np.sqrt(2*np.pi))*(np.exp(-(x-mu)**2/(2*sigma**2)))

def em(dataArray,k,mu,sigma,step = 10):
    n = len(k)
    dataNum = dataArray.size
    gamaArray = np.zeros((n, dataNum))
    for s in range(step):
        for i in range(n):
            for j in range(dataNum):
                Sum = sum([k[t]*converse(dataArray[j], mu[t], sigma[t]) for t in range(n)])
                gamaArray[i][j] = k[i]*converse(dataArray[j], mu[i], sigma[i])/float(Sum)
        for i in range(n):
            mu[i] = np.sum(gamaArray[i]*dataArray)/np.sum(gamaArray[i])
        for i in range(n):
            sigma[i] = np.sqrt(np.sum(gamaArray[i]*(dataArray - mu[i])**2)/np.sum(gamaArray[i]))
        for i in range(n):
            k[i] = np.sum(gamaArray[i])/dataNum
    return [k, mu, sigma]

parser = argparse.ArgumentParser()
parser.add_argument('-bam', '--bam', help='the input_bam_file_path', required=True, type=str)
parser.add_argument('-chr', '--chr', help='the chr', required=True, type=str)
parser.add_argument('-begin', '--begin', help='the begin', required=True, type=int)
parser.add_argument('-end', '--end', help='the end', required=True, type=int)
parser.add_argument('-out', '--outpath', help='the out_contiguous_point_depth_model_path', required=False, type=str)
args = parser.parse_args()
bam = args.bam
chr = args.chr
if "chr" in chr:
    chr = chr[3:]
begin = args.begin
end = args.end
outpath = args.outpath

#bam = "/mnt/X500/farmers/wangshj/atruepeopledepth/170006659BD_normal_sort_markdup_realign_recal.bam"
#chr = 1
#begin = 10318109
#end = 10319065

outdir = re.search(r'[^/]+$', bam).group(0)[:-4]
outdir = "contiguousDepth_" + outdir + "_" + str(chr) + "_" + str(begin) + "_" + str(end)

#outpath = None
if outpath != None:
        out = outpath
        out = out.strip()
        out = out.rstrip("\/")
        isExists=os.path.exists(out)
        if not isExists:
            os.makedirs(out)
else:
    out = outdir
    os.makedirs(out)
out = out + '/'

write_bed = open(out + "tmp_cal_depth.bed", "ab")
to_write_bed = chr + '\t' + str(begin) + '\t' + str(end) + '\n' + 'chr' + chr + '\t' + str(begin) + '\t' + str(end)
write_bed.write(to_write_bed)
write_bed.close()

command = "samtools depth -b " + out + "tmp_cal_depth.bed " + bam + " > " + out + "tmp_contiguous_point_depth.txt"
os.system(command)

file = open(out + "tmp_contiguous_point_depth.txt")
lines = file.readlines()
newfile = open(out + "contiguous_point_depth.txt", "ab")
for line in lines:
    newline = re.search(r'[^\t]+$', line.strip('\n')).group(0)
    newline = newline + "\n"
    newfile.write(newline)
file.close()
newfile.close()

os.remove(out + "tmp_cal_depth.bed")
os.remove(out + "tmp_contiguous_point_depth.txt")

newfile = open(out + "contiguous_point_depth.txt")
lines = newfile.readlines()
newfile.close()
contiguous_point_depth_density = []
aa = []
for line in lines:
    aa.append(int(line.strip('\n')))
i = 0
for a in aa:
    i = i + 1
    for k in range(a):
        contiguous_point_depth_density.append(i / float(len(aa)))
contiguous_point_depth_density = np.array(contiguous_point_depth_density)
k0 = [0.6, 0.4]
mu0 = [0.45, 0.55]
sigma0 = [0.12, 0.12]
step = 6
k, mu, sigma = em(contiguous_point_depth_density, k0, mu0, sigma0, step)

newfile = open(out + "contiguous_point_depth.txt")
lines = newfile.readlines()
newfile.close()
a = []
for line in lines:
    a.append(int(line.strip('\n')))
alength = len(a)
flag = False
for i in range(100, alength - 101):
    if (a[i] < a[i - 100]) and (a[i] < a[i + 100]):
        flag = True
        break

write_parameter = open(out + "parameter.txt", "ab")
if flag:
    first_line = "multiple"
else:
    first_line = "double"
if str(k[0]) == "nan" or str(k[1]) == "nan" or str(mu[0]) == "nan" or str(mu[1]) == "nan" or str(sigma[0]) == "nan" or str(sigma[1]) == "nan":
    first_line = "multiple"
second_line = str(k[0]) + "\t" + str(k[1])
third_line = str(mu[0]) + "\t" + str(mu[1])
last_line = str(sigma[0]) + "\t" + str(sigma[1])
to_write_parameter = first_line + "\n" + second_line + "\n" + third_line + "\n" + last_line
write_parameter.write(to_write_parameter)
write_parameter.close()
