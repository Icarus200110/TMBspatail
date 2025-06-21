# -*- coding: UTF-8 -*
import json
import argparse
import os
import time

parser = argparse.ArgumentParser()
parser.add_argument('-bam', '--bam', help='the input bam files path', required=True, nargs = '+', type=str)
parser.add_argument('-bed', '--bed', help='the input bed file path', required=True, type=str)
parser.add_argument('-out', '--outpath', help='the overall_depth_distribution_model_path', required=False, type=str)

args = parser.parse_args()
bam = args.bam
bed = args.bed
outpath = args.outpath

if outpath != None:
        out = outpath
        out = out.strip()
        out = out.rstrip("\/")
        isExists=os.path.exists(out)
        if not isExists:
            os.makedirs(out)
else:
    out = "overall_depth_distribution_model" + time.strftime('%Y%m%d%H%M%S', time.localtime(time.time()))
    os.makedirs(out)
out = out + '/'

bams = ""
for bam_one in bam:
        position = bam.index(bam_one)
	if position == (len(bam) - 1):
		bams = bams + bam_one
	else:
		bams = bams + bam_one + " "
command = "cnvkit.py batch -n " + bams + " -t " + bed + " -d " + out + " -p " + str(len(bam))
os.system(command)

nu = 0
total = 0
for record in open(out + 'reference.cnn', 'r').readlines():
    try:
        p1, p2, p3, p4, p5, p6, p7 = record.split('\t')
        if p4 == "Antitarget":
            continue
        nu += (int(p3) - int(p2))
        total += float(p6) * (int(p3) - int(p2))
    except:
        continue
mean_depth = total / nu

depth_distribution = {'1': {}, '2': {}, '3': {}, '4': {}, '5': {}, '6': {}, '7': {}, '8': {}, '9': {}, '10': {}, '11': {}, '12': {}, '13': {}, '14': {}, '15': {}, '16': {}, '17': {}, '18': {}, '19': {}, '20': {}, '21': {}, '22': {}, 'X': {}, 'Y': {}}
for record in open(out + 'reference.cnn', 'r').readlines():
    try:
        p1, p2, p3, p4, p5, p6, p7 = record.split('\t')
        if p4 == "Antitarget":
            continue
        for i in range(int(p2), int(p3)):
            depth_distribution[p1][str(i)] = float(p6) / mean_depth
    except:
        continue

with open(out + "depth_distribution.json", "w") as dump_f:
    json.dump(depth_distribution, dump_f)
dump_f.close()
