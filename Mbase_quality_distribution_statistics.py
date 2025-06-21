# coding:utf-8
import os
import time
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-txt', '--txt', help='the input_base_quality_value_txt_file_path,eg:/mnt/GenePlus001/archive/2018/07/30/H52W5BGX7_L001/180730_TPNB500129_H52W5BGX7_L001_HUMCAAXH0PEI-108-mini/Base_quality_value_distribution_by_read_position_1.txt', required=True, type=str)
parser.add_argument('-out', '--outpath', help='the out_base_quality_distribution__model_path', required=False, type=str)

args = parser.parse_args()
txt = args.txt
outpath = args.outpath

if outpath != None:
        out = outpath
        out = out.strip()
        out = out.rstrip("\/")
        isExists=os.path.exists(out)
        if not isExists:
            os.makedirs(out)
else:
    out = "base_quality_distribution_model" + time.strftime('%Y%m%d%H%M%S', time.localtime(time.time()))
    os.makedirs(out)
out = out + '/'

quality_arr = []
for record in open(txt, 'r').readlines():
    try:
        arr = record.split('\t')
        quality_arr.append(arr[-5])
    except:
        break
quality_arr.pop()
del quality_arr[0]
#print len(quality_arr)
#print quality_arr

quality_distribution = []
quality_last = float(quality_arr[-1])
for quality in quality_arr:
    quality_distribution.append(round(float(quality) / quality_last, 3))
quality_distribution = np.array(quality_distribution)

np.savetxt(out + 'base_quality_distribution.txt', quality_distribution)
