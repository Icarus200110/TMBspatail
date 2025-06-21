# -*- coding: UTF-8 -*-
import random
from random import choice

def mixgauss(k,mu,sigma):
    n = len(k)
    rand = random.random()
    Sum = 0
    index = 0
    while (index < n):
        Sum += k[index]
        if (rand < Sum):
            return random.gauss(mu[index], sigma[index])
        else:
            index += 1

def capture_to_chromo(target_region_list,target_region_pos,contiguous_point_depth_densityorparameter):
    #file = open('2222.txt')
    #lines = file.readlines()
    #X = []
    #aa = []
    #for line in lines:
    #    aa.append(int(line.strip('\n')))
    #i = 0
    #for a in aa:
    #    i = i + 1
    #    for k in range(a):
    #        X.append(i / float(len(aa)))
    k = [0.5386058564513958, 0.4613941435486041]
    mu = [0.4643908657876027, 0.5527932249173736]
    sigma = [0.11646926461797634, 0.1237254958583397]
    l = 0
    r = len(target_region_list) - 1
    temp = 0
    while l < r:
        mid = int((l + r) / 2)
        onechrdict_in_list = target_region_list[mid]
        if (onechrdict_in_list["target_list_begin"] <= target_region_pos) & (
            onechrdict_in_list["target_list_end"] >= target_region_pos):
            temp = mid
            break
        elif (onechrdict_in_list["target_list_begin"] > target_region_pos):
            r = mid - 1
        elif (onechrdict_in_list["target_list_end"] < target_region_pos):
            l = mid + 1
    if temp == 0:
        temp = l
    onechrdict_in_list = target_region_list[temp]
    #neutralization = target_region_pos - onechrdict_in_list["target_list_begin"]
    #target_chro_pos = onechrdict_in_list["begin"] + neutralization     # 将其对应到具体染色体上的位置
    #onechrdict_in_list_begin = int(round((onechrdict_in_list["end"] - onechrdict_in_list["begin"]) * 0.5)) + onechrdict_in_list["begin"]
    #target_chro_pos = onechrdict_in_list_begin + int(round((onechrdict_in_list["end"] - onechrdict_in_list_begin) * random.gauss(0, 1) / 4))
    onechrdict_in_list_begin = onechrdict_in_list["begin"]
    if contiguous_point_depth_densityorparameter != None:
        if contiguous_point_depth_densityorparameter[0] == "double":
            second_line_arr = contiguous_point_depth_densityorparameter[1].split("\t")
            k[0] = float(second_line_arr[0])
            k[1] = float(second_line_arr[1])
            third_line_arr = contiguous_point_depth_densityorparameter[2].split("\t")
            mu[0] = float(third_line_arr[0])
            mu[1] = float(third_line_arr[1])
            last_line_arr = contiguous_point_depth_densityorparameter[3].split("\t")
            sigma[0] = float(last_line_arr[0])
            sigma[1] = float(last_line_arr[1])
            target_chro_pos = onechrdict_in_list_begin + int(round((onechrdict_in_list["end"] - onechrdict_in_list_begin) * mixgauss(k,mu,sigma)))
        else:
            target_chro_pos = onechrdict_in_list_begin + int(round((onechrdict_in_list["end"] - onechrdict_in_list_begin) * choice(contiguous_point_depth_densityorparameter)))
    else:
        target_chro_pos = onechrdict_in_list_begin + int(round((onechrdict_in_list["end"] - onechrdict_in_list_begin) * mixgauss(k,mu,sigma)))
    return (target_chro_pos, onechrdict_in_list["chro"])
