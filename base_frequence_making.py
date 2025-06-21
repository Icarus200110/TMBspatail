# -*- coding: UTF-8 -*-

from __future__ import division
import globalname


#var
index = 624
MT = [0]*index
# MT[0] ->seed

def conversion(t):
    return(0xFFFFFFFF & t) #ȡ×ºó»->t

def rotate():
    global index
    for i in range(624):
        y = conversion((MT[i] & 0x80000000) + (MT[(i + 1) % 624] & 0x7fffffff))
        MT[i] = MT[(i + 397) % 624] ^ y >> 1
        if y % 2 != 0:
            MT[i] = MT[i] ^ 0x9908b0df
    index = 0

def reckon():
    global index
    if index >= 624:
        rotate()
    y = MT[index]
    y = y ^ y >> 11
    y = y ^ y << 7 & 2636928640
    y = y ^ y << 15 & 4022730752
    y = y ^ y >> 18
    index = index + 1
    return conversion(y)

def mason(seed):
    MT[0] = seed
    for i in range(1, 624):
        MT[i] = conversion(1812433253 * (MT[i - 1] ^ MT[i - 1] >> 30) + i)
    return reckon()

def simulate_base_frequence(chromo, position):
    chromo_dict = {
        "chr1": 1,
        "chr2": 2,
        "chr3": 3,
        "chr4": 4,
        "chr5": 5,
        "chr6": 6,
        "chr7": 7,
        "chr8": 8,
        "chr9": 9,
        "chr10": 10,
        "chr11": 20,
        "chr12": 30,
        "chr13": 40,
        "chr14": 50,
        "chr15": 60,
        "chr16": 70,
        "chr17": 80,
        "chr18": 90,
        "chr19": 100,
        "chr20": 200,
        "chr21": 300,
        "chr22": 400,
        "chrX": 500,
        "chrY": 600
    }
    chr_pos_number = int(str(chromo_dict[chromo]) + str(position))
    serial_number = globalname.rdm + chr_pos_number
    rd = mason(serial_number) / (2 ** 32 - 1)
    #rd = 0.000000001
    return rd
