# -*- coding: UTF-8 -*-
import globalname
import re
import random
from complementary import complementary
from polymorphic_adjustment import polymorphic_adjust


def if_is_fused(customized_fusion_one, chromo, begin, length):
    if customized_fusion_one["strand"] == "+":
        return (customized_fusion_one["chrom"] == chromo) & (begin < customized_fusion_one["pos"]) & (customized_fusion_one["pos"] < begin + length)
    else:
        return (customized_fusion_one["chrom"] == chromo) & (begin - length < customized_fusion_one["pos"]) & (customized_fusion_one["pos"] < begin)


def load_tofuse_sequence_sides(direction, fuse_direction, hg19index_dict):
    hg_fa_file = file("hg19.fa", 'r')
    begin = fuse_direction["pos"] - direction + 1
    if (begin % 50) == 0:
        hg_fa_file.seek(hg19index_dict[fuse_direction["chrom"]][0] + begin + (begin / 50) - 1, 0)
    else:
        hg_fa_file.seek(hg19index_dict[fuse_direction["chrom"]][0] + begin + (begin / 50), 0)
    tofuse_sequence = hg_fa_file.read(direction + (direction / 51) + 1)
    tofuse_sequence = re.sub('\s', '', tofuse_sequence)
    tofuse_sequence = tofuse_sequence[0:direction]
    tofuse_sequence = tofuse_sequence.upper()
    if globalname.frequence_type == "AF" or globalname.frequence_type == "EAS_AF" or globalname.frequence_type == "AMR_AF" or globalname.frequence_type == "AFR_AF" or globalname.frequence_type == "EUR_AF" or globalname.frequence_type == "SAS_AF":
        tofuse_sequence = polymorphic_adjust(fuse_direction["chrom"], begin, tofuse_sequence)
    return tofuse_sequence
    hg_fa_file.close()


def load_tofuse_sequence_middles(direction, fuse_direction, hg19index_dict):
    hg_fa_file = file("hg19.fa", 'r')
    begin = fuse_direction["pos"]
    if (begin % 50) == 0:
        hg_fa_file.seek(hg19index_dict[fuse_direction["chrom"]][0] + begin + (begin / 50) - 1, 0)
    else:
        hg_fa_file.seek(hg19index_dict[fuse_direction["chrom"]][0] + begin + (begin / 50), 0)
    tofuse_sequence = hg_fa_file.read(direction + (direction / 51) + 1)
    tofuse_sequence = re.sub('\s', '', tofuse_sequence)
    tofuse_sequence = tofuse_sequence[0:direction]
    tofuse_sequence = tofuse_sequence.upper()
    if globalname.frequence_type == "AF" or globalname.frequence_type == "EAS_AF" or globalname.frequence_type == "AMR_AF" or globalname.frequence_type == "AFR_AF" or globalname.frequence_type == "EUR_AF" or globalname.frequence_type == "SAS_AF":
        tofuse_sequence = polymorphic_adjust(fuse_direction["chrom"], begin, tofuse_sequence)
    return tofuse_sequence
    hg_fa_file.close()


def fuse_sequence_making(left, right, fuse_left, fuse_right, hg19index_dict, if_insert):
    #sequence_l = ""              #不知道是否需要这两行
    #sequence_r = ""
    #print fuse_right, fuse_left
    if fuse_left["strand"] == "+":
        sequence_l = load_tofuse_sequence_sides(left, fuse_left, hg19index_dict)
    else:
        sequence_l = load_tofuse_sequence_middles(left, fuse_left, hg19index_dict)
        sequence_l = complementary(sequence_l[::-1])
    if fuse_right["strand"] == "+":
        sequence_r = load_tofuse_sequence_middles(right, fuse_right, hg19index_dict)
    else:
        sequence_r = load_tofuse_sequence_sides(right, fuse_right, hg19index_dict)
        sequence_r = complementary(sequence_r[::-1])
    if if_insert:
        sequence_fused = sequence_l + globalname.insert_fragment + sequence_r
        #print fuse_left, fuse_right
        #print globalname.insert_fragment
    else:
        sequence_fused = sequence_l + sequence_r
        #print fuse_left, fuse_right
    #print globalname.insert_fragment
    #print sequence_fused
    return sequence_fused


def fuse_imitate(chromo, ini_sequence, begin, customized_fusion, hg19index_dict):
    length = len(ini_sequence)
    left = 0
    right = 0
    for fuse in customized_fusion:
        fuse_flag = False
        if if_is_fused(fuse["left"], chromo, begin, length):
            fuse_flag = True
            if fuse["left"]["strand"] == "+":
                left = fuse["left"]["pos"] - begin
            else:
                left = begin - fuse["left"]["pos"]
            left = min(length - 1, max(1, left))
            right = length - left
        elif if_is_fused(fuse["right"], chromo, begin, length):
            fuse_flag = True
            if fuse["right"]["strand"] == "+":
                right = (begin + length) - fuse["right"]["pos"]
            else:
                right = fuse["right"]["pos"] - (begin - length)
            left = length - right
        if fuse_flag & (left > 0) & (right > 0) & (random.random() < fuse["rate"]):
            return fuse_sequence_making(left, right, fuse["left"], fuse["right"], hg19index_dict, fuse["if_insert"])
    return ini_sequence
