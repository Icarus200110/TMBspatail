# -*- coding: UTF-8 -*-
import random
from complementary import complementary
import re
import globalname
from polymorphic_adjustment import polymorphic_adjust

def inversion_imitate(begin, ini_sequence, chromo, customized_inversion, hg19index_dict):
    for inversion in customized_inversion:
        if (chromo == inversion["chrom"]) & (inversion["end"] >= begin) & (inversion["start"] < (begin + len(ini_sequence))) & (random.random() < inversion["rate"]):
            inversion_left = max(inversion["start"] - begin, 0)
            l = ini_sequence[0:inversion_left]
            inversion_right = min(inversion["end"] - (begin + len(ini_sequence)), 0)
            if inversion_right == 0:
                r = ""
            else:
                r = ini_sequence[inversion_right:]
            #symmetric_sum = inversion["start"] + inversion["end"]
            #middle = complementary(ini_sequence[inversion_left:len(ini_sequence) + inversion_right][::-1])
            #middle = complementary(ini_sequence[symmetric_sum - len(ini_sequence) - inversion_right:symmetric_sum - inversion_left][::-1])
            middle_begin = inversion["start"] + inversion["end"] - len(ini_sequence) - inversion_right - begin
            middle_length = len(ini_sequence) + inversion_right - inversion_left
            hg_fa_file = file("hg19.fa", 'r')
            if (middle_begin % 50) == 0:
                hg_fa_file.seek(hg19index_dict[chromo][0] + middle_begin + (middle_begin / 50) - 1, 0)
            else:
                hg_fa_file.seek(hg19index_dict[chromo][0] + middle_begin + (middle_begin / 50), 0)
            middle = hg_fa_file.read(middle_length + (middle_length / 51) + 1)
            middle = re.sub('\s', '', middle)
            middle = middle[0:middle_length]
            middle = middle.upper()
            hg_fa_file.close()
            if globalname.frequence_type == "AF" or globalname.frequence_type == "EAS_AF" or globalname.frequence_type == "AMR_AF" or globalname.frequence_type == "AFR_AF" or globalname.frequence_type == "EUR_AF" or globalname.frequence_type == "SAS_AF":
                middle = polymorphic_adjust(chromo, middle_begin, middle)
            middle = complementary(middle[::-1])
            ini_sequence = l + middle + r
    return ini_sequence
