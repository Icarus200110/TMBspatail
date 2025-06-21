import re
import random
import globalname
from polymorphic_adjustment import polymorphic_adjust

def insert_imitate(hg19index_dict, begin, ini_sequence, chromo, customized_insert):
    ini_sequence_length = len(ini_sequence)
    for insert in customized_insert:
        if (chromo == insert["chrom"]) & (begin <= insert["pos"]) & (begin + ini_sequence_length > insert["pos"]) & (random.random() < insert["rate"]):
            hg_fa_file = file("hg19.fa", 'r')
            double_templet_l_begin = insert["pos"] - ini_sequence_length + 2
            double_templet_l_length = ini_sequence_length - 1
            if (double_templet_l_begin % 50) == 0:
                hg_fa_file.seek(hg19index_dict[chromo][0] + double_templet_l_begin + (double_templet_l_begin / 50) - 1, 0)
            else:
                hg_fa_file.seek(hg19index_dict[chromo][0] + double_templet_l_begin + (double_templet_l_begin / 50), 0)
            double_templet_l = hg_fa_file.read(double_templet_l_length + (double_templet_l_length / 51) + 1)
            double_templet_l = re.sub('\s', '', double_templet_l)
            double_templet_l = double_templet_l[0:double_templet_l_length]
            if globalname.frequence_type == "AF" or globalname.frequence_type == "EAS_AF" or globalname.frequence_type == "AMR_AF" or globalname.frequence_type == "AFR_AF" or globalname.frequence_type == "EUR_AF" or globalname.frequence_type == "SAS_AF":
                double_templet_l = polymorphic_adjust(chromo, double_templet_l_begin, double_templet_l)

            double_templet_r_begin = insert["pos"] + 1
            double_templet_r_length = ini_sequence_length - 1
            if (double_templet_r_begin % 50) == 0:
                hg_fa_file.seek(hg19index_dict[chromo][0] + double_templet_r_begin + (double_templet_r_begin / 50) - 1, 0)
            else:
                hg_fa_file.seek(hg19index_dict[chromo][0] + double_templet_r_begin + (double_templet_r_begin / 50), 0)
            double_templet_r = hg_fa_file.read(double_templet_r_length + (double_templet_r_length / 51) + 1)
            double_templet_r = re.sub('\s', '', double_templet_r)
            double_templet_r = double_templet_r[0:double_templet_r_length]
            if globalname.frequence_type == "AF" or globalname.frequence_type == "EAS_AF" or globalname.frequence_type == "AMR_AF" or globalname.frequence_type == "AFR_AF" or globalname.frequence_type == "EUR_AF" or globalname.frequence_type == "SAS_AF":
                double_templet_r = polymorphic_adjust(chromo, double_templet_r_begin, double_templet_r)

            hg_fa_file.close()

            double_templet = double_templet_l + globalname.insert_fragment + double_templet_r
            random_position = random.randint(1, len(globalname.insert_fragment) + double_templet_l_length)
            ini_sequence = double_templet[random_position - 1:(random_position + ini_sequence_length)]
            ini_sequence = ini_sequence.upper()

            begin = double_templet_l_begin + random_position - 1
            if (random_position > double_templet_l_length) & (random_position + ini_sequence_length <= double_templet_l_length + len(globalname.insert_fragment)):
                begin = -100000000
    return begin, ini_sequence
