# -*- coding: UTF-8 -*-
import re
import random
import globalname
from polymorphic_adjustment import polymorphic_adjust


def load_initial_sequence(hg19index_dict, position, chromo, templet_length, customized_delete, customized_tandem_repeat):
    hg_fa_file = open("hg19.fa", 'r')
    if not chromo in hg19index_dict:
        print("没有找到该染色体！")
    begin = max(1, position - int(templet_length / 2))
    begin = min(begin, hg19index_dict[chromo][1] - templet_length)
    has_delete = False
    has_tandem_repeat = False
    for delete in customized_delete:
        delete_end = delete["end"] + 2
        if (chromo == delete["chrom"]) & (delete_end > begin) & (delete["start"] + 1 < (begin + templet_length)):
            has_delete = True
    for tandem_repeat in customized_tandem_repeat:
        tandem_repeat_end = tandem_repeat["end"] + 1
        if (chromo == tandem_repeat["chrom"]) & (tandem_repeat_end > begin) & (tandem_repeat_end < (begin + templet_length)):
            has_tandem_repeat = True
    if has_delete:
        for delete in customized_delete:
            delete_end = delete["end"] + 2
            if (chromo == delete["chrom"]) & (delete_end > begin) & (delete["start"] + 1 < (begin + templet_length)) & (random.random() < delete["rate"]):
                #print 1
                #print delete["start"] - templet_length
                double_templet_l_begin = max(1, delete["start"] + 1 - templet_length)
                double_templet_l_length = delete["start"] + 1 - double_templet_l_begin
                #double_templet_r_end = delete["end"] + templet_length
                #double_templet_r_length = double_templet_r_end - delete["end"]
                if (double_templet_l_begin % 50) == 0:
                    hg_fa_file.seek(hg19index_dict[chromo][0] + double_templet_l_begin + (double_templet_l_begin / 50) - 1, 0)
                else:
                    hg_fa_file.seek(hg19index_dict[chromo][0] + double_templet_l_begin + (double_templet_l_begin / 50), 0)
                double_templet_l = hg_fa_file.read(double_templet_l_length + (double_templet_l_length / 51) + 1)
                double_templet_l = re.sub('\s', '', double_templet_l)
                double_templet_l = double_templet_l[0:double_templet_l_length]
                if globalname.frequence_type == "AF" or globalname.frequence_type == "EAS_AF" or globalname.frequence_type == "AMR_AF" or globalname.frequence_type == "AFR_AF" or globalname.frequence_type == "EUR_AF" or globalname.frequence_type == "SAS_AF":
                    double_templet_l = polymorphic_adjust(chromo, double_templet_l_begin, double_templet_l)
                #double_templet_l = double_templet_l.upper()

                if (delete_end % 50) == 0:
                    hg_fa_file.seek(hg19index_dict[chromo][0] + delete_end + (delete_end / 50) - 1, 0)
                else:
                    hg_fa_file.seek(hg19index_dict[chromo][0] + delete_end + (delete_end / 50), 0)
                double_templet_r = hg_fa_file.read(templet_length + (templet_length / 51) + 1)
                double_templet_r = re.sub('\s', '', double_templet_r)
                double_templet_r = double_templet_r[0:templet_length]
                if globalname.frequence_type == "AF" or globalname.frequence_type == "EAS_AF" or globalname.frequence_type == "AMR_AF" or globalname.frequence_type == "AFR_AF" or globalname.frequence_type == "EUR_AF" or globalname.frequence_type == "SAS_AF":
                    double_templet_r = polymorphic_adjust(chromo, delete_end, double_templet_r)
                #double_templet_r = double_templet_r.upper()
                hg_fa_file.close()
                #begin = min(position, delete["start"])
                #ini_sequence = double_templet[begin:begin+templet_length]
                if begin >= delete["start"] + 1:
                    #print 11
                    if delete["if_insert"]:
                        ini_sequence = globalname.insert_fragment + double_templet_r
                    else:
                        ini_sequence = double_templet_r
                else:
                    #print 22
                    if delete["if_insert"]:
                        double_templet = double_templet_l + globalname.insert_fragment + double_templet_r
                        ini_sequence = double_templet[(begin - double_templet_l_begin):(begin - double_templet_l_begin + len(globalname.insert_fragment) + templet_length)]
                    else:
                        double_templet = double_templet_l + double_templet_r
                        ini_sequence = double_templet[(begin - double_templet_l_begin):(begin - double_templet_l_begin + templet_length)]
                ini_sequence = ini_sequence.upper()
                return begin, ini_sequence

    if has_tandem_repeat:
        for tandem_repeat in customized_tandem_repeat:
            tandem_repeat_end = tandem_repeat["end"] + 1
            if (chromo == tandem_repeat["chrom"]) & (tandem_repeat_end > begin) & (tandem_repeat_end < (begin + templet_length)) & (random.random() < tandem_repeat["rate"]):
                #print 2
                double_templet_l_begin = max(1, tandem_repeat["start"] - templet_length)
                double_templet_l_length = tandem_repeat["start"] - double_templet_l_begin
                if (double_templet_l_begin % 50) == 0:
                    hg_fa_file.seek(hg19index_dict[chromo][0] + double_templet_l_begin + (double_templet_l_begin / 50) - 1, 0)
                else:
                    hg_fa_file.seek(hg19index_dict[chromo][0] + double_templet_l_begin + (double_templet_l_begin / 50), 0)
                double_templet_l = hg_fa_file.read(double_templet_l_length + (double_templet_l_length / 51) + 1)
                double_templet_l = re.sub('\s', '', double_templet_l)
                double_templet_l = double_templet_l[0:double_templet_l_length]
                if globalname.frequence_type == "AF" or globalname.frequence_type == "EAS_AF" or globalname.frequence_type == "AMR_AF" or globalname.frequence_type == "AFR_AF" or globalname.frequence_type == "EUR_AF" or globalname.frequence_type == "SAS_AF":
                    double_templet_l = polymorphic_adjust(chromo, double_templet_l_begin, double_templet_l)

                tandem_repeat_begin = tandem_repeat["start"]
                tandem_repeat_length = tandem_repeat["end"] - tandem_repeat_begin + 1
                if (tandem_repeat_begin % 50) == 0:
                    hg_fa_file.seek(hg19index_dict[chromo][0] + tandem_repeat_begin + (tandem_repeat_begin / 50) - 1, 0)
                else:
                    hg_fa_file.seek(hg19index_dict[chromo][0] + tandem_repeat_begin + (tandem_repeat_begin / 50), 0)
                double_templet_m = hg_fa_file.read(tandem_repeat_length + (tandem_repeat_length / 51) + 1)
                double_templet_m = re.sub('\s', '', double_templet_m)
                double_templet_m = double_templet_m[0:tandem_repeat_length]
                if globalname.frequence_type == "AF" or globalname.frequence_type == "EAS_AF" or globalname.frequence_type == "AMR_AF" or globalname.frequence_type == "AFR_AF" or globalname.frequence_type == "EUR_AF" or globalname.frequence_type == "SAS_AF":
                    double_templet_m = polymorphic_adjust(chromo, tandem_repeat_begin, double_templet_m)

                if (tandem_repeat_end % 50) == 0:
                    hg_fa_file.seek(hg19index_dict[chromo][0] + tandem_repeat_end + (tandem_repeat_end / 50) - 1, 0)
                else:
                    hg_fa_file.seek(hg19index_dict[chromo][0] + tandem_repeat_end + (tandem_repeat_end / 50), 0)
                double_templet_r = hg_fa_file.read(templet_length + (templet_length / 51) + 1)
                double_templet_r = re.sub('\s', '', double_templet_r)
                double_templet_r = double_templet_r[0:templet_length]
                if globalname.frequence_type == "AF" or globalname.frequence_type == "EAS_AF" or globalname.frequence_type == "AMR_AF" or globalname.frequence_type == "AFR_AF" or globalname.frequence_type == "EUR_AF" or globalname.frequence_type == "SAS_AF":
                    double_templet_r = polymorphic_adjust(chromo, tandem_repeat_end, double_templet_r)
                #double_templet_r = double_templet_r.upper()
                hg_fa_file.close()

                tandem_repeat_copy = tandem_repeat["copy"]
                for i in range(tandem_repeat_copy):
                    double_templet_l = double_templet_l + double_templet_m
                double_templet = double_templet_l + double_templet_r
                random_position = random.randint(1, len(double_templet_l) + 1)
                ini_sequence = double_templet[random_position - 1:(random_position + templet_length)]
                ini_sequence = ini_sequence.upper()
                begin = double_templet_l_begin + random_position - 1
                if (random_position + templet_length) > len(double_templet_l):
                    #begin = tandem_repeat["end"] + templet_length - (len(double_templet) - (random_position + templet_length))
                    begin = tandem_repeat["end"] + 2 * templet_length + random_position - len(double_templet)
                elif (len(double_templet_l) - (tandem_repeat_copy - 1) * len(double_templet_m)) < random_position < len(double_templet_l):
                    begin = -100000000
                return begin, ini_sequence


    if (begin % 50) == 0:
        hg_fa_file.seek(hg19index_dict[chromo][0] + begin + (begin / 50) - 1, 0)
    else:
        hg_fa_file.seek(hg19index_dict[chromo][0] + begin + (begin / 50), 0)
    ini_sequence = hg_fa_file.read(int(templet_length + (templet_length / 51) + 1))
    ini_sequence = re.sub('\s', '', ini_sequence)
    ini_sequence = ini_sequence[0:templet_length]
    ini_sequence = ini_sequence.upper()
    hg_fa_file.close()
    if globalname.frequence_type == "AF" or globalname.frequence_type == "EAS_AF" or globalname.frequence_type == "AMR_AF" or globalname.frequence_type == "AFR_AF" or globalname.frequence_type == "EUR_AF" or globalname.frequence_type == "SAS_AF":
        ini_sequence = polymorphic_adjust(chromo, begin, ini_sequence)
    return begin, ini_sequence

