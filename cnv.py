# -*- coding: UTF-8 -*-

def cnv_multiple(chro, chro_pos, customized_dict):
    flag = 0
    for cnv in customized_dict["cnv"]:
        if chro == cnv["chrom"]:
            if cnv["start"] < chro_pos < cnv["end"]:
                flag = 1
                return cnv["copy"]
    if flag == 0:
        return (1.0)