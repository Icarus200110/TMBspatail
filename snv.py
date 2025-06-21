# -*- coding: UTF-8 -*-
import random

def snv_imitate(begin, ini_sequence, chromo, customized_snv):
    for snv in customized_snv:
        if (chromo == snv["chrom"]) & (begin <= snv["pos"]) & (begin + len(ini_sequence) > snv["pos"]) & (random.random() < snv["rate"]):
                neutralize = snv["pos"] - begin
                l = ini_sequence[0:neutralize]
                middle = snv["alt"]
                r = ini_sequence[neutralize + len(snv["ref"]):]
                ini_sequence = l + middle + r
    return ini_sequence
