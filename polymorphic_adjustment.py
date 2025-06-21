# -*- coding: UTF-8 -*-

from base_frequence_making import simulate_base_frequence
from get_tgp_af import get_tgp_frequence
import globalname
#import json

def polymorphic_adjust(chromo, begin, ini_sequence):
    #json_file = "/mnt/X500/farmers/wangshj/makefrebetter/renqun.json"
    #with open(json_file, "rb") as load_f:  # ´òson?¼þ²¢ת»¯Ϊ?µäʽ
        #load_dict = json.load(load_f)
    #load_f.close()
    ini_sequence_list = list(ini_sequence)
    ini_sequence_list_length = len(ini_sequence_list)
    i = 0
    while i < ini_sequence_list_length:
        ref_length = 1
        simulate_onebase_frequence = simulate_base_frequence(chromo, begin + i)
        #if simulate_onebase_frequence > 0.001:
        #if (begin + i) not in load_dict[chromo[3:]]:
            #i = i + ref_length
            #continue
        if str(begin + i) in globalname.tgp_dict[chromo]:
            this_tgp_frequence = globalname.tgp_dict[chromo][str(begin + i)]
        else:
            this_tgp_frequence = get_tgp_frequence(chromo, begin + i, globalname.frequence_type)
            globalname.tgp_dict[chromo][str(begin + i)] = this_tgp_frequence
        #if get_tgp_frequence(chromo, begin + i, globalname.frequence_type) != "":
            #ref_length, alt, frequence = get_tgp_frequence(chromo, begin + i, globalname.frequence_type)
        if this_tgp_frequence != "":
            ref_length, alt, frequence = this_tgp_frequence
            double = ',' in alt
            if double:
                alt_array = alt.split(',')
                frequence_array = frequence.split(',')
                #if simulate_base_frequence(chromo, begin + i) <= float(frequence_array[0]):
                if simulate_onebase_frequence <= float(frequence_array[0]):
                    ini_sequence_list[i] = alt_array[0]
                #elif simulate_base_frequence(chromo, begin + i) <= (frequence_array[0] + frequence_array[1]):
                elif simulate_onebase_frequence <= (float(frequence_array[0]) + float(frequence_array[1])):
                    ini_sequence_list[i] = alt_array[1]
            else:
                #if simulate_base_frequence(chromo, begin + i) <= frequence:
                if simulate_onebase_frequence <= float(frequence):
                    ini_sequence_list[i] = alt
        if ref_length != 1:
            maximum = min(i + ref_length, ini_sequence_list_length)
            for j in range(i + 1, maximum):
                ini_sequence_list[j] = ""
        i = i + ref_length
    while "" in ini_sequence_list:
        ini_sequence_list.remove("")
    ini_sequence = ''.join(ini_sequence_list)
    return ini_sequence
