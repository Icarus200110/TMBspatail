# -*- coding: UTF-8 -*
import re
import random
from random import choice


def interfere_mutation(original_mutation, sample_serial_num):
    delete_array = ["tandem_repeat", "insert", "cnv", "snv", "inversion", "fusion"]
    tandem_repeat_array = ["delete", "insert", "cnv", "snv", "inversion", "fusion"]
    insert_array = ["tandem_repeat", "delete", "cnv", "snv", "inversion", "fusion"]
    cnv_array = ["tandem_repeat", "insert", "delete", "snv", "inversion", "fusion"]
    snv_array = ["tandem_repeat", "insert", "cnv", "delete", "inversion", "fusion"]
    inversion_array = ["tandem_repeat", "insert", "cnv", "snv", "delete", "fusion"]
    fusion_array = ["tandem_repeat", "insert", "cnv", "snv", "inversion", "delete"]
    if original_mutation == "delete":
        return mutation_making(delete_array[random.randint(0, 5)], sample_serial_num)
    if original_mutation == "tandem_repeat":
        return mutation_making(tandem_repeat_array[random.randint(0, 5)], sample_serial_num)
    if original_mutation == "insert":
        return mutation_making(insert_array[random.randint(0, 5)], sample_serial_num)
    if original_mutation == "cnv":
        return mutation_making(cnv_array[random.randint(0, 5)], sample_serial_num)
    if original_mutation == "snv":
        return mutation_making(snv_array[random.randint(0, 5)], sample_serial_num)
    if original_mutation == "inversion":
        return mutation_making(inversion_array[random.randint(0, 5)], sample_serial_num)
    if original_mutation == "fusion":
        return mutation_making(fusion_array[random.randint(0, 5)], sample_serial_num)


def mutation_making(mutation, sample_serial_num):
    chr_dict = {
        "chr1": [10000, 249240621],
        "chr2": [10000, 243198373],
        "chr3": [60000, 197962430],
        "chr4": [10000, 191044276],
        "chr5": [10000, 180905260],
        "chr6": [60000, 171055067],
        "chr7": [10000, 159137663],
        "chr8": [10000, 146304022],
        "chr9": [10000, 141153431],
        "chr10": [60000, 135524747],
        "chr11": [60000, 134946516],
        "chr12": [60000, 133841895],
        "chr13": [19020000, 115109878],
        "chr14": [10500000, 107289540],
        "chr15": [20000000, 102521392],
        "chr16": [60000, 90294753],
        "chr17": [1, 81195210],
        "chr18": [10000, 78017248],
        "chr19": [60000, 59118983],
        "chr20": [210000, 62965520],
        "chr21": [9411200, 48119895],
        "chr22": [16050000, 51244566],
        "chrX": [60000, 155260560],
        "chrY": [10000, 59363566]
    }
    chr_list = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
                "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX",
                "chrY"]
    if mutation == "delete":
        delete_min = 1
        delete_max = 1000
        chr_rnd = choice(chr_list)
        start_rnd = random.randint(chr_dict[chr_rnd][0], chr_dict[chr_rnd][1])
        delete_len_rnd = random.randint(delete_min - 1, delete_max - 1)
        end_rnd = min(start_rnd + delete_len_rnd, chr_dict[chr_rnd][1])
        if_insert = random.randint(0, 1)
        rate = round(random.random(), 4)
        command_str = " -delete onedelete," + chr_rnd + "," + str(start_rnd) + "," + str(end_rnd) + "," + str(
            rate) + "," + str(if_insert)
        tsv_str = "simu_sample" + str(sample_serial_num + 1) + '\t' + "delete" + '\t' + chr_rnd + '\t' + str(start_rnd) + '\t' + str(end_rnd) + '\t' + str(
            rate) + '\t' + str(if_insert) + '\n'
    if mutation == "tandem_repeat":
        tandem_repeat_min = 1
        tandem_repeat_max = 1000
        copy_min = 2
        copy_max = 10
        copy_num = random.randint(copy_min, copy_max)
        chr_rnd = choice(chr_list)
        start_rnd = random.randint(chr_dict[chr_rnd][0], chr_dict[chr_rnd][1])
        tandem_repeat_len_rnd = random.randint(tandem_repeat_min - 1, tandem_repeat_max - 1)
        end_rnd = min(start_rnd + tandem_repeat_len_rnd, chr_dict[chr_rnd][1])
        rate = round(random.random(), 4)
        command_str = " -tandem_repeat onetandem_repeat," + chr_rnd + "," + str(start_rnd) + "," + str(end_rnd) + "," + str(
            copy_num) + "," + str(rate)
        tsv_str = "simu_sample" + str(sample_serial_num + 1) + '\t' + "tandem_repeat" + '\t' + chr_rnd + '\t' + str(
            start_rnd) + '\t' + str(end_rnd) + '\t' + str(copy_num) + '\n'
    if mutation == "insert":
        chr_rnd = choice(chr_list)
        pos_rnd = random.randint(chr_dict[chr_rnd][0], chr_dict[chr_rnd][1])
        rate = round(random.random(), 4)
        command_str = " -insert oneinsert," + chr_rnd + "," + str(pos_rnd) + "," + str(rate)
        tsv_str = "simu_sample" + str(sample_serial_num + 1) + '\t' + "insert" + '\t' + chr_rnd + '\t' + str(pos_rnd) + '\t' + "unknownInsertLength" + '\n'
    if mutation == "cnv":
        cnv_min = 1
        cnv_max = 1000
        chr_rnd = choice(chr_list)
        start_rnd = random.randint(chr_dict[chr_rnd][0], chr_dict[chr_rnd][1])
        cnv_len_rnd = random.randint(cnv_min - 1, cnv_max - 1)
        end_rnd = min(start_rnd + cnv_len_rnd, chr_dict[chr_rnd][1])
        rate = round(random.uniform(0, 10), 1)
        command_str = " -cnv onecnv," + chr_rnd + "," + str(start_rnd) + "," + str(end_rnd) + "," + str(rate)
        tsv_str = "simu_sample" + str(sample_serial_num + 1) + '\t' + "cnv" + '\t' + chr_rnd + '\t' + str(
            start_rnd) + '\t' + str(end_rnd) + '\n'
    if mutation == "snv":
        substitutition = {
            'A': ['T', 'C', 'G'],
            'T': ['A', 'C', 'G'],
            'C': ['T', 'A', 'G'],
            'G': ['T', 'C', 'A']
        }
        hg19index_dict = {
            "chr1": [5, 249250621],
            "chr2": [254235645, 243199373],
            "chr3": [502299012, 198022430],
            "chr4": [704281897, 191154276],
            "chr5": [899259265, 180915260],
            "chr6": [1083792837, 171115067],
            "chr7": [1258330212, 159138663],
            "chr8": [1420651655, 146364022],
            "chr9": [1569942964, 141213431],
            "chr10": [1713980671, 135534747],
            "chr11": [1852226120, 135006516],
            "chr12": [1989932774, 133851895],
            "chr13": [2126461714, 115169878],
            "chr14": [2243934997, 107349540],
            "chr15": [2353431535, 102531392],
            "chr16": [2458013562, 90354753],
            "chr17": [2550175418, 81195210],
            "chr18": [2632994540, 78077248],
            "chr19": [2712633340, 59128983],
            "chr20": [2772944910, 63025520],
            "chr21": [2837230948, 48129895],
            "chr22": [2886323448, 51304566],
            "chrX": [2938654112, 155270560],
            "chrY": [3097030090, 59373566]
        }
        hg_fa_file = open("hg19.fa", 'r')
        chr_rnd = choice(chr_list)
        pos_rnd = random.randint(chr_dict[chr_rnd][0], chr_dict[chr_rnd][1])
        if (pos_rnd % 50) == 0:
            hg_fa_file.seek(hg19index_dict[chr_rnd][0] + pos_rnd + (pos_rnd / 50) - 1, 0)
        else:
            hg_fa_file.seek(hg19index_dict[chr_rnd][0] + pos_rnd + (pos_rnd / 50), 0)
        ref = hg_fa_file.read(2)
        ref = re.sub('\s', '', ref)
        ref = ref[0:1]
        ref = ref.upper()
        if ref == "N":
            alt = "N"
        else:
            substitute_index = int(round(random.random() * 2))
            alt = substitutition[ref][substitute_index]
        rate = round(random.random(), 4)
        command_str = " -snv onesnv," + chr_rnd + "," + str(pos_rnd) + "," + ref + "," + alt + "," + str(rate)
        tsv_str = "simu_sample" + str(sample_serial_num + 1) + '\t' + "snv" + '\t' + chr_rnd + '\t' + str(
            pos_rnd) + '\t' + ref + '\t' + alt + '\n'
    if mutation == "inversion":
        inversion_min = 1
        inversion_max = 1000
        chr_rnd = choice(chr_list)
        start_rnd = random.randint(chr_dict[chr_rnd][0], chr_dict[chr_rnd][1])
        inversion_len_rnd = random.randint(inversion_min - 1, inversion_max - 1)
        end_rnd = min(start_rnd + inversion_len_rnd, chr_dict[chr_rnd][1])
        rate = round(random.random(), 4)
        command_str = " -inversion oneinversion," + chr_rnd + "," + str(start_rnd) + "," + str(end_rnd) + "," + str(
            rate)
        tsv_str = "simu_sample" + str(sample_serial_num + 1) + '\t' + "inversion" + '\t' + chr_rnd + '\t' + str(
            start_rnd) + '\t' + str(end_rnd) + '\n'
    if mutation == "fusion":
        strand_list = ["+", "-"]
        l_chr_rnd = choice(chr_list)
        l_strand_rnd = choice(strand_list)
        l_pos_rnd = random.randint(chr_dict[l_chr_rnd][0], chr_dict[l_chr_rnd][1])
        r_chr_rnd = choice(chr_list)
        r_strand_rnd = choice(strand_list)
        r_pos_rnd = random.randint(chr_dict[r_chr_rnd][0], chr_dict[r_chr_rnd][1])
        if_insert = random.randint(0, 1)
        rate = round(random.random(), 4)
        command_str = " -fusion onefusion," + l_chr_rnd + "," + l_strand_rnd + "," + str(
            l_pos_rnd) + "," + r_chr_rnd + "," + r_strand_rnd + "," + str(r_pos_rnd) + "," + str(rate) + "," + str(
            if_insert)
        tsv_str = "simu_sample" + str(
            sample_serial_num + 1) + '\t' + "fusion" + '\t' + l_chr_rnd + '\t' + l_strand_rnd + '\t' + str(
            l_pos_rnd) + '\t' + r_chr_rnd + '\t' + r_strand_rnd + '\t' + str(r_pos_rnd) + '\t' + str(
            if_insert) + '\n'
    return command_str, tsv_str
