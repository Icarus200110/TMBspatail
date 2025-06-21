# -*- coding: UTF-8 -*
import re
import os
import time
import argparse
import random
from random import choice
from interfere_mutation import interfere_mutation


parser = argparse.ArgumentParser()
parser.add_argument('-m', '--mutation', help='the mutation', required=False, type=str)
parser.add_argument('-g', '--gradient', help='the gradient', required=False, type=float)
parser.add_argument('-n', '--per_gradient_sample_num', help='the per_gradient_sample_num', required=False, type=int)
parser.add_argument('-mn', '--per_sample_max_mutation_num', help='the per_sample_max_mutation_num', required=False, type=int)
parser.add_argument('-out', '--outpath', help='the outpath', required=False, type=str)
parser.add_argument('-bed', '--bed_file', help='the bedpath with no blank line', required=False, type=str)

args = parser.parse_args()
mutation = args.mutation
gradient = args.gradient
per_gradient_sample_num = args.per_gradient_sample_num
per_sample_max_mutation_num = args.per_sample_max_mutation_num
outpath = args.outpath
bed_file = args.bed_file
if bed_file != None:
    tar_bed_arr = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
    bedfile = open(bed_file)
    try:
        bed = bedfile.read( )
    finally:
        bedfile.close( )
    bed_row_arr = bed.split('\n')
    for bed_row in bed_row_arr:
        every_row_arr = bed_row.split()
        if every_row_arr[0] == 'chr1' or every_row_arr[0] == '1':
            one_arr = []
            one_arr.append(int(every_row_arr[1]))
            one_arr.append(int(every_row_arr[2]))
            tar_bed_arr[0].append(one_arr)
        if every_row_arr[0] == 'chr2' or every_row_arr[0] == '2':
            one_arr = []
            one_arr.append(int(every_row_arr[1]))
            one_arr.append(int(every_row_arr[2]))
            tar_bed_arr[1].append(one_arr)
        if every_row_arr[0] == 'chr3' or every_row_arr[0] == '3':
            one_arr = []
            one_arr.append(int(every_row_arr[1]))
            one_arr.append(int(every_row_arr[2]))
            tar_bed_arr[2].append(one_arr)
        if every_row_arr[0] == 'chr4' or every_row_arr[0] == '4':
            one_arr = []
            one_arr.append(int(every_row_arr[1]))
            one_arr.append(int(every_row_arr[2]))
            tar_bed_arr[3].append(one_arr)
        if every_row_arr[0] == 'chr5' or every_row_arr[0] == '5':
            one_arr = []
            one_arr.append(int(every_row_arr[1]))
            one_arr.append(int(every_row_arr[2]))
            tar_bed_arr[4].append(one_arr)
        if every_row_arr[0] == 'chr6' or every_row_arr[0] == '6':
            one_arr = []
            one_arr.append(int(every_row_arr[1]))
            one_arr.append(int(every_row_arr[2]))
            tar_bed_arr[5].append(one_arr)
        if every_row_arr[0] == 'chr7' or every_row_arr[0] == '7':
            one_arr = []
            one_arr.append(int(every_row_arr[1]))
            one_arr.append(int(every_row_arr[2]))
            tar_bed_arr[6].append(one_arr)
        if every_row_arr[0] == 'chr8' or every_row_arr[0] == '8':
            one_arr = []
            one_arr.append(int(every_row_arr[1]))
            one_arr.append(int(every_row_arr[2]))
            tar_bed_arr[7].append(one_arr)
        if every_row_arr[0] == 'chr9' or every_row_arr[0] == '9':
            one_arr = []
            one_arr.append(int(every_row_arr[1]))
            one_arr.append(int(every_row_arr[2]))
            tar_bed_arr[8].append(one_arr)
        if every_row_arr[0] == 'chr10' or every_row_arr[0] == '10':
            one_arr = []
            one_arr.append(int(every_row_arr[1]))
            one_arr.append(int(every_row_arr[2]))
            tar_bed_arr[9].append(one_arr)
        if every_row_arr[0] == 'chr11' or every_row_arr[0] == '11':
            one_arr = []
            one_arr.append(int(every_row_arr[1]))
            one_arr.append(int(every_row_arr[2]))
            tar_bed_arr[10].append(one_arr)
        if every_row_arr[0] == 'chr12' or every_row_arr[0] == '12':
            one_arr = []
            one_arr.append(int(every_row_arr[1]))
            one_arr.append(int(every_row_arr[2]))
            tar_bed_arr[11].append(one_arr)
        if every_row_arr[0] == 'chr13' or every_row_arr[0] == '13':
            one_arr = []
            one_arr.append(int(every_row_arr[1]))
            one_arr.append(int(every_row_arr[2]))
            tar_bed_arr[12].append(one_arr)
        if every_row_arr[0] == 'chr14' or every_row_arr[0] == '14':
            one_arr = []
            one_arr.append(int(every_row_arr[1]))
            one_arr.append(int(every_row_arr[2]))
            tar_bed_arr[13].append(one_arr)
        if every_row_arr[0] == 'chr15' or every_row_arr[0] == '15':
            one_arr = []
            one_arr.append(int(every_row_arr[1]))
            one_arr.append(int(every_row_arr[2]))
            tar_bed_arr[14].append(one_arr)
        if every_row_arr[0] == 'chr16' or every_row_arr[0] == '16':
            one_arr = []
            one_arr.append(int(every_row_arr[1]))
            one_arr.append(int(every_row_arr[2]))
            tar_bed_arr[15].append(one_arr)
        if every_row_arr[0] == 'chr17' or every_row_arr[0] == '17':
            one_arr = []
            one_arr.append(int(every_row_arr[1]))
            one_arr.append(int(every_row_arr[2]))
            tar_bed_arr[16].append(one_arr)
        if every_row_arr[0] == 'chr18' or every_row_arr[0] == '18':
            one_arr = []
            one_arr.append(int(every_row_arr[1]))
            one_arr.append(int(every_row_arr[2]))
            tar_bed_arr[17].append(one_arr)
        if every_row_arr[0] == 'chr19' or every_row_arr[0] == '19':
            one_arr = []
            one_arr.append(int(every_row_arr[1]))
            one_arr.append(int(every_row_arr[2]))
            tar_bed_arr[18].append(one_arr)
        if every_row_arr[0] == 'chr20' or every_row_arr[0] == '20':
            one_arr = []
            one_arr.append(int(every_row_arr[1]))
            one_arr.append(int(every_row_arr[2]))
            tar_bed_arr[19].append(one_arr)
        if every_row_arr[0] == 'chr21' or every_row_arr[0] == '21':
            one_arr = []
            one_arr.append(int(every_row_arr[1]))
            one_arr.append(int(every_row_arr[2]))
            tar_bed_arr[20].append(one_arr)
        if every_row_arr[0] == 'chr22' or every_row_arr[0] == '22':
            one_arr = []
            one_arr.append(int(every_row_arr[1]))
            one_arr.append(int(every_row_arr[2]))
            tar_bed_arr[21].append(one_arr)
        if every_row_arr[0] == 'chrX' or every_row_arr[0] == 'X':
            one_arr = []
            one_arr.append(int(every_row_arr[1]))
            one_arr.append(int(every_row_arr[2]))
            tar_bed_arr[22].append(one_arr)
        if every_row_arr[0] == 'chrY' or every_row_arr[0] == 'Y':
            one_arr = []
            one_arr.append(int(every_row_arr[1]))
            one_arr.append(int(every_row_arr[2]))
            tar_bed_arr[23].append(one_arr)
if gradient == None:
    if mutation == "cnv":
        gradient = round(random.uniform(0, 10), 1)         #É³ɿ½±´ÊÌ¶È    
    else:
        gradient = round(random.random(), 4)
if per_gradient_sample_num == None:
    per_gradient_sample_num = 5
if per_sample_max_mutation_num == None:
    per_sample_max_mutation_num = 3
if outpath != None:
        out = outpath
        out = out.strip()
        out = out.rstrip("\/")
        isExists = os.path.exists(out)
        if not isExists:
            os.makedirs(out)
        out = out + '/'
else:
    pwd = os.getcwd()
    timestr = time.strftime('%Y%m%d%H%M%S', time.localtime(time.time()))
    if mutation == None:
        out = pwd + '/' + timestr + "normal"
    else:
        out = pwd + '/' + timestr + mutation + str(gradient)
    isExists = os.path.exists(out)
    if not isExists:
        os.makedirs(out)
    out = out + '/'


def makedict(tar_bed_arr):
    chr_dict = dict()
    for i in range(24):
        if tar_bed_arr[i] != []:
            if i == 22:
                chr = "chrX"
            elif i == 23:
                chr = "chrY"
            else:
                chr = "chr" + str(i+1)
            chr_dict[chr] = choice(tar_bed_arr[i])
    return chr_dict

def makelist(chr_dict):
    chr_list = []
    for chr in chr_dict.keys():
        chr_list.append(chr)
    return chr_list

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
#print chr_dict[chr_rnd][1]
chr_list = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY" ]

if mutation == "delete":
    write_tsv = open(out + "delete_" + str(gradient) + ".tsv", "ab")
    write_interfere_tsv = open(out + "interfere.tsv", "ab")
    write_tsv.write("sample_delete_" + str(gradient) + '\t' + "chromo" + '\t' + "start" + '\t' + "end" + '\t' + "if_insert" + '\t' + "if_calling" + '\n')
    delete_min = 1
    delete_max = 1000
    fragment_length_min = 1
    fragment_length_max = 1000
    write_shell = open(out + "tmp_simu_samples.sh", "ab")
    for per_gradient_sample in range(per_gradient_sample_num):
        fragment_length = random.randint(fragment_length_min, fragment_length_max)
        command = "python -B Wangshj.py -fra " + str(fragment_length) + " -delete "
        per_sample_mutation_num = random.randint(1, per_sample_max_mutation_num)
        #per_sample_mutation_num = 3100
        for per_sample_mutation in range(per_sample_mutation_num):
            if bed_file != None:
                chr_dict = makedict(tar_bed_arr)
                chr_list = makelist(chr_dict)
            chr_rnd = choice(chr_list)
            start_rnd = random.randint(chr_dict[chr_rnd][0], chr_dict[chr_rnd][1])
            delete_len_rnd = random.randint(delete_min - 1, delete_max - 1)
            end_rnd = min(start_rnd + delete_len_rnd, chr_dict[chr_rnd][1])
            if_insert = random.randint(0, 1)
            command = command + "onedelete," + chr_rnd + "," + str(start_rnd) + "," + str(end_rnd) + "," + str(gradient) + "," + str(if_insert) + ":"
            if if_insert == 1:
                to_write_tsv = "simu_sample" + str(per_gradient_sample + 1) + '\t' + chr_rnd + '\t' + str(start_rnd) + '\t' + str(end_rnd) + '\t' + str(fragment_length) + '\n'
            else:
                to_write_tsv = "simu_sample" + str(per_gradient_sample + 1) + '\t' + chr_rnd + '\t' + str(start_rnd) + '\t' + str(end_rnd) + '\t' + "0" + '\n'            
            write_tsv.write(to_write_tsv)   
            #print chr_rnd, start_rnd, end_rnd, fragment_length, if_insert
        this_out = out + "simu_sample" + str(per_gradient_sample + 1)
        os.makedirs(this_out)
        command_str, tsv_str = interfere_mutation("delete", per_gradient_sample)
        command = command[:-1] + command_str + " -out " + this_out
        write_interfere_tsv.write(tsv_str)
        if per_gradient_sample != per_gradient_sample_num - 1:
            command = command + "&"
        else:
            command = command + '\n'
        write_shell.write(command)
    delete_self = "rm " + '\"' + "$0" + '\"'
    write_shell.write(delete_self)
    write_shell.close()
    write_tsv.close()
    write_interfere_tsv.close()


if mutation == "tandem_repeat":
    write_tsv = open(out + "tandem_repeat_" + str(gradient) + ".tsv", "ab")
    write_interfere_tsv = open(out + "interfere.tsv", "ab")
    write_tsv.write("sample_tandem_repeat_" + str(gradient) + '\t' + "chromo" + '\t' + "start" + '\t' + "end" + '\t' + "copy" + '\t' + "if_calling" + '\n')
    tandem_repeat_min = 1
    tandem_repeat_max = 1000
    copy_min = 2
    copy_max = 10
    write_shell = open(out + "tmp_simu_samples.sh", "ab")
    for per_gradient_sample in range(per_gradient_sample_num):
        command = "python -B Wangshj.py -tandem_repeat "
        per_sample_mutation_num = random.randint(1, per_sample_max_mutation_num)
        for per_sample_mutation in range(per_sample_mutation_num):
            if bed_file != None:
                chr_dict = makedict(tar_bed_arr)
                chr_list = makelist(chr_dict)
            copy_num = random.randint(copy_min, copy_max)
            chr_rnd = choice(chr_list)
            start_rnd = random.randint(chr_dict[chr_rnd][0], chr_dict[chr_rnd][1])
            tandem_repeat_len_rnd = random.randint(tandem_repeat_min - 1, tandem_repeat_max - 1)
            end_rnd = min(start_rnd + tandem_repeat_len_rnd, chr_dict[chr_rnd][1])
            command = command + "onetandem_repeat," + chr_rnd + "," + str(start_rnd) + "," + str(end_rnd) + "," + str(copy_num) + "," + str(gradient) + ":"
            to_write_tsv = "simu_sample" + str(per_gradient_sample + 1) + '\t' + chr_rnd + '\t' + str(start_rnd) + '\t' + str(end_rnd) + '\t' + str(copy_num) + '\n'
            write_tsv.write(to_write_tsv)
        this_out = out + "simu_sample" + str(per_gradient_sample + 1)
        os.makedirs(this_out)
        command_str, tsv_str = interfere_mutation("tandem_repeat", per_gradient_sample)
        command = command[:-1] + command_str + " -out " + this_out
        write_interfere_tsv.write(tsv_str)
        if per_gradient_sample != per_gradient_sample_num - 1:
            command = command + "&"
        else:
            command = command + '\n'
        write_shell.write(command)
    delete_self = "rm " + '\"' + "$0" + '\"'
    write_shell.write(delete_self)
    write_shell.close()
    write_tsv.close()
    write_interfere_tsv.close()


if mutation == "insert":
    write_tsv = open(out + "insert_" + str(gradient) + ".tsv", "a")
    write_interfere_tsv = open(out + "interfere.tsv", "a")
    write_tsv.write("sample_insert_" + str(gradient) + '\t' + "chromo" + '\t' + "pos" + '\t' + "insert_len" + '\t' + "if_calling" + '\n')
    insert_len_min = 1
    insert_len_max = 1000
    write_shell = open(out + "tmp_simu_samples.sh", "a")
    for per_gradient_sample in range(per_gradient_sample_num):
        insert_len = random.randint(insert_len_min, insert_len_max)
        command = "python -B Wangshj.py -fra " + str(insert_len) + " -insert "
        per_sample_mutation_num = random.randint(1, per_sample_max_mutation_num)
        for per_sample_mutation in range(per_sample_mutation_num):
            if bed_file != None:
                chr_dict = makedict(tar_bed_arr)
                chr_list = makelist(chr_dict)
            chr_rnd = choice(chr_list)
            pos_rnd = random.randint(chr_dict[chr_rnd][0], chr_dict[chr_rnd][1])
            command = command + "oneinsert," + chr_rnd + "," + str(pos_rnd) + "," + str(gradient) + ":"
            to_write_tsv = "simu_sample" + str(per_gradient_sample + 1) + '\t' + chr_rnd + '\t' + str(pos_rnd) + '\t' + str(insert_len) + '\n'
            write_tsv.write(to_write_tsv)
        this_out = out + "simu_sample" + str(per_gradient_sample + 1)
        os.makedirs(this_out)
        command_str, tsv_str = interfere_mutation("insert", per_gradient_sample)
        command = command[:-1] + command_str + " -out " + this_out
        write_interfere_tsv.write(tsv_str)
        if per_gradient_sample != per_gradient_sample_num - 1:
            command = command + "&"
        else:
            command = command + '\n'
        write_shell.write(command)
    delete_self = "rm " + '\"' + "$0" + '\"'
    write_shell.write(delete_self)
    write_shell.close()
    write_tsv.close()
    write_interfere_tsv.close()


if mutation == "cnv":
    write_tsv = open(out + "cnv_" + str(gradient) + ".tsv", "ab")
    write_interfere_tsv = open(out + "interfere.tsv", "ab")
    write_tsv.write("sample_cnv_" + str(gradient) + '\t' + "chromo" + '\t' + "start" + '\t' + "end" + '\t' + "if_calling" + '\n')
    cnv_min = 1
    cnv_max = 1000
    write_shell = open(out + "tmp_simu_samples.sh", "ab")
    for per_gradient_sample in range(per_gradient_sample_num):
        command = "python -B Wangshj.py -cnv "
        per_sample_mutation_num = random.randint(1, per_sample_max_mutation_num)
        for per_sample_mutation in range(per_sample_mutation_num):
            if bed_file != None:
                chr_dict = makedict(tar_bed_arr)
                chr_list = makelist(chr_dict)
            chr_rnd = choice(chr_list)
            start_rnd = random.randint(chr_dict[chr_rnd][0], chr_dict[chr_rnd][1])
            cnv_len_rnd = random.randint(cnv_min - 1, cnv_max - 1)
            end_rnd = min(start_rnd + cnv_len_rnd, chr_dict[chr_rnd][1])
            command = command + "onecnv," + chr_rnd + "," + str(start_rnd) + "," + str(end_rnd) + "," + str(gradient) + ":"
            to_write_tsv = "simu_sample" + str(per_gradient_sample + 1) + '\t' + chr_rnd + '\t' + str(start_rnd) + '\t' + str(end_rnd) + '\n'
            write_tsv.write(to_write_tsv)
        this_out = out + "simu_sample" + str(per_gradient_sample + 1)
        os.makedirs(this_out)
        command_str, tsv_str = interfere_mutation("cnv", per_gradient_sample)
        command = command[:-1] + command_str + " -out " + this_out
        write_interfere_tsv.write(tsv_str)
        if per_gradient_sample != per_gradient_sample_num - 1:
            command = command + "&"
        else:
            command = command + '\n'
        write_shell.write(command)
    delete_self = "rm " + '\"' + "$0" + '\"'
    write_shell.write(delete_self)
    write_shell.close()
    write_tsv.close()
    write_interfere_tsv.close()


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
    #hg_fa_file = file("hg19.fa", 'r')
    hg_fa_file = open("hg19.fa", 'r')
    write_tsv = open(out + "snv_" + str(gradient) + ".tsv", "ab")
    write_interfere_tsv = open(out + "interfere.tsv", "ab")
    write_tsv.write("sample_snv_" + str(gradient) + '\t' + "chromo" + '\t' + "pos" + '\t' + "ref" + '\t' + "alt" + '\t' + "if_calling" + '\n')
    write_shell = open(out + "tmp_simu_samples.sh", "ab")
    for per_gradient_sample in range(per_gradient_sample_num):
        command = "python -B Wangshj.py -snv "
        per_sample_mutation_num = random.randint(1, per_sample_max_mutation_num)
        for per_sample_mutation in range(per_sample_mutation_num):
            if bed_file != None:
                chr_dict = makedict(tar_bed_arr)
                chr_list = makelist(chr_dict)
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
            command = command + "onesnv," + chr_rnd + "," + str(pos_rnd) + "," + ref + "," + alt + "," + str(gradient) + ":"
            to_write_tsv = "simu_sample" + str(per_gradient_sample + 1) + '\t' + chr_rnd + '\t' + str(pos_rnd) + '\t' + ref + '\t' + alt + '\n'
            write_tsv.write(to_write_tsv)
        this_out = out + "simu_sample" + str(per_gradient_sample + 1)
        os.makedirs(this_out)
        command_str, tsv_str = interfere_mutation("snv", per_gradient_sample)
        command = command[:-1] + command_str + " -out " + this_out
        write_interfere_tsv.write(tsv_str)
        if per_gradient_sample != per_gradient_sample_num - 1:
            command = command + "&"
        else:
            command = command + '\n'
        write_shell.write(command)
    delete_self = "rm " + '\"' + "$0" + '\"'
    write_shell.write(delete_self)
    write_shell.close()
    write_tsv.close()
    hg_fa_file.close()
    write_interfere_tsv.close()


if mutation == "inversion":
    write_tsv = open(out + "inversion_" + str(gradient) + ".tsv", "ab")
    write_interfere_tsv = open(out + "interfere.tsv", "ab")
    write_tsv.write("sample_inversion_" + str(gradient) + '\t' + "chromo" + '\t' + "start" + '\t' + "end" + '\t' + "if_calling" + '\n')
    inversion_min = 1
    inversion_max = 1000
    write_shell = open(out + "tmp_simu_samples.sh", "ab")
    for per_gradient_sample in range(per_gradient_sample_num):
        command = "python -B Wangshj.py -inversion "
        per_sample_mutation_num = random.randint(1, per_sample_max_mutation_num)
        for per_sample_mutation in range(per_sample_mutation_num):
            if bed_file != None:
                chr_dict = makedict(tar_bed_arr)
                chr_list = makelist(chr_dict)
            chr_rnd = choice(chr_list)
            start_rnd = random.randint(chr_dict[chr_rnd][0], chr_dict[chr_rnd][1])
            inversion_len_rnd = random.randint(inversion_min - 1, inversion_max - 1)
            end_rnd = min(start_rnd + inversion_len_rnd, chr_dict[chr_rnd][1])
            command = command + "oneinversion," + chr_rnd + "," + str(start_rnd) + "," + str(end_rnd) + "," + str(gradient) + ":"
            to_write_tsv = "simu_sample" + str(per_gradient_sample + 1) + '\t' + chr_rnd + '\t' + str(start_rnd) + '\t' + str(end_rnd) + '\n'
            write_tsv.write(to_write_tsv)
        this_out = out + "simu_sample" + str(per_gradient_sample + 1)
        os.makedirs(this_out)
        command_str, tsv_str = interfere_mutation("inversion", per_gradient_sample)
        command = command[:-1] + command_str + " -out " + this_out
        write_interfere_tsv.write(tsv_str)
        if per_gradient_sample != per_gradient_sample_num - 1:
            command = command + "&"
        else:
            command = command + '\n'
        write_shell.write(command)
    delete_self = "rm " + '\"' + "$0" + '\"'
    write_shell.write(delete_self)
    write_shell.close()
    write_tsv.close()
    write_interfere_tsv.close()


if mutation == "fusion":
    strand_list = ["+", "-"]
    write_tsv = open(out + "fusion_" + str(gradient) + ".tsv", "ab")
    write_interfere_tsv = open(out + "interfere.tsv", "ab")
    write_tsv.write("sample_fusion_" + str(gradient) + '\t' + "chromo" + '\t' + "strand" + '\t' + "pos" + '\t' + "chromo" + '\t' + "strand" + '\t' + "pos" + '\t' + "if_insert" + '\t' + "if_calling" + '\n')
    fragment_length_min = 1
    fragment_length_max = 1000
    write_shell = open(out + "tmp_simu_samples.sh", "ab")
    for per_gradient_sample in range(per_gradient_sample_num):
        fragment_length = random.randint(fragment_length_min, fragment_length_max)
        command = "python -B Wangshj.py -fra " + str(fragment_length) + " -fusion "
        per_sample_mutation_num = random.randint(1, per_sample_max_mutation_num)
        for per_sample_mutation in range(per_sample_mutation_num):
            if bed_file != None:
                chr_dict = makedict(tar_bed_arr)
                chr_list = makelist(chr_dict)
            l_chr_rnd = choice(chr_list)
            l_strand_rnd = choice(strand_list)
            l_pos_rnd = random.randint(chr_dict[l_chr_rnd][0], chr_dict[l_chr_rnd][1])
            r_chr_rnd = choice(chr_list)
            r_strand_rnd = choice(strand_list)
            r_pos_rnd = random.randint(chr_dict[r_chr_rnd][0], chr_dict[r_chr_rnd][1])
            if_insert = random.randint(0, 1)
            command = command + "onefusion," + l_chr_rnd + "," + l_strand_rnd + "," + str(l_pos_rnd) + "," + r_chr_rnd + "," + r_strand_rnd + "," + str(r_pos_rnd) + "," + str(gradient) + "," + str(if_insert) + ":"
            if if_insert == 1:
                to_write_tsv = "simu_sample" + str(per_gradient_sample + 1) + '\t' + l_chr_rnd + '\t' + l_strand_rnd + '\t' + str(l_pos_rnd) + '\t' + r_chr_rnd + '\t' + r_strand_rnd + '\t' + str(r_pos_rnd) + '\t' + str(fragment_length) + '\n'
            else:
                to_write_tsv = "simu_sample" + str(per_gradient_sample + 1) + '\t' + l_chr_rnd + '\t' + l_strand_rnd + '\t' + str(l_pos_rnd) + '\t' + r_chr_rnd + '\t' + r_strand_rnd + '\t' + str(r_pos_rnd) + '\t' + "0" + '\n'
            write_tsv.write(to_write_tsv)
        this_out = out + "simu_sample" + str(per_gradient_sample + 1)
        os.makedirs(this_out)
        command_str, tsv_str = interfere_mutation("fusion", per_gradient_sample)
        command = command[:-1] + command_str + " -out " + this_out
        write_interfere_tsv.write(tsv_str)
        if per_gradient_sample != per_gradient_sample_num - 1:
            command = command + "&"
        else:
            command = command + '\n'
        write_shell.write(command)
    delete_self = "rm " + '\"' + "$0" + '\"'
    write_shell.write(delete_self)
    write_shell.close()
    write_tsv.close()
    write_interfere_tsv.close()


if mutation == None:
    write_tsv = open(out + "normal.tsv", "ab")
    text = "sample_normal" + '\t' + "if_calling" + '\n'
    write_tsv.write(text.encode())
    write_shell = open(out + "tmp_simu_samples.sh", "ab")
    for per_gradient_sample in range(per_gradient_sample_num):
        if bed_file != None:
            command = "python -B Wangshj.py -bed_file " + bed_file + " "
        else: 
            command = "python -B Wangshj.py "
        to_write_tsv = "simu_sample" + str(per_gradient_sample + 1) + '\n'
        write_tsv.write(to_write_tsv.encode())
        this_out = out + "simu_sample" + str(per_gradient_sample + 1)
        os.makedirs(this_out)
        command = command[:-1] + " -out " + this_out
        if per_gradient_sample != per_gradient_sample_num - 1:
            command = command + "&"
        else:
            command = command + '\n'
        write_shell.write(command.encode())
    delete_self = "rm " + '\"' + "$0" + '\"'
    write_shell.write(delete_self.encode())
    write_shell.close()
    write_tsv.close()


os.system("sh " + out + "tmp_simu_samples.sh")
