# -*- coding: utf-8 -*-
"""
Created on Tue May  8 15:11:44 2018

@author: Administrator
"""
from collections import defaultdict
#import pysam
import re
def get_gene(file):
    gene = []
    f = open(file)
    line = f.readline()
    line = line.strip('\r')
    line = line.strip('\n')
    tmp = line.split('\t')

    while line:
        if tmp[3] not in gene:
            gene.append(tmp[3])
        line = f.readline()
        line = line.strip('\r')
        line = line.strip('\n')
        tmp = line.split('\t')

    return gene

class Gene_Element(object):
    def __init__(self, s, e, c):
        self.start = s
        self.end = e
        self.chr = c
def get_depth(source_file):
    s = defaultdict(lambda:defaultdict(list))
    f1 = open(source_file)
    line = f1.readline()
    while line:
        line = line.strip('\r')
        line = line.strip('\n')
        tmp = line.split('\t')
        line = f1.readline()
        if tmp[0] == 'inv':
            chr1 = tmp[1].split(':')[0]
            if 'chr' not in chr1:
                if len(chr1) < 3:
                    chr1 = 'chr' + chr1
            pos1 = tmp[1].split(':')[1]
            chr2 = tmp[2].split(':')[0]
            if 'chr' not in chr2:
                if len(chr2) < 3:
                    chr2 = 'chr' + chr2
            pos2 = tmp[2].split(':')[1]
            if int(pos1) > int(pos2):
                pos1, pos2 = pos2, pos1
            name = '{}_{}'.format(chr1,chr2)
            s['inv'][name].append( Element(int(pos1), int(pos2)) )
        elif tmp[0] == 'del':
            chr1 = tmp[1].split(':')[0]
            if 'chr' not in chr1:
                if len(chr1) < 3:
                    chr1 = 'chr' + chr1
            pos1 = tmp[1].split(':')[1]
            chr2 = tmp[2].split(':')[0]
            if 'chr' not in chr2:
                if len(chr2) < 3:
                    chr2 = 'chr' + chr2
            pos2 = tmp[2].split(':')[1]
            if int(pos1) > int(pos2):
                pos1, pos2 = pos2, pos1
            name = '{}_{}'.format(chr1,chr2)
            s['del'][name].append( Element(int(pos1), int(pos2)) )
        elif tmp[0] == 'dup':
            chr1 = tmp[1].split(':')[0]
            if 'chr' not in chr1:
                if len(chr1) < 3:
                    chr1 = 'chr' + chr1
            pos1 = tmp[1].split(':')[1]
            chr2 = tmp[2].split(':')[0]
            if 'chr' not in chr2:
                if len(chr2) < 3:
                    chr2 = 'chr' + chr2
            pos2 = tmp[2].split(':')[1]
            name = '{}_{}'.format(chr1,chr2)
            if int(pos1) > int(pos2):
                pos1, pos2 = pos2, pos1
            s['dup'][name].append( Element(int(pos1), int(pos2)) )
        elif tmp[0] == 'fusion1':
            chr1 = tmp[1].split(':')[0]
            if 'chr' not in chr1:
                if len(chr1) < 3:
                    chr1 = 'chr' + chr1
            pos1 = tmp[1].split(':')[1]
            chr2 = tmp[2].split(':')[0]
            if 'chr' not in chr2:
                if len(chr2) < 3:
                    chr2 = 'chr' + chr2
            pos2 = tmp[2].split(':')[1]
            name = '{}_{}'.format(chr1,chr2)
            s['fusion'][name].append( Element(int(pos1), int(pos2)) )
        elif tmp[0] == 'fusion2':
            chr1 = tmp[1].split(':')[0]
            if 'chr' not in chr1:
                if len(chr1) < 3:
                    chr1 = 'chr' + chr1
            pos1 = tmp[1].split(':')[1]
            chr2 = tmp[2].split(':')[0]
            if 'chr' not in chr2:
                if len(chr2) < 3:
                    chr2 = 'chr' + chr2
            pos2 = tmp[2].split(':')[1]
            name = '{}_{}'.format(chr1,chr2)
            s['fusion'][name].append( Element(int(pos1), int(pos2)) )
        elif tmp[0] == 'fusion3':
            chr1 = tmp[1].split(':')[0]
            if 'chr' not in chr1:
                if len(chr1) < 3:
                    chr1 = 'chr' + chr1
            pos1 = tmp[1].split(':')[1]
            chr2 = tmp[2].split(':')[0]
            if 'chr' not in chr2:
                if len(chr2) < 3:
                    chr2 = 'chr' + chr2
            pos2 = tmp[2].split(':')[1]
            name = '{}_{}'.format(chr1,chr2)
            s['fusion'][name].append( Element(int(pos1), int(pos2)) )
        elif tmp[0] == 'fusion4':
            chr1 = tmp[1].split(':')[0]
            if 'chr' not in chr1:
                if len(chr1) < 3:
                    chr1 = 'chr' + chr1
            pos1 = tmp[1].split(':')[1]
            chr2 = tmp[2].split(':')[0]
            if 'chr' not in chr2:
                if len(chr2) < 3:
                    chr2 = 'chr' + chr2
            pos2 = tmp[2].split(':')[1]
            name = '{}_{}'.format(chr1,chr2)
            s['fusion'][name].append( Element(int(pos1), int(pos2)) )
        else:
            continue
    f1.close()
    f2 = open('test.bed','a')
    for sv_type in s:
        for name in s[sv_type]:
            for e in s[sv_type][name]:
                chr1 = name.split('_')[0]
                chr1 = name.split('_')[1]
                f2.writelines('{}\t{}\t{}\n'.format(chr1,e.start-150,e.start+150))
                f2.writelines('{}\t{}\t{}\n'.format(chr2,e.end-150,e.end+150))
    f2.close()






def get_pos(gene, file = 'all_gene.bed', file_need = 'gene.bed' ):
    f = open(file)
    line = f.readline()
    line = line.strip('\r')
    line = line.strip('\n')
    tmp = line.split('\t')
    data = defaultdict(Gene_Element)
   
    line = f.readline()

    while(line):
        line = line.strip('\r')
        line = line.strip('\n')
        tmp = line.split('\t')
        line = f.readline()
        gene_name = tmp[3]
        chrom = tmp[0]
        start = tmp[1]
        end = tmp[2]

        if gene_name in gene:
            data[gene_name] = Gene_Element(start, end, chrom)
    f.close()
    f_n = open(file_need)
    line = f_n.readline()
    line = line.strip('\r')
    line = line.strip('\n')
    tmp = line.split('\t')
   
    line = f_n.readline()
    while(line):
        line = line.strip('\r')
        line = line.strip('\n')
        tmp = line.split('\t')
        line = f_n.readline()
        gene_name = tmp[3]
        chrom = tmp[0]
        if gene_name in data.keys():
            start = tmp[1] if int(tmp[1]) < int(data[gene_name].start) else data[gene_name].start
            end = tmp[2] if int(tmp[2]) > int(data[gene_name].end) else data[gene_name].end
            del data[gene_name]
            data[gene_name] = Gene_Element(start, end, chrom)
        else:
            start = tmp[1]
            end = tmp[2]
            data[gene_name] = Gene_Element(start, end, chrom)
    f_n.close()
    return data
def write_bed(data):
    f= open('gene.bed','a')
    for gene in data:
        info = data[gene]
        f.writelines('{}\t{}\t{}\n'.format(info.chr, info.start, info.end))
def form_bedpe_to_mutation(file, frequency, if_insert):
    inversion = ""
    delete = ""
    tandem_repeat = ""
    fusion = ""
    r = get_bedpe_src_data(file)
    for sv_type in r:
        for name in r[sv_type]:
            for e in r[sv_type][name]:
                chr1 = name.split('_')[0]
                chr2 = name.split('_')[1]
                pos1 = e.start
                pos2 = e.end
                if sv_type == 'inv':
                    result = 'inv,{},{},{},{},{}:'.format(chr1, pos1,pos2, frequency, if_insert)
                    result = result.strip('\r')
                    inversion = inversion + result
                elif sv_type == 'del':
                    result = 'del,{},{},{},{},{}:'.format(chr1, pos1,pos2, frequency, if_insert)
                    result = result.strip('\r')
                    delete = delete + result
                elif sv_type == 'dup':
                    result = 'dup,{},{},{},{},{}:'.format(chr1, pos1,pos2, 2, frequency)
                    result = result.strip('\r')
                    tandem_repeat = tandem_repeat + result
                elif sv_type == 'fusion':
                            result = 'fusion1,{},+,{},{},-,{},{},{}:'.format(chr1, pos1,chr2,pos2,frequency, if_insert)
                            result = result.strip('\r')
                            fusion = fusion + result
    if inversion == '' :
        inversion = None
    else:
        inversion = inversion[:-1]

    if delete == '' :
        delete = None
    else:
        delete = delete[:-1]

    if tandem_repeat == '':
        tandem_repeat = None
    else:
        tandem_repeat = tandem_repeat[:-1]

    if fusion == '':
        fusion = None
    else:
        fusion = fusion[:-1]

    return inversion, delete, tandem_repeat, fusion
def get_mutation(file, frequency='0.1', if_insert='0'):
    if 'bedpe' in file:
        inversion, delete, tandem_repeat, fusion = form_bedpe_to_mutation(file,frequency,if_insert)
    elif 'txt' in file:
        inversion, delete, tandem_repeat, fusion = form_txt_to_mutation(file)
    return inversion, delete, tandem_repeat, fusion
def form_txt_to_mutation(file):
    inversion = ""
    delete = ""
    tandem_repeat = ""
    fusion = ""
    f = open(file)
    line = f.readline()
    line = line.strip('\r')
    line = line.strip('\n')
    tmp = line.split('\t')
    line = f.readline()
    print(tmp)
    if 'sv_type' in tmp:
        sv_index = tmp.index('sv_type')
    else:
        print('do not have sv_type')
    if 'pos1' in tmp:
        pos1_index = tmp.index('pos1')
    else:
        print('do not have pos1')
    if 'pos2' in tmp:
        pos2_index = tmp.index('pos2')
    else:
        print('do not have person data')

    if 'if_insert' in tmp:
        insert_index = tmp.index('if_insert')
    else:
        print('do not have insert info')
    if 'frequency' in tmp:
        frequency_index = tmp.index('frequency')
    else:
        print('do not have frequency info')


    while line:
        line = line.strip('\r')
        line = line.strip('\n')
        tmp = line.split('\t')
        line = f.readline()
        sv_type = tmp[sv_index]
        pos1_info = tmp[pos1_index]
        pos2_info = tmp[pos2_index]
        if_insert = tmp[insert_index]
        frequency = tmp[frequency_index]
        if sv_type == 'inv':
            chr1 = pos1_info.split(':')[0]
            if 'chr' not in chr1:
                if len(chr1) < 3:
                    chr1 = 'chr' + chr1
                else:
                    continue
            pos1 = pos1_info.split(':')[1]
            #chr2 = tmp[2].split(':')[0]
            pos2 = pos2_info.split(':')[1]
            if int(pos1) > int(pos2):
                pos1, pos2 = pos2, pos1
            result = 'inv,{},{},{},{},{}:'.format(chr1, pos1,pos2, frequency, if_insert)
            result = result.strip('\r')
            inversion = inversion + result
        elif sv_type == 'del':
            chr1 = pos1_info.split(':')[0]
            if 'chr' not in chr1:
                    if len(chr1) < 3:
                        chr1 = 'chr' + chr1
                    else:
                        continue
            pos1 = pos1_info.split(':')[1]
            #chr2 = tmp[2].split(':')[0]
            pos2 = pos2_info.split(':')[1]
            if int(pos1) > int(pos2):
                pos1, pos2 = pos2, pos1
            result = 'del,{},{},{},{},{}:'.format(chr1, pos1,pos2, frequency, if_insert)
            result = result.strip('\r')
            delete = delete + result
        elif sv_type == 'dup':
            chr1 = pos1_info.split(':')[0]
            if 'chr' not in chr1:
                if len(chr1) < 3:
                    chr1 = 'chr' + chr1
                else:
                    continue
            pos1 = pos1_info.split(':')[1]
            #chr2 = tmp[2].split(':')[0]
            pos2 = pos2_info.split(':')[1]
            if int(pos1) > int(pos2):
                pos1, pos2 = pos2, pos1
            result = 'dup,{},{},{},{},{}:'.format(chr1, pos1,pos2, 2, frequency)
            result = result.strip('\r')
            tandem_repeat = tandem_repeat + result
        elif sv_type == 'fusion1':
            chr1 = pos1_info.split(':')[0]
            if 'chr' not in chr1:
                if len(chr1) < 3:
                    chr1 = 'chr' + chr1
                else:
                    continue
            pos1 = pos1_info.split(':')[1]
            chr2 = pos2_info.split(':')[0]
            if 'chr' not in chr2:
                if len(chr2) < 3:
                    chr2 = 'chr' + chr2
                else:
                    continue
            pos2 = pos2_info.split(':')[1]
            result = 'fusion1,{},+,{},{},-,{},{},{}:'.format(chr1, pos1,chr2,pos2,frequency, if_insert)
            result = result.strip('\r')
            fusion = fusion + result
        elif sv_type == 'fusion2':
            chr1 = pos1_info.split(':')[0]
            if 'chr' not in chr1:
                if len(chr1) < 3:
                    chr1 = 'chr' + chr1
                else:
                    continue
            pos1 = pos1_info.split(':')[1]
            chr2 = pos2_info.split(':')[0]
            if 'chr' not in chr2:
                if len(chr2) < 3:
                    chr2 = 'chr' + chr2
                else:
                    continue
            pos2 = pos2_info.split(':')[1]
            result = 'fusion2,{},-,{},+,{},{},{},{}:'.format(chr1, pos1,chr2,pos2,frequency, if_insert)
            fusion = fusion + result
        elif sv_type == 'fusion3':
            chr1 = pos1_info.split(':')[0]
            if 'chr' not in chr1:
                if len(chr1) < 3:
                    chr1 = 'chr' + chr1
                else:
                    continue
            pos1 = pos1_info.split(':')[1]
            chr2 = pos2_info.split(':')[0]
            if 'chr' not in chr2:
                if len(chr2) < 3:
                    chr2 = 'chr' + chr2
                else:
                    continue
            pos2 = pos2_info.split(':')[1]
            result = 'fusion3,{},+,{},{},+,{},{},{}:'.format(chr1, pos1,chr2,pos2,frequency, if_insert)
            result = result.strip('\r')
            fusion = fusion + result
        elif sv_type == 'fusion4':
            chr1 = pos1_info.split(':')[0]
            if 'chr' not in chr1:
                if len(chr1) < 3:
                    chr1 = 'chr' + chr1
                else:
                    continue
            pos1 = pos1_info.split(':')[1]
            chr2 =pos2_info.split(':')[0]
            if 'chr' not in chr2:
                if len(chr2) < 3:
                    chr2 = 'chr' + chr2
                else:
                    line = f.readline()
                    continue
            pos2 = pos2_info.split(':')[1]
            result = 'fusion4,{},+,{},{},+,{},{},{}:'.format(chr1, pos1,chr2,pos2,frequency,if_insert)
            result = result.strip('\r')
            fusion = fusion + result
        else:
            continue
    f.close()
    if inversion == '' :
        inversion = None
    else:
        inversion = inversion[:-1]

    if delete == '' :
        delete = None
    else:
        delete = delete[:-1]

    if tandem_repeat == '':
        tandem_repeat = None
    else:
        tandem_repeat = tandem_repeat[:-1]

    if fusion == '':
        fusion = None
    else:
        fusion = fusion[:-1]

    return inversion, delete, tandem_repeat, fusion

class Config(object):
    LEFT_INVERSION = 0
    RIGHT_INVERSION = 1
    DELETION = 2
    DUPLICATION = 3
    '''
    TRANSLOCATION_0:in order inversion 0
    TRANSLOCATION_1:reverse order inversion 1
    TRANSLOCATION_2:in order not inversion 2
    TRANSLOCATION_3:in order inversion 3
    '''
    FUSION_1 = 5
    FUSION_2 = 6
    FUSION_3 = 7
    FUSION_4 = 8

class Element(object):
    def __init__(self, s, e, p=''):
        self.start = s
        self.end = e
        self.hit = False
        self.if_precise = p
def get_vcf_res_data(source_file='output.vcf'):
    r = defaultdict(lambda:defaultdict(list))
    for c in pysam.VariantFile(source_file):
            info = c.info
            sv_type = info['SVTYPE']
            if 'DEL' in sv_type:
                if info['IMPRECISE']:
                    precise = 'imprecise'
                else:
                    precise = 'precise'
                chr1 = c.chrom
                chr2 = c.contig
                if 'chr' not in chr1:
                    chr1 = 'chr{}'.format(chr1)
                if 'chr' not in chr2:
                    chr2 = 'chr{}'.format(chr2)
                name = '{}_{}'.format(chr1,chr2)
                pos1 = c.start
                pos2 = c.stop
                r['del'][name].append(Element(int(pos1), int(pos2), precise) )
            elif 'INV' in sv_type:
                if info['IMPRECISE']:
                    precise = 'imprecise'
                else:
                    precise = 'precise'
                chr1 = c.chrom
                chr2 = c.contig
                if 'chr' not in chr1:
                    chr1 = 'chr{}'.format(chr1)
                if 'chr' not in chr2:
                    chr2 = 'chr{}'.format(chr2)
                name = '{}_{}'.format(chr1,chr2)
                pos1 = c.start
                pos2 = c.stop
                r['inv'][name].append(Element(int(pos1), int(pos2), precise) )
            elif 'DUP' in sv_type:
                if info['IMPRECISE']:
                    precise = 'imprecise'
                else:
                    precise = 'precise'
                chr1 = c.chrom
                chr2 = c.contig
                if 'chr' not in chr1:
                    chr1 = 'chr{}'.format(chr1)
                if 'chr' not in chr2:
                    chr2 = 'chr{}'.format(chr2)
                name = '{}_{}'.format(chr1,chr2)
                pos1 = c.start
                pos2 = c.stop
                r['dup'][name].append(Element(int(pos1), int(pos2), precise) )
            elif 'INS' in sv_type:
                if info['IMPRECISE']:
                    precise = 'imprecise'
                else:
                    precise = 'precise'
                chr1 = c.chrom
                chr2 = c.contig
                if 'chr' not in chr1:
                    chr1 = 'chr{}'.format(chr1)
                if 'chr' not in chr2:
                    chr2 = 'chr{}'.format(chr2)
                name = '{}_{}'.format(chr1,chr2)
                pos1 = c.start
                pos2 = c.stop
                r['ins'][name].append(Element(int(pos1), int(pos2), precise) )
            elif 'BND' in sv_type:
                if info['IMPRECISE']:
                    precise = 'imprecise'
                else:
                    precise = 'precise'
                chr1 = c.chrom
                a = c.alts[0].split(':')[0]
                chr2 = re.sub('[^0-9]','', a)
                if 'chr' not in chr1:
                    chr1 = 'chr{}'.format(chr1)
                if 'chr' not in chr2:
                    chr2 = 'chr{}'.format(chr2)
                chr1, chr2 = chr2, chr1
                name = '{}_{}'.format(chr1,chr2)
                pos1 = c.start
                a = c.alts[0].split(':')[1]
                pos2 = re.sub('[^0-9]','', a)
                pos1, pos2 = pos2, pos1
                r['fusion'][name].append(Element(int(pos1), int(pos2), precise) )
    return r


def get_txt_res_data(result_file='res.txt'):
    f = open(result_file)
    r = defaultdict(lambda:defaultdict(list))
    line = f.readline()
    while line:
        line = line.strip('\r')
        line = line.strip('\n')
        tmp = line.split('\t')
        line = f.readline()
        if_precise = tmp[0]
        if if_precise != 'precise' and if_precise != 'imprecise':
            continue
        tmp = tmp[1:]
        if int(tmp[0]) == Config.LEFT_INVERSION or tmp[0] == Config.RIGHT_INVERSION:
            chr1 = tmp[1].split(':')[0]
            if 'chr' not in chr1:
                if len(chr1) < 3:
                    chr1 = 'chr' + chr1
            pos1 = tmp[1].split(':')[1]
            chr2 = tmp[2].split(':')[0]
            if 'chr' not in chr2:
                if len(chr2) < 3:
                    chr2 = 'chr' + chr2
            pos2 = tmp[2].split(':')[1]
            if int(pos1) > int(pos2):
                pos1, pos2 = pos2, pos1
            name = '{}_{}'.format(chr1,chr2)
            r['inv'][name].append(Element(int(pos1), int(pos2), if_precise) )
        elif int(tmp[0]) == Config.DELETION:
            chr1 = tmp[1].split(':')[0]
            if 'chr' not in chr1:
                if len(chr1) < 3:
                    chr1 = 'chr' + chr1
            pos1 = tmp[1].split(':')[1]
            chr2 = tmp[2].split(':')[0]
            if 'chr' not in chr2:
                if len(chr2) < 3:
                    chr2 = 'chr' + chr2

            pos2 = tmp[2].split(':')[1]
            if int(pos1) > int(pos2):
                pos1, pos2 = pos2, pos1
            name = '{}_{}'.format(chr1,chr2)
            r['del'][name].append( Element(int(pos1), int(pos2), if_precise) )
        elif int(tmp[0]) == Config.DUPLICATION:
            chr1 = tmp[1].split(':')[0]
            if 'chr' not in chr1:
                if len(chr1) < 3:
                    chr1 = 'chr' + chr1
            pos1 = tmp[1].split(':')[1]
            chr2 = tmp[2].split(':')[0]
            if 'chr' not in chr2:
                if len(chr2) < 3:
                    chr2 = 'chr' + chr2
            pos2 = tmp[2].split(':')[1]
            name = '{}_{}'.format(chr1,chr2)
            if int(pos1) > int(pos2):
                pos1, pos2 = pos2, pos1
            r['dup'][name].append( Element(int(pos1), int(pos2), if_precise) )
        elif int(tmp[0]) == Config.FUSION_1:

            chr1 = tmp[1].split(':')[0]
            if 'chr' not in chr1:
                if len(chr1) < 3:
                    chr1 = 'chr' + chr1
            pos1 = tmp[1].split(':')[1]
            chr2 = tmp[2].split(':')[0]
            if 'chr' not in chr2:
                if len(chr2) < 3:
                    chr2 = 'chr' + chr2
            pos2 = tmp[2].split(':')[1]
            name = '{}_{}'.format(chr2,chr1)
            r['fusion'][name].append ( Element(int(pos2), int(pos1), if_precise) )
        elif int(tmp[0]) == Config.FUSION_2:
            chr1 = tmp[1].split(':')[0]
            if 'chr' not in chr1:
                if len(chr1) < 3:
                    chr1 = 'chr' + chr1
            pos1 = tmp[1].split(':')[1]
            chr2 = tmp[2].split(':')[0]
            if 'chr' not in chr2:
                if len(chr2) < 3:
                    chr2 = 'chr' + chr2
            pos2 = tmp[2].split(':')[1]
            name = '{}_{}'.format(chr2,chr1)
            r['fusion'][name].append ( Element(int(pos2), int(pos1), if_precise) )
        elif int(tmp[0]) == Config.FUSION_3:
            chr1 = tmp[1].split(':')[0]
            if 'chr' not in chr1:
                if len(chr1) < 3:
                    chr1 = 'chr' + chr1
            pos1 = tmp[1].split(':')[1]
            chr2 = tmp[2].split(':')[0]
            if 'chr' not in chr2:
                if len(chr2) < 3:
                    chr2 = 'chr' + chr2
            pos2 = tmp[2].split(':')[1]
            name = '{}_{}'.format(chr2,chr1)
            r['fusion'][name].append ( Element(int(pos2), int(pos1), if_precise) )
        elif int(tmp[0]) == Config.FUSION_4:
            chr1 = tmp[1].split(':')[0]
            if 'chr' not in chr1:
                if len(chr1) < 3:
                    chr1 = 'chr' + chr1
            pos1 = tmp[1].split(':')[1]
            chr2 = tmp[2].split(':')[0]
            if 'chr' not in chr2:
                if len(chr2) < 3:
                    chr2 = 'chr' + chr2
            pos2 = tmp[2].split(':')[1]
            name = '{}_{}'.format(chr2,chr1)
            r['fusion'][name].append ( Element(int(pos2), int(pos1), if_precise) )
        else:
            continue
    f.close()
    return r
def get_bedpe_src_data(source_file):
    s = defaultdict(lambda:defaultdict(list))
    f1 = open(source_file)
    line = f1.readline()
    count = 0
    while line:
        count += 1
        line = line.strip('\r')
        line = line.strip('\n')
        tmp = line.split('\t')
        line = f1.readline()
        if 'DUP' in tmp[6]:
            chr1 = tmp[0]
            if 'chr' not in chr1:
                if len(chr1) < 3:
                    chr1 = 'chr' + chr1
            chr2 = chr1
            pos1 = int(tmp[1])
            pos2 = int(tmp[4])
        if pos1 > pos2:
            pos1,pos2 = pos2,pos1
            name = '{}_{}'.format(chr1,chr2)
            s['dup'][name].append( Element(int(pos1), int(pos2)) )
        elif 'DEL' in tmp[6]:
            chr1 = tmp[0]
            if 'chr' not in chr1:
                if len(chr1) < 3:
                    chr1 = 'chr' + chr1
            chr2 = chr1
            pos1 = int(tmp[1])
            pos2 = int(tmp[4])
        if pos1 > pos2:
            pos1,pos2 = pos2,pos1
            name = '{}_{}'.format(chr1,chr2)
            s['del'][name].append( Element(int(pos1), int(pos2)) )
        elif 'INV' in tmp[6]:
            chr1 = tmp[0]
            if 'chr' not in chr1:
                if len(chr1) < 3:
                    chr1 = 'chr' + chr1
            chr2 = chr1
            pos1 = int(tmp[1])
            pos2 = int(tmp[4])
        if pos1 > pos2:
            pos1,pos2 = pos2,pos1
            name = '{}_{}'.format(chr1,chr2)
            s['inv'][name].append( Element(int(pos1), int(pos2)) )
        elif 'INR' in tmp[6]:
            chr1 = tmp[0]
            chr2 = tmp[3]
            pos1 = int(tmp[1])
            pos2 = int(tmp[4])

        count += 1

        if chr1.isdigit() and chr2.isdigit():
            if int(chr1) > int(chr2):
                chr1,chr2 = chr2,chr1
                pos1,pos2 = pos2,pos1
                chr1 = 'chr' + chr1
                chr2 = 'chr' + chr2
            else:
                if chr1 > chr2:
                    chr1,chr2 = chr2,chr1
                    pos1,pos2 = pos2,pos1
                if chr1.isdigit() :
                    chr1 = 'chr' + chr1
                if chr2.isdigit():
                    chr2 = 'chr' + chr2
            name = '{}_{}'.format(chr1,chr2)
            s['fusion'][name].append( Element(int(pos1), int(pos2)) )

    return s
def  get_txt_src_data(source_file):
    s = defaultdict(lambda:defaultdict(list))
    f1 = open(source_file)
    line = f1.readline()
    while line:
        line = line.strip('\r')
        line = line.strip('\n')
        tmp = line.split('\t')
        line = f1.readline()
        if tmp[0] == 'inv':
            chr1 = tmp[1].split(':')[0]
            if 'chr' not in chr1:
                if len(chr1) < 3:
                    chr1 = 'chr' + chr1
            pos1 = tmp[1].split(':')[1]
            chr2 = tmp[2].split(':')[0]
            if 'chr' not in chr2:
                if len(chr2) < 3:
                    chr2 = 'chr' + chr2
            pos2 = tmp[2].split(':')[1]
            if int(pos1) > int(pos2):
                pos1, pos2 = pos2, pos1
            name = '{}_{}'.format(chr1,chr2)
            s['inv'][name].append( Element(int(pos1), int(pos2)) )
        elif tmp[0] == 'del':
            chr1 = tmp[1].split(':')[0]
            if 'chr' not in chr1:
                if len(chr1) < 3:
                    chr1 = 'chr' + chr1
            pos1 = tmp[1].split(':')[1]
            chr2 = tmp[2].split(':')[0]
            if 'chr' not in chr2:
                if len(chr2) < 3:
                    chr2 = 'chr' + chr2
            pos2 = tmp[2].split(':')[1]
            if int(pos1) > int(pos2):
                pos1, pos2 = pos2, pos1
            name = '{}_{}'.format(chr1,chr2)
            s['del'][name].append( Element(int(pos1), int(pos2)) )
        elif tmp[0] == 'dup':
            chr1 = tmp[1].split(':')[0]
            if 'chr' not in chr1:
                if len(chr1) < 3:
                    chr1 = 'chr' + chr1
            pos1 = tmp[1].split(':')[1]
            chr2 = tmp[2].split(':')[0]
            if 'chr' not in chr2:
                if len(chr2) < 3:
                    chr2 = 'chr' + chr2
            pos2 = tmp[2].split(':')[1]
            name = '{}_{}'.format(chr1,chr2)
            if int(pos1) > int(pos2):
                pos1, pos2 = pos2, pos1
            s['dup'][name].append( Element(int(pos1), int(pos2)) )
        elif tmp[0] == 'fusion1':
            chr1 = tmp[1].split(':')[0]
            if 'chr' not in chr1:
                if len(chr1) < 3:
                    chr1 = 'chr' + chr1
            pos1 = tmp[1].split(':')[1]
            chr2 = tmp[2].split(':')[0]
            if 'chr' not in chr2:
                if len(chr2) < 3:
                    chr2 = 'chr' + chr2
            pos2 = tmp[2].split(':')[1]
            name = '{}_{}'.format(chr1,chr2)
            s['fusion'][name].append( Element(int(pos1), int(pos2)) )
        elif tmp[0] == 'fusion2':
            chr1 = tmp[1].split(':')[0]
            if 'chr' not in chr1:
                if len(chr1) < 3:
                    chr1 = 'chr' + chr1
            pos1 = tmp[1].split(':')[1]
            chr2 = tmp[2].split(':')[0]
            if 'chr' not in chr2:
                if len(chr2) < 3:
                    chr2 = 'chr' + chr2
            pos2 = tmp[2].split(':')[1]
            name = '{}_{}'.format(chr1,chr2)
            s['fusion'][name].append( Element(int(pos1), int(pos2)) )
        elif tmp[0] == 'fusion3':
            chr1 = tmp[1].split(':')[0]
            if 'chr' not in chr1:
                if len(chr1) < 3:
                    chr1 = 'chr' + chr1
            pos1 = tmp[1].split(':')[1]
            chr2 = tmp[2].split(':')[0]
            if 'chr' not in chr2:
                if len(chr2) < 3:
                    chr2 = 'chr' + chr2
            pos2 = tmp[2].split(':')[1]
            name = '{}_{}'.format(chr1,chr2)
            s['fusion'][name].append( Element(int(pos1), int(pos2)) )
        elif tmp[0] == 'fusion4':
            chr1 = tmp[1].split(':')[0]
            if 'chr' not in chr1:
                if len(chr1) < 3:
                    chr1 = 'chr' + chr1
            pos1 = tmp[1].split(':')[1]
            chr2 = tmp[2].split(':')[0]
            if 'chr' not in chr2:
                if len(chr2) < 3:
                    chr2 = 'chr' + chr2
            pos2 = tmp[2].split(':')[1]
            name = '{}_{}'.format(chr1,chr2)
            s['fusion'][name].append( Element(int(pos1), int(pos2)) )
        else:
            continue
    f1.close()

    return s

def get_data(source_file='del.txt', result_file='res.txt'):
    s = defaultdict(lambda:defaultdict(list))
    r = defaultdict(lambda:defaultdict(list))

    if 'txt' in result_file:
        r = get_txt_res_data(result_file)
    elif 'vcf' in result_file or 'bcf' in  result_file:
        r = get_vcf_res_data(result_file)
   
    if 'txt' in source_file:
        s = get_txt_src_data(source_file)
    elif 'bedpe' in source_file:
        s = get_bedpe_src_data(source_file)



    for sv_type in s:
        for name in s[sv_type]:
            s[sv_type][name].sort(key = (lambda e: e.end))
            s[sv_type][name].sort(key = (lambda e: e.start))

    for sv_type in r:
        for name in r[sv_type]:
            r[sv_type][name].sort(key = (lambda e: e.end))
            r[sv_type][name].sort(key = (lambda e: e.start))
    '''       
    for sv_type in s:
        for name in s[sv_type]:
            for e in s[sv_type][name]:
                f_w.writelines('sv_type:{},{},start:{},end:{}\n'.format(sv_type, name, e.start, e.end))
    
    for sv_type in r:
        for name in r[sv_type]:
            for e in r[sv_type][name]:
                f_w.writelines('sv_type:{},{},start:{},end:{}\n'.format(sv_type, name, e.start, e.end))
    '''


    for sv_type in s:
        if sv_type not in r.keys():
            continue
        for name in s[sv_type]:
            if name not in r[sv_type].keys():
                continue
            for e1 in r[sv_type][name]:
                for e2 in s[sv_type][name]:
                    if e1.if_precise == 'precise':
                        if e1.start in range (e2.start - 10, e2.start + 10) and e1.end in range(e2.end - 10 , e2.end + 10):
                            e1.hit = True
                            e2.hit = True
                            e2.if_precise = 'precise'
                    elif e1.if_precise == 'imprecise':
                        if e1.start in range (e2.start - 500, e2.start + 500) and e1.end in range(e2.end - 500 , e2.end + 500):
                            e1.hit = True
                            e2.hit = True
                            e2.if_precise = 'imprecise'

    fdr_count = 0
    for sv_type in r:
        for name in r[sv_type]:
            for e in r[sv_type][name]:
                if e.hit == True:
                    continue
                else:
                    fdr_count += 1
                    print('sv_type:{},{},pos1:{},pos2:{}\n'.format(sv_type, name, e.start, e.end))
    pre_hit_count = 0
    impre_hit_count = 0
    for sv_type in s:
        for name in s[sv_type]:
            for e in s[sv_type][name]:
                if e.hit == True:
                    if e.if_precise == 'precise':
                        pre_hit_count += 1
                    elif e.if_precise == 'imprecise':
                        impre_hit_count += 1


    print('precise hit count:{},imprecise hit count:{},fdr count:{}'.format(pre_hit_count, impre_hit_count, fdr_count))
   
    #f_w.close()




