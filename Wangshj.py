#!C:\Python27\python.exe
# -*- coding: UTF-8 -*
from random import choice
from mutation import get_mutation
import argparse
from get_depth_loop import create_loop
from insert import insert_imitate
from inversion import inversion_imitate
from insertmaking import insert_making
import os
import multiprocessing
import globalname
import random
import math
from loadjsondata import load_json_file
from cnv import cnv_multiple
from targetregioncapture import loadtargetregioncapture
from capturetochromo import capture_to_chromo
#from readfasta import load_hg19_data
from initialsequence import load_initial_sequence
from snv import snv_imitate
from fuse import fuse_imitate
from readsProducing import reads_making
from complementary import complementary
import time
import sys
from templet_fit_random import make_templet_length
from bamtocnv import graBam
import numpy as np
import os
import re



parser = argparse.ArgumentParser()
parser.add_argument('-fre', '--frequence_type', help='the frequence_type, you can set as AF EAS_AF AMR_AF AFR_AF EUR_AF or SAS_AF', required=False, type=str)
parser.add_argument('-dep_fit', '--depth_pattern_fitting', help='the depth_pattern_fitting', required=False, type=int)
parser.add_argument('-fra', '--fragment_length', help='the fragment_length', required=False, type=int)
parser.add_argument('-fib', '--fragment_insert_base', help='the fragment_insert_base', required=False, type=str)
parser.add_argument('-depth', '--depth', help='the depth', required=False, type=float)
parser.add_argument('-pe', '--pair_end', help='the pair_end', required=False, type=int)
parser.add_argument('-readlen', '--readlen', help='the readlen', required=False, type=int)
#parser.add_argument('-genome', '--genome', help='the genome', required=False, type=str)
parser.add_argument('-err', '--seq_error_rate', help='the seq_error_rate', required=False, type=float)
parser.add_argument('-dup', '--duplication_rate', help='the duplication_rate', required=False, type=float)
parser.add_argument('-ind1', '--random_index1', help='the random_index1', required=False, type=int)
parser.add_argument('-ind2', '--random_index2', help='the random_index2', required=False, type=int)
parser.add_argument('-sc', '--sqing_cycles', help='the sqing_cycles', required=False, type=int)
parser.add_argument('-ada1', '--read1_adapter', help='the read1_adapter', required=False, type=str)
parser.add_argument('-ada2', '--read2_adapter', help='the read2_adapter', required=False, type=str)
parser.add_argument('-tem_mi', '--template_len_min', help='the template_len_min', required=False, type=int)
parser.add_argument('-tem_ma', '--template_len_max', help='the template_len_max', required=False, type=int)
parser.add_argument('-nor_mi', '--normal_base_qual_min', help='the normal_base_qual_min', required=False, type=int)
parser.add_argument('-nor_ma', '--normal_base_qual_max', help='the norm	al_base_qual_max', required=False, type=int)
parser.add_argument('-err_mi', '--seq_error_qual_min', help='the seq_error_qual_min', required=False, type=int)
parser.add_argument('-err_ma', '--seq_error_qual_max', help='the seq_error_qual_max', required=False, type=int)
parser.add_argument('-fusion', '--fusion', help='the fusion', required=False, type=str)
parser.add_argument('-snv', '--snv', help='the snv, example:-snv onesnv,chr12,25398284,C,T,0.001', required=False, type=str)
parser.add_argument('-cnv', '--cnv', help='the cnv', required=False, type=str)
parser.add_argument('-inversion', '--inversion', help='the inversion', required=False, type=str)
parser.add_argument('-delete', '--delete', help='the delete', required=False, type=str)
parser.add_argument('-tandem_repeat', '--tandem_repeat', help='the tandem_repeat', required=False, type=str)
parser.add_argument('-insert', '--insert', help='the insert', required=False, type=str)
parser.add_argument('-sv_file', '--sv_file', help='the sv_file', required=False, type=str)
parser.add_argument('-bed_file', '--bed_file', help='the bed_file, example:-bed_file targetregioncapture.bed', required=False, type=str)
parser.add_argument('-p', '--multiprocess', help='the multiprocess', required=False, type=int)
parser.add_argument('-out', '--outpath', help='the outpath, example:-out /mnt/X500/farmers/wangshj/simu_samples/SimulationSample025', required=False, type=str)
parser.add_argument('-bamtocnv', '--bamtocnv', help='example:-bamtocnv /mnt/NL200_1/prod/workspace/IFA20170803002/OncoH/output/170005964BD/normal/5_recal_bam/170005964BD_normal_sort_markdup_realign_recal.bam,0.5:0.6,1:2,100000000:110000000,100300000:110300000', required=False, type=str)
parser.add_argument('-contiguous_point_depth_model', '--contiguous_point_depth_model', help='the contiguous_point_depth_model', required=False, type=str)
parser.add_argument('-template_length_distribution_model', '--templateLengthDistributionModel_path', help='the templateLengthDistributionModel_path', required=False, type=str)
parser.add_argument('-base_quality_distribution_model', '--baseQualityDistributionModel_path', help='the baseQualityDistributionModel_path', required=False, type=str)
parser.add_argument('-overall_depth_distribution_model', '--overallDepthDistributionModel_path', help='the overallDepthDistributionModel_path', required=False, type=str)


def multi_process_write(appro_templet_amount, n, insert_fragment, rdm, customized_dict, hg19index_dict, target_region_size, target_region_list, templetmin, templetmax, contiguous_point_depth_densityorparameter, templateLengthDistributionModel_path,out,depth_dict):
	one = 2 * n - 1
	another = 2 * n
	globalname.fastq1 = open(out + "makeit" + str(one) + ".fq", "a")
	globalname.fastq2 = open(out + "makeit" + str(another) + ".fq", "a")
	globalname.insert_fragment = insert_fragment
	globalname.rdm = rdm
	globalname.frequence_type = customized_dict["frequence_type"]
	appro_templet_amount = int(appro_templet_amount)
	templet_length_array = make_templet_length(appro_templet_amount, templetmax, templetmin, templateLengthDistributionModel_path)
	for num in range(appro_templet_amount):
		target_region_pos = int(math.floor(random.random() * target_region_size))  # rand_target_pos[num]
		#target_region_pos = int(math.floor(random.normalvariate(0, 1) * target_region_size))
		templet_length = templet_length_array[num]
		#rand_templet_length = random.random()  # 后来加的
		#templet_length = int(rand_templet_length * templetmin + (1 - rand_templet_length) * templetmax)  # rand_templet_length[num]
		# print "初定模板长度：", templet_length
		target_chro_pos, onechr_in_list = capture_to_chromo(target_region_list, target_region_pos, contiguous_point_depth_densityorparameter)
		# print "初定染色体截取位置及染色体：", target_chro_pos, onechr_in_list
		loop = 1
		if customized_dict["depth_pattern_fitting"]:
			loop = create_loop(onechr_in_list, target_chro_pos, templet_length, depth_dict)
		if 'cnv' in customized_dict:
			cnvmultiple = cnv_multiple(onechr_in_list, target_chro_pos, customized_dict)
			#loop = int(round((random.random() * cnvmultiple) * 2)) * loop   # 有可能出现0的问题待解决
			cnvmultiple_floor = int(cnvmultiple)
			cnvmultiple_ceil = cnvmultiple_floor + 1
			cnvmultiple_decimal = cnvmultiple - cnvmultiple_floor
			cnvmultiple_adjust_rdm = random.random()
			if cnvmultiple_adjust_rdm < cnvmultiple_decimal:
				loop = cnvmultiple_ceil * loop
			else:
				loop = cnvmultiple_floor * loop
			# print "loop次数：", loop
		for l in range(loop):
			once_target_chro_pos = target_chro_pos + int(round((random.random() - 0.5) * 0.4 * templet_length))
			once_templet_length = templet_length + int(round((random.random() - 0.5) * 0.1 * templet_length))
			# print "本次模板长度及染色体截取位置：", once_templet_length, once_target_chro_pos
			# onechr_in_list = ">chr1"
			begin, ini_sequence = load_initial_sequence(hg19index_dict, once_target_chro_pos, onechr_in_list, once_templet_length, customized_dict["delete"], customized_dict["tandem_repeat"])
			# print "开始截取位置及取到的序列：", begin, ini_sequence
			begin, ini_sequence = insert_imitate(hg19index_dict, begin, ini_sequence, onechr_in_list, customized_dict["insert"])
			ini_sequence = snv_imitate(begin, ini_sequence, onechr_in_list, customized_dict["snv"])
			# print "snv模拟后序列：",ini_sequence
			ini_sequence = inversion_imitate(begin, ini_sequence, onechr_in_list, customized_dict["inversion"], hg19index_dict)
			ini_sequence = fuse_imitate(onechr_in_list, ini_sequence, begin, customized_dict["fusion"], hg19index_dict)
			# print "基因融合模拟后序列：",ini_sequence
			if random.random() < 0.5:
				ini_sequence = complementary(ini_sequence[::-1])
			reads_making(ini_sequence, customized_dict["config"])
		percent=float(num)*100/float(appro_templet_amount)
		sys.stdout.write("%.4f"%percent)
		sys.stdout.write("%\r")
		sys.stdout.flush()
	globalname.fastq1.close()
	globalname.fastq2.close()


#def main(reference_data, customized_data, output_dir):

def generate_simdata(id, TMBspatial_outpath, TMBspatial_snv=None, TMBspatial_cnv=None, redepth=None):
	#print "os.path.dirname(os.path.realpath(__file__))=%s" % os.path.dirname(os.path.realpath(__file__))



	args, unknown = parser.parse_known_args()
	frequence_type = args.frequence_type
	depth_pattern_fitting = args.depth_pattern_fitting
	fragment_length = args.fragment_length
	fragment_insert_base = args.fragment_insert_base
	if redepth != None:
		depth = redepth
	else:
		depth = args.depth
	pair_end = args.pair_end
	readlen = args.readlen
	#genome = args.genome
	seq_error_rate = args.seq_error_rate
	duplication_rate = args.duplication_rate
	random_index1 = args.random_index1
	random_index2 = args.random_index2
	sqing_cycles = args.sqing_cycles
	read1_adapter = args.read1_adapter
	read2_adapter = args.read2_adapter
	template_len_min = args.template_len_min
	template_len_max = args.template_len_max
	normal_base_qual_min = args.normal_base_qual_min
	normal_base_qual_max = args.normal_base_qual_max
	seq_error_qual_min = args.seq_error_qual_min
	seq_error_qual_max = args.seq_error_qual_max
	fusion = args.fusion
	if TMBspatial_snv != None:
		snv = TMBspatial_snv
	else:
		snv = args.snv
	if TMBspatial_cnv != None:
		cnv = TMBspatial_cnv
	else:
		cnv = args.cnv
	inversion = args.inversion
	delete = args.delete
	tandem_repeat = args.tandem_repeat
	insert = args.insert
	sv_file = args.sv_file
	bed_file = args.bed_file
	multiprocess = args.multiprocess
	outpath = args.outpath
	bamtocnv = args.bamtocnv
	contiguous_point_depth_model = args.contiguous_point_depth_model
	templateLengthDistributionModel_path = args.templateLengthDistributionModel_path
	baseQualityDistributionModel_path = args.baseQualityDistributionModel_path
	overallDepthDistributionModel_path = args.overallDepthDistributionModel_path
	if baseQualityDistributionModel_path != None:
		baseQualityDistribution_arr = np.loadtxt(baseQualityDistributionModel_path + "/base_quality_distribution.txt")
		globalname.reads_baseQualityDistribution_arr = baseQualityDistribution_arr.tolist()
	else:
		globalname.reads_baseQualityDistribution_arr = None
	if contiguous_point_depth_model != None:
		contiguous_point_depth_model = contiguous_point_depth_model.strip()
		contiguous_point_depth_model = contiguous_point_depth_model.rstrip("\/")
		contiguous_point_depth_model = contiguous_point_depth_model + "/"
		contiguous_point_depth_file_path = contiguous_point_depth_model + "contiguous_point_depth.txt"
		contiguous_point_depth_parameter_path = contiguous_point_depth_model + "parameter.txt"
	if bamtocnv != None:
		s_array = bamtocnv.split(",")
		bam = s_array[0]
		if outpath != None:
			bam_path = outpath
			bam_path = bam_path.strip()
			bam_path = bam_path.rstrip("\/")
			isExists=os.path.exists(bam_path)
			if not isExists:
				os.makedirs(bam_path)
		else:
			bam_path = TMBspatial_outpath + '/' + id
			os.makedirs(bam_path)
		bam_path = bam_path + '/'
		ratio = []
		chrom = []
		start = []
		end = []
		ratio_array = s_array[1].split(":")
		chrom_array = s_array[2].split(":")
		start_array = s_array[3].split(":")
		end_array = s_array[4].split(":")
		for ratio_one in ratio_array:
			ratio.append(float(ratio_one))
		for chrom_one in chrom_array:
			chrom.append(chrom_one)
		for start_one in start_array:
			start.append(int(start_one))
		for end_one in end_array:
			end.append(int(end_one))
		graBam(bam, bam_path, ratio, chrom, start, end)
		sys.exit()
	if outpath != None:
		out = outpath
		out = out.strip()
		out = out.rstrip("\/")
		isExists=os.path.exists(out)
		if not isExists:
			os.makedirs(out)
	else:
		#out = time.strftime('%Y%m%d%H%M%S', time.localtime(time.time())) + '/'
		out = TMBspatial_outpath + '/' + id
		os.makedirs(out)
	out = out + '/'
	if sv_file != None:
		inversion, delete, tandem_repeat, fusion = get_mutation(sv_file)
	#print frequence_type, depth_pattern_fitting, fragment_length, depth, pair_end, readlen, assembly, seq_error_rate, duplication_rate, random_index1, random_index2, read1_adapter, read2_adapter, template_len_min, template_len_max, normal_base_qual_min, normal_base_qual_max, seq_error_qual_min, seq_error_qual_max, fusion, snv, cnv, inversion, delete, tandem_repeat, insert, sv_file
   # print inversion, delete, tandem_repeat, fusion



	starttime = time.time()
	chr1 = []
	chr2 = []
	chr3 = []
	chr4 = []
	chr5 = []
	chr6 = []
	chr7 = []
	chr8 = []
	chr9 = []
	chr10 = []
	chr11 = []
	chr12 = []
	chr13 = []
	chr14 = []
	chr15 = []
	chr16 = []
	chr17 = []
	chr18 = []
	chr19 = []
	chr20 = []
	chr21 = []
	chr22 = []
	chrX = []
	chrY = []
	#customized_jsondata = "customizeddata.json"  # "C:\Users\Administrator\PycharmProjects\digitalsimulator\wangmiao_customizeddata.json"
	#customized_dict = load_json_file(customized_jsondata)  # »ñֵäʽÅÖÎ¼þ
	customized_dict = {u'cnv': [], u'snv': [], u'inversion': [], u'insert': [], u'tandem_repeat': [], u'depth_pattern_fitting': False, u'fragment_length': 5, u'fusion': [], u'frequence_type': u'A', u'config': {u'read1_adapter': u'', u'reads_length': 80, u'template_len': {u'max': 600, u'min': 0}, u'read2_adapter': u'', u'pair_end': True, u'normal_base_qual': {u'max': 36, u'min': 33}, u'seq_error_qual': {u'max': 19, u'min': 7}, u'depth': 200, u'duplication_rate': 3.0, u'sqing_cycles': 0, u'reads1_index': True, u'reads2_index': True, u'ref_genome': u'hg19', u'seq_error_rate': 0.001}, u'delete': []}
	if frequence_type != None:
		customized_dict["frequence_type"] = frequence_type
	if depth_pattern_fitting != None:
		if depth_pattern_fitting == 0:
			customized_dict["depth_pattern_fitting"] = False
		else:
			customized_dict["depth_pattern_fitting"] = True
	if fragment_length != None:
		customized_dict["fragment_length"] = fragment_length
	if depth != None:
		customized_dict["config"]["depth"] = depth
	if pair_end != None:
		if pair_end == 0:
			customized_dict["config"]["pair_end"] = False
		else:
			customized_dict["config"]["pair_end"] = True
	if readlen != None:
		customized_dict["config"]["reads_length"] = readlen
	#if genome != None:
		#customized_dict["config"]["ref_genome"] = genome
	if seq_error_rate != None:
		customized_dict["config"]["seq_error_rate"] = seq_error_rate
	if duplication_rate != None:
		customized_dict["config"]["duplication_rate"] = duplication_rate
	if random_index1 != None:
		if random_index1 == 0:
			customized_dict["config"]["reads1_index"] = False
		else:
			customized_dict["config"]["reads1_index"] = True
	if random_index2 != None:
		if random_index2 == 0:
			customized_dict["config"]["reads2_index"] = False
		else:
			customized_dict["config"]["reads2_index"] = True
	if sqing_cycles != None:
		customized_dict["config"]["sqing_cycles"] = sqing_cycles
	if read1_adapter != None:
		customized_dict["config"]["read1_adapter"] = read1_adapter
	if read2_adapter != None:
		customized_dict["config"]["read2_adapter"] = read2_adapter
	if template_len_min != None:
		customized_dict["config"]["template_len"]["min"] = template_len_min
	if template_len_max != None:
		customized_dict["config"]["template_len"]["max"] = template_len_max
	if normal_base_qual_min != None:
		customized_dict["config"]["normal_base_qual"]["min"] = normal_base_qual_min
	if normal_base_qual_max != None:
		customized_dict["config"]["normal_base_qual"]["max"] = normal_base_qual_max
	if seq_error_qual_min != None:
		customized_dict["config"]["seq_error_qual"]["min"] = seq_error_qual_min
	if seq_error_qual_max != None:
		customized_dict["config"]["seq_error_qual"]["max"] = seq_error_qual_max
	if fusion != None:
		fusion_arr = fusion.split(":")
		for fusion_one in fusion_arr:
			fusion_one_arr = fusion_one.split(",")
			#print fusion_one_arr[1]
			#print eval(fusion_one_arr[1])
			if eval(fusion_one_arr[1]):
				flag = True
				for one_in_this_chr in eval(fusion_one_arr[1]):
					if (int(fusion_one_arr[3]) - 50 >= one_in_this_chr[1]) & (int(fusion_one_arr[3]) - 50 <= one_in_this_chr[2]):
						one_in_this_chr[2] = int(fusion_one_arr[3]) + 50
						flag = False
					elif (int(fusion_one_arr[3]) + 50 >= one_in_this_chr[1]) & (int(fusion_one_arr[3]) + 50 <= one_in_this_chr[2]):
						one_in_this_chr[1] = int(fusion_one_arr[3]) - 50
						flag = False
				if flag:
					 eval(fusion_one_arr[1]).append([fusion_one_arr[1], int(fusion_one_arr[3]) - 50, int(fusion_one_arr[3]) + 50])
			else:
				eval(fusion_one_arr[1]).append([fusion_one_arr[1], int(fusion_one_arr[3]) - 50, int(fusion_one_arr[3]) + 50])
			if eval(fusion_one_arr[4]):
				flag = True
				for one_in_this_chr in eval(fusion_one_arr[4]):
					if (int(fusion_one_arr[6]) - 50 >= one_in_this_chr[1]) & (int(fusion_one_arr[6]) - 50 <= one_in_this_chr[2]):
						one_in_this_chr[2] = int(fusion_one_arr[6]) + 50
						flag = False
					elif (int(fusion_one_arr[6]) + 50 >= one_in_this_chr[1]) & (int(fusion_one_arr[6]) + 50 <= one_in_this_chr[2]):
						one_in_this_chr[1] = int(fusion_one_arr[6]) - 50
						flag = False
				if flag:
					 eval(fusion_one_arr[4]).append([fusion_one_arr[4], int(fusion_one_arr[6]) - 50, int(fusion_one_arr[6]) + 50])
			else:
				eval(fusion_one_arr[4]).append([fusion_one_arr[4], int(fusion_one_arr[6]) - 50, int(fusion_one_arr[6]) + 50])
			dic = {'name': fusion_one_arr[0],
				   'left': {'chrom': fusion_one_arr[1],
							'strand': fusion_one_arr[2],
							'pos': int(fusion_one_arr[3])},
				   'right': {'chrom': fusion_one_arr[4],
							 'strand': fusion_one_arr[5],
							 'pos': int(fusion_one_arr[6])},
				   'rate': float(fusion_one_arr[7]),
				   'if_insert': int(fusion_one_arr[8]) != 0}
			customized_dict["fusion"].append(dic)
	if snv != None:
		snv_arr = snv.split(":")
		for snv_one in snv_arr:
			snv_one_arr = snv_one.split(",")
			if eval(snv_one_arr[1]):
				flag = True
				for one_in_this_chr in eval(snv_one_arr[1]):
					if (int(snv_one_arr[2]) - 50 >= one_in_this_chr[1]) & (int(snv_one_arr[2]) - 50 <= one_in_this_chr[2]):
						one_in_this_chr[2] = int(snv_one_arr[2]) + 50
						flag = False
					elif (int(snv_one_arr[2]) + 50 >= one_in_this_chr[1]) & (int(snv_one_arr[2]) + 50 <= one_in_this_chr[2]):
						one_in_this_chr[1] = int(snv_one_arr[2]) - 50
						flag = False
				if flag:
					 eval(snv_one_arr[1]).append([snv_one_arr[1], int(snv_one_arr[2]) - 50, int(snv_one_arr[2]) + 50])
			else:
				eval(snv_one_arr[1]).append([snv_one_arr[1], int(snv_one_arr[2]) - 50, int(snv_one_arr[2]) + 50])
			dic = {'name': snv_one_arr[0],
				   'chrom': snv_one_arr[1],
				   'pos': int(snv_one_arr[2]),
				   'ref': snv_one_arr[3],
				   'alt': snv_one_arr[4],
				   'rate': float(snv_one_arr[5])}
			customized_dict["snv"].append(dic)
	if cnv != None:
		cnv_arr = cnv.split(":")
		for cnv_one in cnv_arr:
			cnv_one_arr = cnv_one.split(",")
			if eval(cnv_one_arr[1]):
				flag = True
				for one_in_this_chr in eval(cnv_one_arr[1]):
					if (int(cnv_one_arr[2]) - 50 >= one_in_this_chr[1]) & (int(cnv_one_arr[2]) - 50 <= one_in_this_chr[2]):
						one_in_this_chr[2] = int(cnv_one_arr[2]) + 50
						flag = False
					elif (int(cnv_one_arr[2]) + 50 >= one_in_this_chr[1]) & (int(cnv_one_arr[2]) + 50 <= one_in_this_chr[2]):
						one_in_this_chr[1] = int(cnv_one_arr[2]) - 50
						flag = False
				if flag:
					 eval(cnv_one_arr[1]).append([cnv_one_arr[1], int(cnv_one_arr[2]) - 50, int(cnv_one_arr[2]) + 50])
			else:
				eval(cnv_one_arr[1]).append([cnv_one_arr[1], int(cnv_one_arr[2]) - 50, int(cnv_one_arr[2]) + 50])
			if eval(cnv_one_arr[1]):
				flag = True
				for one_in_this_chr in eval(cnv_one_arr[1]):
					if (int(cnv_one_arr[3]) - 50 >= one_in_this_chr[1]) & (int(cnv_one_arr[3]) - 50 <= one_in_this_chr[2]):
						one_in_this_chr[2] = int(cnv_one_arr[3]) + 50
						flag = False
					elif (int(cnv_one_arr[3]) + 50 >= one_in_this_chr[1]) & (int(cnv_one_arr[3]) + 50 <= one_in_this_chr[2]):
						one_in_this_chr[1] = int(cnv_one_arr[3]) - 50
						flag = False
				if flag:
					 eval(cnv_one_arr[1]).append([cnv_one_arr[1], int(cnv_one_arr[3]) - 50, int(cnv_one_arr[3]) + 50])
			else:
				eval(cnv_one_arr[1]).append([cnv_one_arr[1], int(cnv_one_arr[3]) - 50, int(cnv_one_arr[3]) + 50])
			dic = {'name': cnv_one_arr[0],
				   'chrom': cnv_one_arr[1],
				   'start': int(cnv_one_arr[2]),
				   'end': int(cnv_one_arr[3]),
				   'copy': float(cnv_one_arr[4])}
			customized_dict["cnv"].append(dic)
	if inversion != None:
		inversion_arr = inversion.split(":")
		for inversion_one in inversion_arr:
			inversion_one_arr = inversion_one.split(",")
			if eval(inversion_one_arr[1]):
				flag = True
				for one_in_this_chr in eval(inversion_one_arr[1]):
					if (int(inversion_one_arr[2]) - 50 >= one_in_this_chr[1]) & (int(inversion_one_arr[2]) - 50 <= one_in_this_chr[2]):
						one_in_this_chr[2] = int(inversion_one_arr[2]) + 50
						flag = False
					elif (int(inversion_one_arr[2]) + 50 >= one_in_this_chr[1]) & (int(inversion_one_arr[2]) + 50 <= one_in_this_chr[2]):
						one_in_this_chr[1] = int(inversion_one_arr[2]) - 50
						flag = False
				if flag:
					 eval(inversion_one_arr[1]).append([inversion_one_arr[1], int(inversion_one_arr[2]) - 50, int(inversion_one_arr[2]) + 50])
			else:
				eval(inversion_one_arr[1]).append([inversion_one_arr[1], int(inversion_one_arr[2]) - 50, int(inversion_one_arr[2]) + 50])
			if eval(inversion_one_arr[1]):
				flag = True
				for one_in_this_chr in eval(inversion_one_arr[1]):
					if (int(inversion_one_arr[3]) - 50 >= one_in_this_chr[1]) & (int(inversion_one_arr[3]) - 50 <= one_in_this_chr[2]):
						one_in_this_chr[2] = int(inversion_one_arr[3]) + 50
						flag = False
					elif (int(inversion_one_arr[3]) + 50 >= one_in_this_chr[1]) & (int(inversion_one_arr[3]) + 50 <= one_in_this_chr[2]):
						one_in_this_chr[1] = int(inversion_one_arr[3]) - 50
						flag = False
				if flag:
					 eval(inversion_one_arr[1]).append([inversion_one_arr[1], int(inversion_one_arr[3]) - 50, int(inversion_one_arr[3]) + 50])
			else:
				eval(inversion_one_arr[1]).append([inversion_one_arr[1], int(inversion_one_arr[3]) - 50, int(inversion_one_arr[3]) + 50])
			dic = {'name': inversion_one_arr[0],
				   'chrom': inversion_one_arr[1],
				   'start': int(inversion_one_arr[2]),
				   'end': int(inversion_one_arr[3]),
				   'rate': float(inversion_one_arr[4])}
			customized_dict["inversion"].append(dic)
	if delete != None:
		delete_arr = delete.split(":")
		for delete_one in delete_arr:
			delete_one_arr = delete_one.split(",")
			if eval(delete_one_arr[1]):
				flag = True
				for one_in_this_chr in eval(delete_one_arr[1]):
					if (int(delete_one_arr[2]) - 50 >= one_in_this_chr[1]) & (int(delete_one_arr[2]) - 50 <= one_in_this_chr[2]):
						one_in_this_chr[2] = int(delete_one_arr[2]) + 50
						flag = False
					elif (int(delete_one_arr[2]) + 50 >= one_in_this_chr[1]) & (int(delete_one_arr[2]) + 50 <= one_in_this_chr[2]):
						one_in_this_chr[1] = int(delete_one_arr[2]) - 50
						flag = False
				if flag:
					 eval(delete_one_arr[1]).append([delete_one_arr[1], int(delete_one_arr[2]) - 50, int(delete_one_arr[2]) + 50])
			else:
				eval(delete_one_arr[1]).append([delete_one_arr[1], int(delete_one_arr[2]) - 50, int(delete_one_arr[2]) + 50])
			if eval(delete_one_arr[1]):
				flag = True
				for one_in_this_chr in eval(delete_one_arr[1]):
					if (int(delete_one_arr[3]) - 50 >= one_in_this_chr[1]) & (int(delete_one_arr[3]) - 50 <= one_in_this_chr[2]):
						one_in_this_chr[2] = int(delete_one_arr[3]) + 50
						flag = False
					elif (int(delete_one_arr[3]) + 50 >= one_in_this_chr[1]) & (int(delete_one_arr[3]) + 50 <= one_in_this_chr[2]):
						one_in_this_chr[1] = int(delete_one_arr[3]) - 50
						flag = False
				if flag:
					 eval(delete_one_arr[1]).append([delete_one_arr[1], int(delete_one_arr[3]) - 50, int(delete_one_arr[3]) + 50])
			else:
				eval(delete_one_arr[1]).append([delete_one_arr[1], int(delete_one_arr[3]) - 50, int(delete_one_arr[3]) + 50])
			dic = {'name': delete_one_arr[0],
				   'chrom': delete_one_arr[1],
				   'start': int(delete_one_arr[2]),
				   'end': int(delete_one_arr[3]),
				   'rate': float(delete_one_arr[4]),
				   'if_insert': int(delete_one_arr[5]) != 0}
			customized_dict["delete"].append(dic)
	if tandem_repeat != None:
		tandem_repeat_arr = tandem_repeat.split(":")
		for tandem_repeat_one in tandem_repeat_arr:
			tandem_repeat_one_arr = tandem_repeat_one.split(",")
			if eval(tandem_repeat_one_arr[1]):
				flag = True
				for one_in_this_chr in eval(tandem_repeat_one_arr[1]):
					if (int(tandem_repeat_one_arr[2]) - 50 >= one_in_this_chr[1]) & (int(tandem_repeat_one_arr[2]) - 50 <= one_in_this_chr[2]):
						one_in_this_chr[2] = int(tandem_repeat_one_arr[2]) + 50
						flag = False
					elif (int(tandem_repeat_one_arr[2]) + 50 >= one_in_this_chr[1]) & (int(tandem_repeat_one_arr[2]) + 50 <= one_in_this_chr[2]):
						one_in_this_chr[1] = int(tandem_repeat_one_arr[2]) - 50
						flag = False
				if flag:
					 eval(tandem_repeat_one_arr[1]).append([tandem_repeat_one_arr[1], int(tandem_repeat_one_arr[2]) - 50, int(tandem_repeat_one_arr[2]) + 50])
			else:
				eval(tandem_repeat_one_arr[1]).append([tandem_repeat_one_arr[1], int(tandem_repeat_one_arr[2]) - 50, int(tandem_repeat_one_arr[2]) + 50])
			if eval(tandem_repeat_one_arr[1]):
				flag = True
				for one_in_this_chr in eval(tandem_repeat_one_arr[1]):
					if (int(tandem_repeat_one_arr[3]) - 50 >= one_in_this_chr[1]) & (int(tandem_repeat_one_arr[3]) - 50 <= one_in_this_chr[2]):
						one_in_this_chr[2] = int(tandem_repeat_one_arr[3]) + 50
						flag = False
					elif (int(tandem_repeat_one_arr[3]) + 50 >= one_in_this_chr[1]) & (int(tandem_repeat_one_arr[3]) + 50 <= one_in_this_chr[2]):
						one_in_this_chr[1] = int(tandem_repeat_one_arr[3]) - 50
						flag = False
				if flag:
					 eval(tandem_repeat_one_arr[1]).append([tandem_repeat_one_arr[1], int(tandem_repeat_one_arr[3]) - 50, int(tandem_repeat_one_arr[3]) + 50])
			else:
				eval(tandem_repeat_one_arr[1]).append([tandem_repeat_one_arr[1], int(tandem_repeat_one_arr[3]) - 50, int(tandem_repeat_one_arr[3]) + 50])
			dic = {'name': tandem_repeat_one_arr[0],
				   'chrom': tandem_repeat_one_arr[1],
				   'start': int(tandem_repeat_one_arr[2]),
				   'end': int(tandem_repeat_one_arr[3]),
				   'copy': int(tandem_repeat_one_arr[4]),
				   'rate': float(tandem_repeat_one_arr[5])}
			customized_dict["tandem_repeat"].append(dic)
	if insert != None:
		insert_arr = insert.split(":")
		for insert_one in insert_arr:
			insert_one_arr = insert_one.split(",")
			if eval(insert_one_arr[1]):
				flag = True
				for one_in_this_chr in eval(insert_one_arr[1]):
					if (int(insert_one_arr[2]) - 50 >= one_in_this_chr[1]) & (int(insert_one_arr[2]) - 50 <= one_in_this_chr[2]):
						one_in_this_chr[2] = int(insert_one_arr[2]) + 50
						flag = False
					elif (int(insert_one_arr[2]) + 50 >= one_in_this_chr[1]) & (int(insert_one_arr[2]) + 50 <= one_in_this_chr[2]):
						one_in_this_chr[1] = int(insert_one_arr[2]) - 50
						flag = False
				if flag:
					 eval(insert_one_arr[1]).append([insert_one_arr[1], int(insert_one_arr[2]) - 50, int(insert_one_arr[2]) + 50])
			else:
				eval(insert_one_arr[1]).append([insert_one_arr[1], int(insert_one_arr[2]) - 50, int(insert_one_arr[2]) + 50])
			dic = {'name': insert_one_arr[0],
				   'chrom': insert_one_arr[1],
				   'pos': int(insert_one_arr[2]),
				   'rate': float(insert_one_arr[3])}
			customized_dict["insert"].append(dic)
	if (fusion == None) & (snv == None) & (cnv == None) & (inversion == None) & (delete == None) & (tandem_repeat == None ) & (insert == None):
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
		chr_total = 3
		for one_in_chr_total in range(chr_total):
			chr_rnd = choice(chr_list)
			start_rnd = random.randint(chr_dict[chr_rnd][0], chr_dict[chr_rnd][1])
			end_rnd = min(start_rnd + 100, chr_dict[chr_rnd][1])
			eval(chr_rnd).append([chr_rnd, start_rnd, end_rnd])



	write_bed = open(out + "tmp_targetregioncapture.bed", "a")
	if chr1:
		for one_in_this_chr in chr1:
			to_write_bed = one_in_this_chr[0] + '\t' + str(one_in_this_chr[1]) + '\t' + str(one_in_this_chr[2]) + '\n'
			write_bed.write(to_write_bed)
	if chr2:
		for one_in_this_chr in chr2:
			to_write_bed = one_in_this_chr[0] + '\t' + str(one_in_this_chr[1]) + '\t' + str(one_in_this_chr[2]) + '\n'
			write_bed.write(to_write_bed)
	if chr3:
		for one_in_this_chr in chr3:
			to_write_bed = one_in_this_chr[0] + '\t' + str(one_in_this_chr[1]) + '\t' + str(one_in_this_chr[2]) + '\n'
			write_bed.write(to_write_bed)
	if chr4:
		for one_in_this_chr in chr4:
			to_write_bed = one_in_this_chr[0] + '\t' + str(one_in_this_chr[1]) + '\t' + str(one_in_this_chr[2]) + '\n'
			write_bed.write(to_write_bed)
	if chr5:
		for one_in_this_chr in chr5:
			to_write_bed = one_in_this_chr[0] + '\t' + str(one_in_this_chr[1]) + '\t' + str(one_in_this_chr[2]) + '\n'
			write_bed.write(to_write_bed)
	if chr6:
		for one_in_this_chr in chr6:
			to_write_bed = one_in_this_chr[0] + '\t' + str(one_in_this_chr[1]) + '\t' + str(one_in_this_chr[2]) + '\n'
			write_bed.write(to_write_bed)
	if chr7:
		for one_in_this_chr in chr7:
			to_write_bed = one_in_this_chr[0] + '\t' + str(one_in_this_chr[1]) + '\t' + str(one_in_this_chr[2]) + '\n'
			write_bed.write(to_write_bed)
	if chr8:
		for one_in_this_chr in chr8:
			to_write_bed = one_in_this_chr[0] + '\t' + str(one_in_this_chr[1]) + '\t' + str(one_in_this_chr[2]) + '\n'
			write_bed.write(to_write_bed)
	if chr9:
		for one_in_this_chr in chr9:
			to_write_bed = one_in_this_chr[0] + '\t' + str(one_in_this_chr[1]) + '\t' + str(one_in_this_chr[2]) + '\n'
			write_bed.write(to_write_bed)
	if chr10:
		for one_in_this_chr in chr10:
			to_write_bed = one_in_this_chr[0] + '\t' + str(one_in_this_chr[1]) + '\t' + str(one_in_this_chr[2]) + '\n'
			write_bed.write(to_write_bed)
	if chr11:
		for one_in_this_chr in chr11:
			to_write_bed = one_in_this_chr[0] + '\t' + str(one_in_this_chr[1]) + '\t' + str(one_in_this_chr[2]) + '\n'
			write_bed.write(to_write_bed)
	if chr12:
		for one_in_this_chr in chr12:
			to_write_bed = one_in_this_chr[0] + '\t' + str(one_in_this_chr[1]) + '\t' + str(one_in_this_chr[2]) + '\n'
			write_bed.write(to_write_bed)
	if chr13:
		for one_in_this_chr in chr13:
			to_write_bed = one_in_this_chr[0] + '\t' + str(one_in_this_chr[1]) + '\t' + str(one_in_this_chr[2]) + '\n'
			write_bed.write(to_write_bed)
	if chr14:
		for one_in_this_chr in chr14:
			to_write_bed = one_in_this_chr[0] + '\t' + str(one_in_this_chr[1]) + '\t' + str(one_in_this_chr[2]) + '\n'
			write_bed.write(to_write_bed)
	if chr15:
		for one_in_this_chr in chr15:
			to_write_bed = one_in_this_chr[0] + '\t' + str(one_in_this_chr[1]) + '\t' + str(one_in_this_chr[2]) + '\n'
			write_bed.write(to_write_bed)
	if chr16:
		for one_in_this_chr in chr16:
			to_write_bed = one_in_this_chr[0] + '\t' + str(one_in_this_chr[1]) + '\t' + str(one_in_this_chr[2]) + '\n'
			write_bed.write(to_write_bed)
	if chr17:
		for one_in_this_chr in chr17:
			to_write_bed = one_in_this_chr[0] + '\t' + str(one_in_this_chr[1]) + '\t' + str(one_in_this_chr[2]) + '\n'
			write_bed.write(to_write_bed)
	if chr18:
		for one_in_this_chr in chr18:
			to_write_bed = one_in_this_chr[0] + '\t' + str(one_in_this_chr[1]) + '\t' + str(one_in_this_chr[2]) + '\n'
			write_bed.write(to_write_bed)
	if chr19:
		for one_in_this_chr in chr19:
			to_write_bed = one_in_this_chr[0] + '\t' + str(one_in_this_chr[1]) + '\t' + str(one_in_this_chr[2]) + '\n'
			write_bed.write(to_write_bed)
	if chr20:
		for one_in_this_chr in chr20:
			to_write_bed = one_in_this_chr[0] + '\t' + str(one_in_this_chr[1]) + '\t' + str(one_in_this_chr[2]) + '\n'
			write_bed.write(to_write_bed)
	if chr21:
		for one_in_this_chr in chr21:
			to_write_bed = one_in_this_chr[0] + '\t' + str(one_in_this_chr[1]) + '\t' + str(one_in_this_chr[2]) + '\n'
			write_bed.write(to_write_bed)
	if chr22:
		for one_in_this_chr in chr22:
			to_write_bed = one_in_this_chr[0] + '\t' + str(one_in_this_chr[1]) + '\t' + str(one_in_this_chr[2]) + '\n'
			write_bed.write(to_write_bed)
	if chrX:
		for one_in_this_chr in chrX:
			to_write_bed = one_in_this_chr[0] + '\t' + str(one_in_this_chr[1]) + '\t' + str(one_in_this_chr[2]) + '\n'
			write_bed.write(to_write_bed)
	if chrY:
		for one_in_this_chr in chrY:
			to_write_bed = one_in_this_chr[0] + '\t' + str(one_in_this_chr[1]) + '\t' + str(one_in_this_chr[2]) + '\n'
			write_bed.write(to_write_bed)
	write_bed.close()
  
  

	#for zidian in customized_dict:
		#print zidian
	#for zidian in customized_dict["config"]:
		#print zidian
		# if zidian["conference_date"] == "":
		#     zidian["conference_date"] = "0001-01-01T00:00:00Z"
		# if zidian["publication_date"] == "":
		#     zidian["publication_date"] = "0001-01-01T00:00:00Z"
	# after = before
	# with open("wangmiao_customizeddata.json", 'wb') as f:
	#     data = json.dump(after, f)
	#print customized_dict



	hg19index_jsondata = "hg19indexdata.json"  # "C:\Users\Administrator\PycharmProjects\digitalsimulator\hg19indexdata.json"
	hg19index_dict = load_json_file(hg19index_jsondata)
	if overallDepthDistributionModel_path != None:
		depth_dict = load_json_file(overallDepthDistributionModel_path + "/depth_distribution.json")
	else:
		depth_dict = load_json_file("overall_depth.json")
	hgbed = 'hg19.bed'  # 'C:\Users\Administrator\PycharmProjects\digitalsimulator\\hg19.bed'
	if bed_file != None:
		targetregioncapturebed = bed_file
	else:
		targetregioncapturebed = out + 'tmp_targetregioncapture.bed'  # 'C:\Users\Administrator\PycharmProjects\digitalsimulator\\wangmiao_targetregioncapture.bed'
	target_region_size, target_region_list = loadtargetregioncapture(hgbed, targetregioncapturebed)
	templetmax = customized_dict["config"]["template_len"]["max"]
	templetmin = customized_dict["config"]["template_len"]["min"]
	templetmid = (templetmax + templetmin) / 2
	appro_templet_amount = int(math.floor((customized_dict["config"]["depth"] * 10 * target_region_size) / (
	templetmid * customized_dict["config"]["duplication_rate"] * 6)))    #为减小数据量做了修改，加了个“*6”
	print(appro_templet_amount)
	# print appro_templet_amount
	# appro_templet_amount = 100            #½âÄ´æÌºóýfragment_length = customized_dict["fragment_length"]
	if fragment_insert_base != None:
		insert_fragment = fragment_insert_base
	else:
		fragment_length = customized_dict["fragment_length"]
		insert_fragment = insert_making(fragment_length)
	globalname.tgp_dict = {'chr1': {}, 'chr2': {}, 'chr3': {}, 'chr4': {}, 'chr5': {}, 'chr6': {}, 'chr7': {}, 'chr8': {}, 'chr9': {}, 'chr10': {}, 'chr11': {}, 'chr12': {}, 'chr13': {}, 'chr14': {}, 'chr15': {}, 'chr16': {}, 'chr17': {}, 'chr18': {}, 'chr19': {}, 'chr20': {}, 'chr21': {}, 'chr22': {}, 'chrX': {}, 'chrY': {}}
	rdm = int(time.time())
	# print rdm
	print(insert_fragment)



	contiguous_point_depth_densityorparameter = None
	if contiguous_point_depth_model != None:
		contiguous_point_depth_file = open(contiguous_point_depth_file_path)
		lines = contiguous_point_depth_file.readlines()
		contiguous_point_depth_density = []
		contiguous_point_depths = []
		for line in lines:
			contiguous_point_depths.append(int(line.strip('\n')))
		i = 0
		contiguous_point_depths_len = float(len(contiguous_point_depths))
		for contiguous_point_depth in contiguous_point_depths:
			i = i + 1
			for k in range(contiguous_point_depth):
				contiguous_point_depth_density.append(i / contiguous_point_depths_len)
		contiguous_point_depth_file.close()

		contiguous_point_depth_parameter = open(contiguous_point_depth_parameter_path)
		lines = contiguous_point_depth_parameter.readlines()
		contiguous_point_depth_parameter_arr = []
		for line in lines:
			line = line.rstrip("\n")
			contiguous_point_depth_parameter_arr.append(line)
		contiguous_point_depth_parameter.close()

		if contiguous_point_depth_parameter_arr[0] == "multiple":
			contiguous_point_depth_densityorparameter = contiguous_point_depth_density
		else:
			contiguous_point_depth_densityorparameter = contiguous_point_depth_parameter_arr


	record = []
	if multiprocess != None:
		process_num = multiprocess
	else:
		process_num = 10
	for i in range(process_num):
		process = multiprocessing.Process(target=multi_process_write, args=(appro_templet_amount / process_num, i+1, insert_fragment, rdm, customized_dict, hg19index_dict, target_region_size, target_region_list, templetmin, templetmax, contiguous_point_depth_densityorparameter, templateLengthDistributionModel_path,out,depth_dict))
		process.start()
		record.append(process)

	for process in record:
		process.join()



	midtime = time.time()
	#print midtime - starttime
	sys.stdout.write("99.9999%\r")
	sys.stdout.flush()

	os.remove(out + 'tmp_targetregioncapture.bed')

	fq_emerged1 = open(out + 'make1.fq', 'a')
	for i in range(process_num):
		fastq_to_be_merged = out + "makeit" + str(2 * (i + 1) - 1) + ".fq"
		fq_i = open(fastq_to_be_merged)
		while True:
			s = fq_i.read()  # 16 * 1024
			if not s:
				break
			fq_emerged1.write(s)
		fq_i.close()
		os.remove(fastq_to_be_merged)
	fq_emerged1.close()

	if customized_dict["config"]["pair_end"]:
		fq_emerged2 = open(out + 'make2.fq', 'a')
		for i in range(process_num):
			fastq_to_be_merged = out + "makeit" + str(2 * (i + 1)) + ".fq"
			fq_i = open(fastq_to_be_merged)
			while True:
				s = fq_i.read()  # 16 * 1024
				if not s:
					break
				fq_emerged2.write(s)
			fq_i.close()
			os.remove(fastq_to_be_merged)
		fq_emerged2.close()
	else:
		for i in range(process_num):
			fastq_to_be_merged = out + "makeit" + str(2 * (i + 1)) + ".fq"
			os.remove(fastq_to_be_merged)

	#globalname.fastq1.close()
	#globalname.fastq2.close()


	write_shell = open(out + "tmp_bwa.sh", "a")
	samout = out
	# if customized_dict["config"]["pair_end"]:
	# 	command = "bwa mem -t 50  /data/home/std_17/mnt/NL200/refgenomes/tmp_hg19/hg19.fa " + samout + "make1.fq " + samout + "make2.fq > " + out + "m.sam" + '\n' + "samtools view -bS " + samout + "m.sam > " + samout + "mm.bam" + '\n' + "samtools sort " + samout + "mm.bam > " + samout + "mmm.bam" + '\n' + "samtools index " + samout + "mmm.bam" + '\n' + "rm " + '\"' + "$0" + '\"'
	# else:
	#command = "bwa mem -t 50  /data/home/std_17/mnt/NL200/refgenomes/tmp_hg19/hg19.fa " + samout + "make1.fq > " + samout + "m.sam" + '\n' + "samtools view -bS " + samout + "m.sam > " + samout + "mm.bam" + '\n' + "samtools sort " + samout + "mm.bam > " + samout + "mmm.bam" + '\n' + "samtools index " + samout + "mmm.bam" + '\n' + "rm " + '\"' + "$0" + '\"'
	command = "bwa mem -t 50 /data/home/std_17/mnt/NL200/refgenomes/tmp_hg19/hg19.fa " + samout + "make1.fq | " + \
          "samtools view -bS - > " + samout + "m.bam" + '\n' + \
          "samtools sort -@ 10 -o " + samout + "m_sorted.bam " + samout + "m.bam" + '\n' + \
          "samtools markdup -@ 10 -s " + samout + "m_sorted.bam " + samout + "m_dedup.bam" + '\n' + \
          "samtools index " + samout + "m_dedup.bam" + '\n' + \
          "rm " + samout + "m.bam " + samout + "m_sorted.bam" + '\n' + \
          "rm \"" + "$0" + "\""
	#command = "bwa mem -t 50  /mnt/NL200/refgenomes/tmp_hg19/hg19.fa " + out + "make1.fq " + out + "make2.fq > " + out + "m.sam " + "-r "@RG\tID:$sample\tSM:$sample\tLB:WES\tPL:Illumina"" + '\n' + "samtools view -bS " + out + "m.sam > " + out + "mm.bam" + '\n' + "samtools sort " + out + "mm.bam > " + out + "mmm.bam" + '\n' + "samtools index " + out + "mmm.bam" + '\n' + "rm " + '\"' + "$0" + '\"'
	write_shell.write(command)
	write_shell.close()
	os.system("sh " + out + "tmp_bwa.sh")


	endtime = time.time()
	#print endtime - midtime
	#print endtime - starttime
	sys.stdout.write("100%!finish!\n")


def bulksimulator(outfile=None):
	if outfile!= None:
		out = outfile
	else:
		out = time.strftime('%Y%m%d%H%M%S', time.localtime(time.time()))
	generate_simdata('',out)


if __name__ == '__main__':
	bulksimulator()