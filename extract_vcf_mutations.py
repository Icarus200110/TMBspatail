#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
提取VCF文件中的变异信息

该脚本用于处理VCF目录下的所有VCF文件，提取每个样本的变异信息，
包括染色体、位置、参考碱基、变异碱基和等位基因频率。
将每个样本的所有变异信息合并为一行，按照指定格式输出到结果文件中。

输出格式：DOXX:id,chr,POS,REF,ALT,AF,id,chr,POS,REF,ALT,AF...
其中，DO在文件名中提取，id从0递增
"""

import os
import re
import glob

# VCF文件目录
vcf_dir = 'd:/Project/python/GSDcreator/VCF'
# 输出文件路径
output_file = 'd:/Project/python/GSDcreator/mutation_list.txt'

# 正则表达式，用于从文件名中提取DO编号
do_pattern = re.compile(r'DO(\d+)')

# 存储所有样本的变异信息
all_mutations = {}

# 获取所有VCF文件（不包括.gz和.tbi文件）
vcf_files = [f for f in glob.glob(os.path.join(vcf_dir, '*.vcf')) if not f.endswith('.gz') and not f.endswith('.tbi')]

# 处理每个VCF文件
for vcf_file in vcf_files:
    # 从文件名中提取DO编号
    filename = os.path.basename(vcf_file)
    do_match = do_pattern.search(filename)
    
    if do_match:
        do_number = do_match.group(1)  # 提取DO编号
        do_key = f'DO{do_number}'  # 创建DO键
        
        # 初始化该样本的变异列表
        if do_key not in all_mutations:
            all_mutations[do_key] = []
        
        # 读取VCF文件并提取变异信息
        with open(vcf_file, 'r') as f:
            for line in f:
                # 跳过注释行
                if line.startswith('#'):
                    continue
                
                # 解析变异行
                fields = line.strip().split('\t')
                if len(fields) >= 10:  # 确保有足够的字段
                    chrom = fields[0]  # 染色体
                    pos = fields[1]    # 位置
                    ref = fields[3]    # 参考碱基
                    alt = fields[4]    # 变异碱基
                    
                    # 提取等位基因频率 (AF)
                    # 格式字段通常在第9个位置，样本数据在第10个位置
                    format_fields = fields[8].split(':')
                    sample_fields = fields[9].split(':')
                    
                    # 查找AF在格式字段中的索引
                    af_index = -1
                    for i, field in enumerate(format_fields):
                        if field == 'AF':
                            af_index = i
                            break
                    
                    # 如果找到AF字段，提取其值
                    af = ''
                    if af_index != -1 and af_index < len(sample_fields):
                        af = sample_fields[af_index]
                    
                    # 将变异信息添加到样本的变异列表中
                    all_mutations[do_key].append((chrom, pos, ref, alt, af))

# 将结果写入输出文件
with open(output_file, 'w') as out_f:
    # 遍历所有样本
    for do_key in sorted(all_mutations.keys()):
        mutations = all_mutations[do_key]
        
        # 如果样本有变异，则写入文件
        if mutations:
            # 为每个变异分配一个ID（从0开始递增）
            mutation_entries = [do_key]
            for i, (chrom, pos, ref, alt, af) in enumerate(mutations):
                mutation_entries.append(f"{i},{chrom},{pos},{ref},{alt},{af}")
            
            # 将所有变异信息合并为一行
            out_f.write(','.join(mutation_entries) + '\n')

print(f"处理完成！共处理了 {len(vcf_files)} 个VCF文件。")
print(f"结果已保存到 {output_file}")