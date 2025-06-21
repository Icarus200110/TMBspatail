#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
计算每个细胞中SNP的频率，并生成频率矩阵
"""

import os
import numpy as np

# 输入和输出文件路径
input_ad_file = 'd:/Project/python/GSDcreator/snv_test/cellSNP.tag.AD.mtx'
input_dp_file = 'd:/Project/python/GSDcreator/snv_test/cellSNP.tag.DP.mtx'
output_freq_file = 'd:/Project/python/GSDcreator/snv_test/cellSNP.tag.AF.mtx'

# 读取AD矩阵文件的头部信息和数据
ad_header = []
ad_data = []
with open(input_ad_file, 'r') as f:
    for line in f:
        if line.startswith('%'):
            ad_header.append(line)
        else:
            ad_data.append(line)

# 读取DP矩阵文件的头部信息和数据
dp_header = []
dp_data = []
with open(input_dp_file, 'r') as f:
    for line in f:
        if line.startswith('%'):
            dp_header.append(line)
        else:
            dp_data.append(line)

# 解析AD数据，提取行索引、列索引和值
ad_entries = {}
for i, line in enumerate(ad_data):
    if i == 0:  # 矩阵维度信息行
        continue
    parts = line.strip().split()
    if len(parts) == 3:
        row, col, value = int(parts[0]), int(parts[1]), int(parts[2])
        key = (row, col)
        ad_entries[key] = value

# 解析DP数据，提取行索引、列索引和值
dp_entries = {}
for i, line in enumerate(dp_data):
    if i == 0:  # 矩阵维度信息行
        continue
    parts = line.strip().split()
    if len(parts) == 3:
        row, col, value = int(parts[0]), int(parts[1]), int(parts[2])
        key = (row, col)
        dp_entries[key] = value

# 计算频率
freq_entries = []
for key in ad_entries:
    if key in dp_entries and dp_entries[key] > 0:
        row, col = key
        ad_value = ad_entries[key]
        dp_value = dp_entries[key]
        freq = ad_value / dp_value
        freq_entries.append((row, col, freq))

# 获取矩阵维度信息
matrix_dims = ad_data[0].strip().split()
matrix_dims[2] = str(len(freq_entries))  # 更新非零元素数量

# 写入频率矩阵
with open(output_freq_file, 'w') as f:
    # 写入头部信息
    for line in ad_header:
        f.write(line)
    # 写入矩阵维度信息
    f.write('\t'.join(matrix_dims) + '\n')
    # 写入频率数据
    for row, col, freq in sorted(freq_entries, key=lambda x: (x[0], x[1])):
        f.write(f"{row}\t{col}\t{freq:.6f}\n")

print(f'处理完成！')
print(f'频率矩阵已保存为: {output_freq_file}')

# 读取VCF文件获取SNP信息
vcf_file = 'd:/Project/python/GSDcreator/snv_test/cellSNP.base.vcf'
snp_info = []
with open(vcf_file, 'r') as f:
    for line in f:
        if line.startswith('#CHROM'):
            continue
        elif not line.startswith('#'):
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = fields[1]
            ref = fields[3]
            alt = fields[4]
            snp_info.append((chrom, pos, ref, alt))

# 打印SNP信息
print("\nSNP信息:")
for i, (chrom, pos, ref, alt) in enumerate(snp_info, 1):
    print(f"SNP {i}: 染色体 {chrom}, 位置 {pos}, 参考碱基 {ref}, 变异碱基 {alt}")

# 计算每个SNP在所有细胞中的平均频率
snp_avg_freq = {}
for row in range(1, len(snp_info) + 1):
    freqs = []
    for entry in freq_entries:
        if entry[0] == row:  # 如果行索引匹配当前SNP
            freqs.append(entry[2])  # 添加频率值
    if freqs:
        snp_avg_freq[row] = sum(freqs) / len(freqs)
    else:
        snp_avg_freq[row] = 0

print("\n每个SNP在所有细胞中的平均频率:")
for i, (snp_idx, avg_freq) in enumerate(snp_avg_freq.items(), 1):
    if i <= len(snp_info):
        chrom, pos, ref, alt = snp_info[i-1]
        print(f"SNP {i} (染色体 {chrom}, 位置 {pos}): {avg_freq:.6f}")