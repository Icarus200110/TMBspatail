#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
解压VCF目录中的所有.gz文件
"""

import gzip
import os
import glob

vcf_dir = 'd:\\Project\\python\\GSDcreator\\VCF'

for gz_file in glob.glob(os.path.join(vcf_dir, '*.vcf.gz')):
    output_file = gz_file.replace('.gz', '')
    try:
        with gzip.open(gz_file, 'rb') as f_in:
            with open(output_file, 'wb') as f_out:
                f_out.write(f_in.read())
        print(f'Successfully decompressed: {os.path.basename(gz_file)}')
    except Exception as e:
        print(f'Error decompressing {os.path.basename(gz_file)}: {str(e)}')