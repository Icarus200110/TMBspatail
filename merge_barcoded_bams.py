import os
import pysam
import argparse

def merge_barcoded_bam_files(input_dir, output_bam_path):
    """
    将已经添加了条形码的多个BAM文件合并为一个单一的BAM文件
    
    Args:
        input_dir: 包含已添加条形码的BAM文件的目录
        output_bam_path: 输出的合并BAM文件路径
    """
    # 获取目录中的所有BAM文件
    bam_files = []
    for filename in os.listdir(input_dir):
        if filename.endswith('.bam') and not filename.endswith('.bai'):
            bam_files.append(os.path.join(input_dir, filename))
    
    if not bam_files:
        print(f"错误：在目录 {input_dir} 中找不到任何BAM文件")
        return
    
    print(f"找到 {len(bam_files)} 个BAM文件")
    
    # 打开第一个BAM文件以获取头信息
    first_bam = pysam.AlignmentFile(bam_files[0], "rb")
    header = first_bam.header.to_dict()
    
    # 创建输出BAM文件
    output_bam = pysam.AlignmentFile(output_bam_path, "wb", header=header)
    
    # 处理每个BAM文件
    processed_count = 0
    total_reads = 0
    
    for bam_file in bam_files:
        # 打开BAM文件
        input_bam = pysam.AlignmentFile(bam_file, "rb")
        
        # 遍历BAM文件中的每条记录
        file_reads = 0
        for read in input_bam:
            # 将记录写入输出文件
            output_bam.write(read)
            file_reads += 1
        
        input_bam.close()
        processed_count += 1
        total_reads += file_reads
        print(f"已处理文件 {os.path.basename(bam_file)}，读取了 {file_reads} 条记录")
    
    # 关闭输出文件
    output_bam.close()
    
    # 对BAM文件进行排序
    sorted_bam_path = output_bam_path + ".sorted.bam"
    pysam.sort("-o", sorted_bam_path, output_bam_path)
    
    # 用排序后的文件替换原文件
    os.replace(sorted_bam_path, output_bam_path)
    print(f"已对 {output_bam_path} 进行排序")
    
    # 创建索引
    pysam.index(output_bam_path)
    print(f"已为 {output_bam_path} 创建索引")
    
    print(f"处理完成，共合并了 {processed_count} 个BAM文件，总共 {total_reads} 条读取记录")

def main():
    parser = argparse.ArgumentParser(description="将已添加条形码的多个BAM文件合并为一个")
    parser.add_argument("--input_dir", type=str, required=True, help="包含已添加条形码的BAM文件的目录")
    parser.add_argument("--output_bam", type=str, required=True, help="输出的合并BAM文件路径")
    args = parser.parse_args()
    
    merge_barcoded_bam_files(args.input_dir, args.output_bam)

if __name__ == "__main__":
    main()