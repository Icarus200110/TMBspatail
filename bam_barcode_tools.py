import os
import csv
import pysam
import argparse
import sys

def read_cell_coordinates(cell_csv_path):
    """
    读取细胞坐标信息
    
    Args:
        cell_csv_path: 细胞坐标CSV文件路径
    
    Returns:
        字典，键为细胞ID，值为坐标元组(x, y)
    """
    cell_coords = {}
    with open(cell_csv_path, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) >= 2:
                cell_id = row[0]
                # 解析坐标字符串，格式为"(x, y)"
                coord_str = row[1].strip('()')
                x, y = map(int, coord_str.split(','))
                cell_coords[cell_id] = (x, y)
    return cell_coords

def read_barcodes(barcode_csv_path):
    """
    读取10X Genomics格式的条形码信息
    
    Args:
        barcode_csv_path: 条形码CSV文件路径
    
    Returns:
        字典，键为坐标元组(x, y)，值为条形码
    """
    barcodes = {}
    with open(barcode_csv_path, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) >= 3:
                barcode = row[0]
                x = int(row[1])
                y = int(row[2])
                barcodes[(x, y)] = barcode+','+row[1]+','+row[2]
    return barcodes

def add_barcode_to_bam(input_bam_path, output_bam_path, barcode):
    """
    将条形码信息添加到BAM文件中
    
    Args:
        input_bam_path: 输入BAM文件路径
        output_bam_path: 输出BAM文件路径
        barcode: 要添加的条形码
    """
    # 打开输入BAM文件
    input_bam = pysam.AlignmentFile(input_bam_path, "rb")
    
    # 创建输出BAM文件，使用与输入相同的头信息
    output_bam = pysam.AlignmentFile(output_bam_path, "wb", header=input_bam.header)
    
    # 遍历BAM文件中的每条记录
    for read in input_bam:
        # 添加条形码标签
        read.set_tag("CB", barcode, "Z")  # CB是10X Genomics使用的条形码标签
        
        # 将修改后的记录写入输出文件
        output_bam.write(read)
    
    # 关闭文件
    input_bam.close()
    output_bam.close()
    
    print(f"已将条形码 {barcode} 添加到 {output_bam_path}")

def process_cell_folder(cell_folder, cell_id, barcode, temp_dir=None):
    """
    处理单个细胞文件夹，为其中的BAM文件添加条形码
    
    Args:
        cell_folder: 细胞文件夹路径
        cell_id: 细胞ID
        barcode: 要添加的条形码
        temp_dir: 临时文件目录，如果为None则使用原目录
    """
    # 查找BAM文件
    bam_path = os.path.join(cell_folder, "m_dedup.bam")
    if not os.path.exists(bam_path):
        print(f"警告：细胞 {cell_id} 的BAM文件不存在: {bam_path}")
        return
    
    # 检查索引文件
    bai_path = os.path.join(cell_folder, "m_dedup.bam.bai")
    if not os.path.exists(bai_path):
        print(f"警告：细胞 {cell_id} 的BAM索引文件不存在: {bai_path}")
    
    # 设置输出路径
    if temp_dir:
        os.makedirs(temp_dir, exist_ok=True)
        output_bam_path = os.path.join(temp_dir, f"{cell_id}.bam")
    else:
        output_bam_path = os.path.join(cell_folder, "m_dedup_with_barcode.bam")
    
    # 添加条形码
    add_barcode_to_bam(bam_path, output_bam_path, barcode)
    
    # 创建索引
    pysam.index(output_bam_path)
    print(f"已为 {output_bam_path} 创建索引")

def add_barcodes_command(args):
    """
    添加条形码命令的实现
    """
    # 读取细胞坐标信息
    cell_csv_path = os.path.join(args.data_dir, "cell.csv")
    if not os.path.exists(cell_csv_path):
        print(f"错误：细胞坐标文件不存在: {cell_csv_path}")
        return
    
    cell_coords = read_cell_coordinates(cell_csv_path)
    print(f"已读取 {len(cell_coords)} 个细胞坐标信息")
    
    # 读取条形码信息
    if not os.path.exists(args.barcode_file):
        print(f"错误：条形码文件不存在: {args.barcode_file}")
        return
    
    barcodes = read_barcodes(args.barcode_file)
    print(f"已读取 {len(barcodes)} 个条形码信息")
    
    # 处理每个细胞文件夹
    processed_count = 0
    for cell_id, coords in cell_coords.items():
        # 检查是否有对应的条形码
        if coords not in barcodes:
            print(f"警告：细胞 {cell_id} 的坐标 {coords} 没有对应的条形码，跳过处理")
            continue
        
        barcode = barcodes[coords]
        
        # 处理细胞文件夹
        cell_folder = os.path.join(args.data_dir, cell_id)
        if not os.path.exists(cell_folder):
            print(f"警告：细胞文件夹不存在: {cell_folder}")
            continue
        
        process_cell_folder(cell_folder, cell_id, barcode, args.output_dir)
        processed_count += 1
    
    print(f"处理完成，共处理了 {processed_count} 个细胞的BAM文件")

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

def merge_command(args):
    """
    合并BAM文件命令的实现
    """
    merge_barcoded_bam_files(args.input_dir, args.output_bam)

def add_and_merge_command(args):
    """
    先添加条形码然后合并BAM文件命令的实现
    """
    # 首先添加条形码
    add_barcodes_command(args)
    
    # 然后合并文件
    # 如果指定了输出目录，则使用该目录作为合并输入
    input_dir = args.output_dir if args.output_dir else os.path.join(args.data_dir, "with_barcodes")
    merge_barcoded_bam_files(input_dir, args.output_bam)

def main():
    # 创建主解析器
    parser = argparse.ArgumentParser(description="BAM文件条形码工具集")
    subparsers = parser.add_subparsers(dest="command", help="子命令")
    
    # 添加条形码子命令
    add_parser = subparsers.add_parser("add", help="将10X Genomics格式的条形码添加到BAM文件中")
    add_parser.add_argument("--data_dir", type=str, required=True, help="数据目录，包含细胞文件夹")
    add_parser.add_argument("--barcode_file", type=str, required=True, help="10X Genomics格式的条形码文件路径")
    add_parser.add_argument("--output_dir", type=str, help="输出目录，如果不指定则在原目录中创建新文件")
    
    # 合并BAM文件子命令
    merge_parser = subparsers.add_parser("merge", help="将已添加条形码的多个BAM文件合并为一个")
    merge_parser.add_argument("--input_dir", type=str, required=True, help="包含已添加条形码的BAM文件的目录")
    merge_parser.add_argument("--output_bam", type=str, required=True, help="输出的合并BAM文件路径")
    
    # 添加条形码并合并子命令
    add_merge_parser = subparsers.add_parser("add-merge", help="先添加条形码然后合并BAM文件")
    add_merge_parser.add_argument("--data_dir", type=str, required=True, help="数据目录，包含细胞文件夹")
    add_merge_parser.add_argument("--barcode_file", type=str, required=True, help="10X Genomics格式的条形码文件路径")
    add_merge_parser.add_argument("--output_dir", type=str, help="添加条形码的输出目录，如果不指定则在原目录中创建新文件")
    add_merge_parser.add_argument("--output_bam", type=str, required=True, help="最终合并后的BAM文件路径")
    
    args = parser.parse_args()
    
    # 根据子命令执行相应的功能
    if args.command == "add":
        add_barcodes_command(args)
    elif args.command == "merge":
        merge_command(args)
    elif args.command == "add-merge":
        add_and_merge_command(args)
    else:
        parser.print_help()
        sys.exit(1)

if __name__ == "__main__":
    main()