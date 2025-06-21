import os
import csv
import pysam
import argparse

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

def main():
    parser = argparse.ArgumentParser(description="将10X Genomics格式的条形码添加到BAM文件中")
    parser.add_argument("--data_dir", type=str, required=True, help="数据目录，包含细胞文件夹")
    parser.add_argument("--barcode_file", type=str, required=True, help="10X Genomics格式的条形码文件路径")
    parser.add_argument("--output_dir", type=str, help="输出目录，如果不指定则在原目录中创建新文件")
    args = parser.parse_args()
    
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

if __name__ == "__main__":
    main()