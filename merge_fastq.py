import os
import glob

def merge_fastq_files(input_dir, output_prefix):
    """
    合并所有细胞文件夹中的make1.fq和make2.fq文件
    
    Args:
        input_dir: 包含所有细胞文件夹的目录
        output_prefix: 输出文件的前缀
    """
    # 确保输入目录存在
    if not os.path.exists(input_dir):
        print(f"输入目录不存在: {input_dir}")
        return
    
    # 获取所有细胞文件夹
    cell_dirs = [d for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d)) and d.startswith('0')]
    cell_dirs.sort()  # 按照名称排序
    
    # 输出文件路径
    output_make1 = os.path.join(input_dir, f"{output_prefix}_make1.fq")
    output_make2 = os.path.join(input_dir, f"{output_prefix}_make2.fq")
    
    # 合并make1.fq文件
    print(f"开始合并make1.fq文件到: {output_make1}")
    with open(output_make1, 'w') as outfile:
        for cell_dir in cell_dirs:
            make1_path = os.path.join(input_dir, cell_dir, 'make1.fq')
            if os.path.exists(make1_path):
                print(f"处理: {make1_path}")
                with open(make1_path, 'r') as infile:
                    outfile.write(infile.read())
    
    # 合并make2.fq文件
    print(f"开始合并make2.fq文件到: {output_make2}")
    with open(output_make2, 'w') as outfile:
        for cell_dir in cell_dirs:
            make2_path = os.path.join(input_dir, cell_dir, 'make2.fq')
            if os.path.exists(make2_path):
                print(f"处理: {make2_path}")
                with open(make2_path, 'r') as infile:
                    outfile.write(infile.read())
    
    print("合并完成!")

if __name__ == "__main__":
    # 指定包含所有细胞文件夹的目录
    input_directory = r"d:\Project\python\GSDcreator\20250310172408"
    
    # 指定输出文件的前缀
    output_prefix = "merged"
    
    # 执行合并
    merge_fastq_files(input_directory, output_prefix)