# 将spatial_barcodes_new.txt转换为10X Genomics格式

input_file = 'd:/Project/python/GSDcreator/spatial_barcodes_new.txt'
output_file = 'd:/Project/python/GSDcreator/spatial_barcodes_10x.csv'

with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    # 写入CSV头部
    outfile.write('barcode,x,y\n')
    
    # 处理每一行数据
    for line in infile:
        line = line.strip()
        if line:
            parts = line.split('\t')
            if len(parts) == 3:
                barcode, row, col = parts
                # 10X Genomics格式：barcode,x,y
                outfile.write(f'{barcode},{col},{row}\n')

print(f'转换完成。输出已保存到 {output_file}')