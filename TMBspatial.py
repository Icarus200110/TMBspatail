import argparse
import os

from re import I
import time

import numpy as np

import cv2

import matplotlib.pyplot as plt

from pandas.compat.pyarrow import pa

import TMBspatial_cell

from Wangshj import generate_simdata

from merge_fastq import merge_fastq_files



parser  = argparse.ArgumentParser(description='Generate filled matrix with different shapes or from an image.')

parser.add_argument('--rows', type=int, required=False, help='Number of rows in the 2D coordinate system.',default=50)

parser.add_argument('--cols', type=int, required=False, help='Number of columns in the 2D coordinate system.',default=50)

parser.add_argument('--shape', type=str, choices=['Circle','Square', 'image'], required=False,

                    help='Organization shape',default='image')

parser.add_argument('--image_path', type=str, help='Path to the image.',default='./data/img.png')

parser.add_argument('--out_path', type=str, help='Path to simulating data.',default=None)

# 添加barcode文件参数

parser.add_argument('--barcode_file', type=str, help='Path to the spatial barcode file.',default='spatial_barcodes_new.txt')




parser.add_argument('--TMBspatial_snv', help='the snv, example:-snv onesnv,chr12,25398284,C,T,0.001', required=False, type=str)

parser.add_argument('--TMBspatial_cnv', help='the cnv', required=False, type=str)

parser.add_argument('--TMBspatial_mutation_mode', help='TMBspatial_mutation_mode', required=False, type=str,choices=['boundary','gradient','nested'],default='gradient')

parser.add_argument('--gradient_direction', help='the gradient direction', required=False, type=str,choices=['horizontal','vertical','radial'],default='radial')




args, unknown = parser.parse_known_args()





def generate_circle(rows, cols, shrink_ratio=0):

    """

    :param rows: 矩阵的行数

    :param cols: 矩阵的列数

    :param shrink_ratio: 缩小比例，默认为 0.15

    :return: 填充圆形矩阵

    """

    center_row = rows // 2

    center_col = cols // 2

    radius = min(center_row, center_col) * (1 - shrink_ratio)

    filled_matrix = np.zeros((rows, cols), dtype=np.float32)

    for i in range(rows):

        for j in range(cols):

            distance = np.sqrt((i - center_row) ** 2 + (j - center_col) ** 2)

            if distance <= radius:

                filled_matrix[i, j] = 1

    return filled_matrix



def generate_square(rows, cols, shrink_ratio=0):

    """

    :param rows: 矩阵的行数

    :param cols: 矩阵的列数

    :param shrink_ratio: 缩小比例，默认为 0.15

    :return: 填充正方形矩阵

    """

    center_row = rows // 2

    center_col = cols // 2

    side_length = min(rows, cols) // 2 * (1 - shrink_ratio)

    filled_matrix = np.zeros((rows, cols), dtype=int)

    top = int(center_row - side_length)

    bottom = int(center_row + side_length)

    left = int(center_col - side_length)

    right = int(center_col + side_length)

    filled_matrix[top:bottom + 1, left:right + 1] = 1

    return filled_matrix



def generate_from_image(rows, cols, image_path, shrink_ratio=0):

    """

    :param rows: 矩阵的行数

    :param cols: 矩阵的列数

    :param image_path: 图片的路径

    :param shrink_ratio: 缩小比例，默认为 0.15

    :return: 填充后的矩阵

    """

    try:

        # 读取图片

        image = cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)

        if image is None:

            raise ValueError("无法读取图片，请检查图片路径。")

        # 调整图片大小

        resized_image = cv2.resize(image, (cols, rows))

        # 使用Canny边缘检测

        edges = cv2.Canny(resized_image, 50, 150)


        # 进行形态学膨胀操作，连接断开的边缘

        kernel = np.ones((3, 3), np.uint8)

        dilated_edges = cv2.dilate(edges, kernel, iterations=1)


        # 查找轮廓

        contours, _ = cv2.findContours(dilated_edges.copy(), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)


        # 计算缩放后的轮廓

        M = cv2.moments(contours[0])

        if M["m00"] != 0:

            cX = int(M["m10"] / M["m00"])

            cY = int(M["m01"] / M["m00"])
        else:

            cX, cY = 0, 0


        scaled_contours = []

        for contour in contours:

            scaled_contour = []

            for point in contour:

                x = point[0][0]

                y = point[0][1]

                new_x = int(cX + (x - cX) * (1 - shrink_ratio))

                new_y = int(cY + (y - cY) * (1 - shrink_ratio))

                scaled_contour.append([[new_x, new_y]])

            scaled_contours.append(np.array(scaled_contour))


        filled_matrix = np.zeros((rows, cols), dtype=int)
        # 填充轮廓
        cv2.drawContours(filled_matrix, scaled_contours, -1, 1, thickness=cv2.FILLED)
        return filled_matrix

    except Exception as e:

        print(f"读取图片时出错: {e}")

        return np.zeros((rows, cols), dtype=int)


def get_all_coordinates(matrix):

    """

    读取矩阵中值为 1 的所有坐标

    :param matrix: 输入的矩阵

    :return: 包含值为 1 的坐标的列表，每个坐标以 (行索引, 列索引) 的形式表示

    """

    coordinates = []

    rows, cols = matrix.shape

    for i in range(rows):

        for j in range(cols):

            if matrix[i, j] == 1:

                coordinates.append((i+1, j+1))
    return coordinates


def visualize_matrix(matrix, output_path=None):

    """

    可视化矩阵

    :param matrix: 矩阵

    :param save_image: 是否保存图像

    :param output_path: 图像保存路径

    """

    plt.figure(figsize=(10, 8))

    plt.imshow(matrix, cmap='Blues', interpolation='nearest', vmin=0, vmax=np.max(matrix))

    #plt.imshow(matrix, cmap='Blues', interpolation='bicubic', vmin=0, vmax=np.max(matrix))

    if output_path is not None:

        name = output_path.split('/')[-1]

    plt.colorbar(label='SNP Frequency')

    plt.xlabel('y')

    plt.ylabel('x')

    plt.axis('on')
    

    if output_path is not None:

        plt.savefig(output_path)

        print(f"图像已保存至: {output_path}")
    else:

        plt.show()


def generate_matrix():

    # 获取用户输入的二维坐标系大小

    rows = args.rows

    cols = args.cols

    matrix = None
    

    if args.shape =='circle':

        matrix = generate_circle(rows, cols)

    elif args.shape == 'square':

        matrix = generate_square(rows, cols)

    elif args.shape == 'image':

        image_path = args.image_path

        matrix = generate_from_image(rows, cols, image_path)
    

    if matrix is None:

        # 如果matrix未被赋值，返回空矩阵

        matrix = np.zeros((rows, cols), dtype=int)
        

    return matrix


def boundary_mutation(matrix, boundary_type='diagonal', normal_copy=2, variant_copy=4):

    """

    分界变异模式

    :param matrix: 原始细胞矩阵

    :param boundary_type: 分界类型 ('diagonal'-对角线, 'horizontal'-水平, 'vertical'-垂直)

    :param normal_copy: 正常区域的拷贝数 (默认为2)

    :param variant_copy: 变异区域的拷贝数 (默认为4)

    :return: 拷贝数变异矩阵

    """

    rows, cols = matrix.shape

    mutation_matrix = np.ones_like(matrix, dtype=float) * normal_copy
    

    for i in range(rows):

        for j in range(cols):

            if boundary_type == 'diagonal':

                if i > j:  # 对角线分界

                    mutation_matrix[i,j] = variant_copy

            elif boundary_type == 'horizontal':

                if i > rows//2:  # 水平分界

                    mutation_matrix[i,j] = variant_copy

            elif boundary_type == 'vertical':

                if j > cols//2:  # 垂直分界

                    mutation_matrix[i,j] = variant_copy
    

    return np.where(matrix == 1, mutation_matrix, 0)


def gradient_mutation(matrix, gradient_direction='horizontal', max_copy=9, center_row=None, center_col=None, min_copy=0):

    """

    渐变变异模式

    :param matrix: 原始细胞矩阵

    :param gradient_direction: 渐变方向 ('horizontal'-水平, 'vertical'-垂直, 'radial'-圆心扩散)

    :param max_copy: 最大拷贝数 (默认为9)

    :param center_row: 中心行坐标 (None表示自动计算中心)

    :param center_col: 中心列坐标 (None表示自动计算中心)

    :param min_copy: 最小拷贝数 (默认为0)

    :return: 拷贝数变异矩阵 (值大于min_copy且保留一位小数)

    """

    max_copy = float(max_copy)

    rows, cols = matrix.shape

    mutation_matrix = np.ones_like(matrix, dtype=float) * min_copy
    

    # 设置中心点

    if center_row is None:

        center_row = rows // 2

    if center_col is None:

        center_col = cols // 2
        

    max_radius = min(center_row, center_col, rows-center_row, cols-center_col)
    

    for i in range(rows):

        for j in range(cols):

            if gradient_direction == 'horizontal':

                # 水平渐变 (从左到右，拷贝数从min_copy增加到max_copy)

                copy_num = min_copy + (max_copy - min_copy) * j/cols

                mutation_matrix[i,j] = round(max(min_copy + 0.1, copy_num), 1)

            elif gradient_direction == 'vertical':

                # 垂直渐变 (从上到下，拷贝数从min_copy增加到max_copy)

                copy_num = min_copy + (max_copy - min_copy) * i/rows

                mutation_matrix[i,j] = round(max(min_copy + 0.1, copy_num), 1)

            elif gradient_direction == 'radial':

                # 圆心扩散 (从中心向外，拷贝数从max_copy减少到min_copy)

                distance = np.sqrt((i - center_row)**2 + (j - center_col)**2)

                copy_num = max_copy - (max_copy - min_copy) * (distance/max_radius)

                mutation_matrix[i,j] = round(max(min_copy + 0.1, copy_num), 1)
    

    return np.where(matrix == 1, mutation_matrix, 0)


def nested_mutation(matrix, inner_copy=4, middle_copy=1, outer_copy=2, center_row=None, center_col=None):

    """

    嵌套变异模式

    :param matrix: 原始细胞矩阵

    :param inner_copy: 核心区域拷贝数 (默认为4)

    :param middle_copy: 中间区域拷贝数 (默认为1)

    :param outer_copy: 外围区域拷贝数 (默认为2)

    :param center_row: 中心行坐标 (None表示自动计算中心)

    :param center_col: 中心列坐标 (None表示自动计算中心)

    :return: 拷贝数变异矩阵

    """

    rows, cols = matrix.shape

    mutation_matrix = np.ones_like(matrix, dtype=float) * outer_copy
    

    # 设置中心点

    if center_row is None:

        center_row = rows // 2

    if center_col is None:

        center_col = cols // 2
        

    max_radius = min(center_row, center_col, rows-center_row, cols-center_col)
    

    for i in range(rows):

        for j in range(cols):

            distance = np.sqrt((i - center_row)**2 + (j - center_col)**2)

            # 根据距离中心的远近设置不同的拷贝数

            if distance < max_radius * 0.3:

                mutation_matrix[i,j] = inner_copy  # 核心区域拷贝数

            elif distance < max_radius * 0.6:

                mutation_matrix[i,j] = middle_copy  # 中间区域拷贝数
    

    return np.where(matrix == 1, mutation_matrix, 0)


def mutation_matrix(matrix, mode='boundary', boundary_type='diagonal', gradient_direction='radial', max_copy=6):

    """

    生成变异矩阵 (入口函数)

    :param matrix: 原始细胞矩阵

    :param mode: 变异模式 ('boundary'-分界, 'gradient'-渐变, 'nested'-嵌套)

    :param boundary_type: 分界类型 ('diagonal'-对角线, 'horizontal'-水平, 'vertical'-垂直)

    :param gradient_direction: 渐变方向 ('horizontal'-水平, 'vertical'-垂直, 'radial'-圆心扩散)

    :param max_copy: 最大拷贝数 (默认为4)

    :return: 拷贝数变异矩阵 (2表示正常拷贝数)

    """

    if mode == 'boundary':

        return boundary_mutation(matrix, boundary_type)

    elif mode == 'gradient':

        return gradient_mutation(matrix, gradient_direction, max_copy)

    elif mode == 'nested':

        return nested_mutation(matrix)
    else:

        return np.where(matrix == 1, 2, 0)  # 默认返回正常拷贝数



def snv_matrix(matrix, mode='boundary', boundary_type='diagonal', gradient_direction='radial', max_copy=1.0):

    """

    生成SNV变异概率矩阵 (入口函数)

    :param matrix: 原始细胞矩阵

    :param mode: 变异模式 ('boundary'-分界, 'gradient'-渐变, 'nested'-嵌套)

    :param boundary_type: 分界类型 ('diagonal'-对角线, 'horizontal'-水平, 'vertical'-垂直)

    :param gradient_direction: 渐变方向 ('horizontal'-水平, 'vertical'-垂直, 'radial'-圆心扩散)

    :param max_copy: 最大变异概率 (默认为1.0)

    :return: SNV变异概率矩阵 (值范围0-1)

    """

    # 使用mutation_matrix生成基础矩阵

    prob_matrix = mutation_matrix(matrix, mode=mode, boundary_type=boundary_type, 

                               gradient_direction=gradient_direction, max_copy=max_copy)
    

    # 将拷贝数转换为变异概率 (归一化到0-1范围)

    # 对于boundary模式，将4转换为max_copy，将2转换为0

    if mode == 'boundary':

        for i in range(prob_matrix.shape[0]):

            for j in range(prob_matrix.shape[1]):

                if prob_matrix[i,j] == 4:  # 变异区域

                    prob_matrix[i,j] = max_copy

                elif prob_matrix[i,j] == 2:  # 正常区域

                    prob_matrix[i,j] = 0

    # 对于nested模式，将middle_copy(1)转换为0（正常区域），将inner_copy(4)和outer_copy(2)转换为变异概率

    elif mode == 'nested':

        for i in range(prob_matrix.shape[0]):

            for j in range(prob_matrix.shape[1]):

                if prob_matrix[i,j] == 1:  # 中间区域（正常区域）

                    prob_matrix[i,j] = 0

                elif prob_matrix[i,j] == 4:  # 内部区域（变异区域）

                    prob_matrix[i,j] = max_copy

                elif prob_matrix[i,j] == 2:  # 外部区域（变异区域）

                    prob_matrix[i,j] = max_copy * 0.5  # 外部区域变异概率设为最大值的一半
    

    return prob_matrix


def generate_organization(matrix):

    """

    根据矩阵生成组织，可选包含拷贝数变异信息

    :param matrix: 原始细胞矩阵

    :param mutation_matrix:

    :return: Organization对象

    """

    organization = TMBspatial_cell.Organization(matrix.shape[0], matrix.shape[1])

    all_coordinates = get_all_coordinates(matrix)

    i = 0

    a_snv_matrix = None

    cnv = args.TMBspatial_cnv

    snv = args.TMBspatial_snv


    if snv:

        snv_arr = snv.split(",")

        a_snv_matrix = snv_matrix(matrix, mode=args.TMBspatial_mutation_mode, gradient_direction=args.gradient_direction, max_copy=float(snv_arr[5]) if len(snv_arr) > 5 else 1.0)
    


    if cnv:

        cnv_arr = cnv.split(":")

        cnv_matrix_arr = []

        for cnv_one in cnv_arr:

            cnv_one_arr = cnv_one.split(",")

            cnv_matrix = mutation_matrix(matrix, mode=args.TMBspatial_mutation_mode,gradient_direction=args.gradient_direction,max_copy=cnv_one_arr[4])

            cnv_matrix_arr.append(cnv_matrix)
    else:

        cnv_matrix_arr = None
    
    

    for coordinate in all_coordinates:

        id = str(i).zfill(8)


        # 创建细胞

        cell = TMBspatial_cell.Cell(id, coordinate)

        #snv变异

        if a_snv_matrix is not None:

            snv_arr[5] = str(a_snv_matrix[coordinate[0]-1, coordinate[1]-1])

            asnv=",".join(snv_arr)

            cell.set_snv(asnv)


        # 使用set_cnv方法设置拷贝数

        if cnv_matrix_arr is not None:

            cell_cnv = ''

            for (cnv_matrix,acnv) in zip(cnv_matrix_arr,cnv_arr):

                # 获取拷贝数，默认为2

                copy_number = 2

                # 注意坐标从1开始，而矩阵索引从0开始

                copy_number = cnv_matrix[coordinate[0]-1, coordinate[1]-1]

                cnv_one_arr = acnv.split(",")

                cnv_one_arr[4] = str(copy_number)

                cell_cnv += ",".join(cnv_one_arr)+':'

            cell_cnv = cell_cnv[:-1]

            cell.set_cnv(cell_cnv)
  

        organization.add({id:cell})

        i += 1

        organization.set_matrix(matrix)

    return organization


def read_barcodes(barcode_file):

    """

    读取barcode文件

    :param barcode_file: barcode文件路径

    :return: 坐标到barcode的映射字典

    """

    barcodes = {}

    if barcode_file and os.path.exists(barcode_file):

        try:

            with open(barcode_file, 'r') as f:

                for line in f:

                    parts = line.strip().split('\t')

                    if len(parts) >= 3:

                        barcode = parts[0]

                        row = int(parts[1])

                        col = int(parts[2])

                        barcodes[(row, col)] = barcode

            print(f"成功读取{len(barcodes)}个barcode")

        except Exception as e:

            print(f"读取barcode文件出错: {e}")

    return barcodes


def insert_barcode_to_make2(out_dir, cell_coord, barcodes):

    """

    在make2.fq文件中插入barcode

    :param out_dir: 输出目录

    :param cell_coord: 细胞坐标

    :param barcodes: 坐标到barcode的映射字典

    """

    make2_path = os.path.join(out_dir, 'make2.fq')

    if not os.path.exists(make2_path):

        print(f"make2.fq文件不存在，跳过barcode插入: {make2_path}")
        return
    

    # 如果没有对应的barcode，直接返回

    if cell_coord not in barcodes:

        print(f"坐标 {cell_coord} 没有对应的barcode，跳过插入")
        return
    

    barcode = barcodes[cell_coord]

    if len(barcode) < 16:

        print(f"barcode长度不足16位，跳过插入: {barcode}")
        return
    

    # 读取文件内容

    with open(make2_path, 'r') as f:

        lines = f.readlines()
    

    # 只修改第一个序列，每4行为一组（标题行、序列行、加号行、质量行）

    if len(lines) >= 4:

        header = lines[0]

        sequence = lines[1].strip()

        plus = lines[2]

        quality = lines[3].strip()
        

        # 插入barcode

        barcode1 = barcode[8:16]  # 取barcode的9-16位作为Barcode1

        barcode2 = barcode[0:8]   # 取barcode的1-8位作为Barcode2
        

        # 确保序列长度足够

        if len(sequence) >= 40:

            # 在33-40位置插入Barcode2

            sequence = sequence[:32] + barcode2 + sequence[40:]
            

            # 调整质量值

            #quality = quality[:32] + 'I' * 8 + quality[40:]
            

            # 如果序列长度足够，在71-78位置插入Barcode1

            if len(sequence) >= 78:

                sequence = sequence[:70] + barcode1 + sequence[78:]

                #quality = quality[:70] + 'I' * 8 + quality[78:]
        

        # 更新文件内容

        lines[1] = sequence + '\n'

        #lines[3] = quality + '\n'
        

        # 写回文件

        with open(make2_path, 'w') as f:

            f.writelines(lines)
        

        print(f"barcode插入完成: {make2_path}")
    else:

        print(f"文件格式不正确，跳过barcode插入: {make2_path}")


def out_TMBspatial_data(organization):

    out = args.out_path

    #用户为指定数据输出路径则使用当前时间做文件夹名保存文件

    if out is None:

        out = time.strftime('%Y%m%d%H%M%S', time.localtime(time.time()))

    if not os.path.exists(out):

        os.makedirs(out)

    organization.save_to_csv(out)

    matrix = organization.matrix

    # 保存原始组织形态的可视化

    visualize_matrix(matrix, out+'/original_matrix.png')

    a_snv_matrix = None

    cnv = args.TMBspatial_cnv

    snv = args.TMBspatial_snv


    if snv:

        snv_arr = snv.split(",")

        a_snv_matrix = snv_matrix(matrix, mode=args.TMBspatial_mutation_mode, gradient_direction=args.gradient_direction, max_copy=float(snv_arr[5]) if len(snv_arr) > 5 else 1.0)

        visualize_matrix(a_snv_matrix, out+'/snv_matrix.png')
    

    if cnv:

        cnv_arr = cnv.split(":")

        i = 0

        for cnv_one in cnv_arr:

            cnv_one_arr = cnv_one.split(",")

            cnv_matrix = mutation_matrix(matrix, mode=args.TMBspatial_mutation_mode,gradient_direction=args.gradient_direction,max_copy=cnv_one_arr[4])

            visualize_matrix(cnv_matrix, out+f'/cnv_matrix{i}.png')

            i += 1

    
    

    #读取barcode文件

    barcodes = read_barcodes(args.barcode_file)

    # 为每个细胞生成数据并插入barcode

    for v in organization.cell_dict.values():

        # 生成模拟数据

        generate_simdata(v.id, out, TMBspatial_snv=v.snv,TMBspatial_cnv=v.cnv)
    

        # 如果有barcode文件，则插入barcode

        if args.barcode_file:

            insert_barcode_to_make2(out+'/'+v.id, v.coordinate, barcodes)
    

    merge_fastq_files(out,'merge')
    
    return




def main():

    #生成组织矩阵

    #获得细胞坐标

    matrix = generate_matrix()

    organization = generate_organization(matrix)

    out_TMBspatial_data(organization)



if __name__ == "__main__":
    main()