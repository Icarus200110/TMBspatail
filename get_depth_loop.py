import random


def get_depth(chromo, position, depth_dict):
    chromo = chromo[3:]
    position = str(position)
    try:
        depth = depth_dict[chromo][position]
        return depth
    except:
        return 0


def create_loop(chromo, position, length, depth_dict):
    depth_num = 0
    depth_all = 0
    for i in range(length):
        depth_one = get_depth(chromo, position + i, depth_dict)
        depth_all = depth_all + depth_one
        if depth_one != 0:
            depth_num = depth_num + 1
    #print depth_num, depth_all
    if depth_num == 0:
        return 1
    else:
        depth_average = depth_all / depth_num
        depth_average_floor = int(depth_average)
        depth_average_ceil = depth_average_floor + 1
        depth_average_decimal = depth_average - depth_average_floor
        #print depth_average_decimal
        depth_adjust_rdm = random.random()
        #print depth_adjust_rdm
        if depth_adjust_rdm < depth_average_decimal:
            return depth_average_ceil
        else:
            return depth_average_floor

