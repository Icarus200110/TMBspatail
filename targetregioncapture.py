# -*- coding: UTF-8 -*-

def loadtargetregioncapture(hgbed , targetregioncapturebed):
    bed_object = open(hgbed)
    try:
         hg_19 = bed_object.read( )
    finally:
         bed_object.close( )
    bed_object = open(targetregioncapturebed)
    try:
         target_region_capture = bed_object.read( )
    finally:
         bed_object.close( )
    hg_19_rows = hg_19.split('\n')
    target_region_capture_rows = target_region_capture.split('\n')
    #print (hg_19_rows)
    #print (target_region_capture_rows)
    target_region_list = [ ]          #目标捕获区域放在一个每项都是一个字典的列表中
    target_list_begin = 1          #目标捕获区域列表的起始位置
    target_region_size = 0          #目标捕获区域的大小初始化为０，注意并非指列表的元素个数
    for target_region_capture_row in target_region_capture_rows:
        every_row_cols = target_region_capture_row.split()
        #print (every_row_cols)
        #print target_region_capture_rows.index(target_region_capture_row)
        if len(every_row_cols) < 3:
            continue
        if len(every_row_cols[0])>3:
            chro = every_row_cols[0]
        else:
            chro = "chr" + every_row_cols[0]
        if hg_19.find(chro) == -1:
            print("目标捕获区域该值有误!")
            continue
        for hg_19_row in hg_19_rows:          #该for循环得到某一条染色体长度
            every_hgrow_cols = hg_19_row.split( )
            if every_hgrow_cols[0] == chro:
                chro_len = int(every_hgrow_cols[2])          #chro_len指hg19中某一条染色体的长度
        #print (chro_len)
        beg = int(every_row_cols[1])
        en = int(every_row_cols[2])
        #begin = beg-(en-beg)                            #为减少数据量改动了右式begin = beg-(en-beg)*50
        #end = en+(en-beg)                               #为减少数据量改动了右式end = en+(en-beg)*50
        begin = beg-(en-beg)*50
        end = en+(en-beg)*50
        #print (target_region_capture_row)
        #print (begin)
        #print (end)
        begin = max(1, min(chro_len, begin))
        end = max(1, min(chro_len, end))
        length = abs(end - begin) + 1
        #print length
        contigname = "undefined"
        if len(every_row_cols) >= 4:
            contigname = every_row_cols[3]
        #print contigname
        target_region_size = target_list_begin + length - 1
        #print target_region_size
        onechrdict_in_list = {
            'chro' : chro ,
            'contigname' : contigname ,
            'begin' : begin ,
            'end' : end ,
            'target_list_begin' : target_list_begin ,
            'target_list_end' : target_region_size
            }
        #print onechrdict_in_list
        target_region_list.append(onechrdict_in_list)
        target_list_begin += length
    return target_region_size , target_region_list
