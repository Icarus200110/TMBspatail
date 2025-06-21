# -*- coding: UTF-8 -*-
import json

def load_json_file(json_file):
    with open(json_file, "rb") as load_f:  # 打开json文件并转化为字典格式
         load_dict = json.load(load_f)
    load_f.close()
    return load_dict








