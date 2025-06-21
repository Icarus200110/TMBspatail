# -*- coding: UTF-8 -*-
import random


def insert_making(fragment_length):
    fragment_array = []
    for i in range(fragment_length):
        base_array = ['A', 'T', 'C', 'G']
        rand_base = int(round(random.random() * 3))
        fragment_array.append(base_array[rand_base])
    fragment = ''.join(fragment_array)
    return fragment
