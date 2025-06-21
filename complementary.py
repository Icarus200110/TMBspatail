# -*- coding: UTF-8 -*-

def complementary(seq_to_complement):
    complementary_dict = {
        'A': 'T',
        'C': 'G',
        'T': 'A',
        'G': 'C',
        'N': 'N'
    }
    have_complemented = ""
    for i in seq_to_complement:
        have_complemented = have_complemented + complementary_dict[i]
    return have_complemented