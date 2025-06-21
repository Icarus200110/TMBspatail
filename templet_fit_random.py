# -*- coding: UTF-8 -*
import numpy as np
import scipy as sp
import joblib
from sklearn.preprocessing import PolynomialFeatures


def make_templet_length(appro_templet_amount, templetmax, templetmin, templateLengthDistributionModel_path):
    if templateLengthDistributionModel_path == None:
        templet_length_array = []

        alphabet = []
        for i in range(1000):
            alphabet.append(i + 1)
        pmf = []
        for i in range(1, 167):
            value1 = []
            value2 = []
            value1.append(
                0.0000379718774 * (i ** 1) + 0.00118295334 * (i ** 2) + 0.0184003103 * (i ** 3) - 0.000636816316 * (
                    i ** 4) + 0.0000114074319 * (i ** 5) - 0.000000101030748 * (i ** 6) + 0.000000000431231351 * (
                    i ** 7) - 0.000000000000715259673 * (i ** 8) + 47.18052515)
            value2.append(value1)
            pmf.append(value2)
        for i in range(167, 551):
            value1 = []
            value2 = []
            value1.append(
                575.031991 * (i ** 1) - 2.31071674 * (i ** 2) + 0.00395637192 * (i ** 3) - 0.00000295084638 * (
                    i ** 4) + 0.000000000690420165 * (i ** 5) - 37722.14711989)
            value2.append(value1)
            pmf.append(value2)
        for i in range(551, 1001):
            value1 = []
            value2 = []
            value1.append(
                0.000606129888 * (i ** 3) - 0.277263744 * (i ** 2) - 0.000737390103 * (i ** 1) - 0.000000501123911 * (
                    i ** 4) + 0.000000000148277086 * (i ** 5) + 24019.54295762)
            value2.append(value1)
            pmf.append(value2)
        alphabet = np.array(alphabet)
        pmf = np.array(pmf)
        X = sp.zeros([appro_templet_amount, 1])
        for i in sp.r_[0:appro_templet_amount]:
            X[i] = GenerateSample(alphabet, pmf)
            templet_length_array.append(
                int(round(float(list(X[i])[0]) / 1000 * (templetmax - templetmin + 1) + (templetmin - 1))))
        return templet_length_array
    else:
        templet_length_array = []

        parameter_file = open(templateLengthDistributionModel_path + "/parameter.txt",'r')
        try:
            parameter_file_str = parameter_file.read()
        finally:
            parameter_file.close()
        parameter_file_arr = parameter_file_str.split('\n')
        inflectionPoint_arr = parameter_file_arr[0].split('\t')
        degree_arr = parameter_file_arr[1].split('\t')
        #print inflectionPoint_arr
        #print degree_arr

        left = joblib.load(templateLengthDistributionModel_path + "/model_left.m")
        middle = joblib.load(templateLengthDistributionModel_path + "/model_middle.m")
        right = joblib.load(templateLengthDistributionModel_path + "/model_right.m")

        left_featurizer = PolynomialFeatures(degree=int(degree_arr[0]))
        middle_featurizer = PolynomialFeatures(degree=int(degree_arr[1]))
        right_featurizer = PolynomialFeatures(degree=int(degree_arr[2]))

        alphabet = []
        for i in range(1000):
            alphabet.append(i + 1)
        pmf = []
        for i in range(1, int(inflectionPoint_arr[0])):
            argument1 = []
            argument2 = []
            argument1.append(i)
            argument2.append(argument1)
            pmf.append(left.predict(left_featurizer.fit_transform(argument2)))
        for i in range(int(inflectionPoint_arr[0]), int(inflectionPoint_arr[1])):
            argument1 = []
            argument2 = []
            argument1.append(i)
            argument2.append(argument1)
            pmf.append(middle.predict(middle_featurizer.fit_transform(argument2)))
        for i in range(int(inflectionPoint_arr[1]), 1001):
            argument1 = []
            argument2 = []
            argument1.append(i)
            argument2.append(argument1)
            pmf.append(right.predict(right_featurizer.fit_transform(argument2)))
        alphabet = np.array(alphabet)
        pmf = np.array(pmf)
        X = sp.zeros([appro_templet_amount, 1])
        for i in sp.r_[0:appro_templet_amount]:
            X[i] = GenerateSample(alphabet, pmf)
            templet_length_array.append(int(round(float(list(X[i])[0]) / 1000 * (templetmax - templetmin + 1) + (templetmin - 1))))
        return templet_length_array

def GenerateSample(alphabet, pmf):
    numSymbols = len(alphabet)
    assert(numSymbols == len(pmf))
    reference_pmf = sp.ones(numSymbols)/numSymbols
    c=max(pmf)*numSymbols
    condition_satisfied = False
    while (not condition_satisfied):
        Y = spread(alphabet)
        fY = pmf[Y==alphabet]
        gY = reference_pmf[Y==alphabet]
        U = sp.rand()
        if ( U <= fY/(c*gY)):
            X=Y
            condition_satisfied = True
        else:
            condition_satisfied = False
    return X

def spread(alphabet):
    n = len(alphabet)
    U = sp.rand()
    X = sp.int32(n*U) + 1
    return alphabet[X-1]
