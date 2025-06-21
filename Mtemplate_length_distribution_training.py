# coding:utf-8
import os
import re
import argparse
import sys
from sklearn.externals import joblib
reload(sys)
sys.setdefaultencoding( "utf-8" )
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures

parser = argparse.ArgumentParser()
parser.add_argument('-csv', '--csv', help='the input_insert_csv_file_path', required=True, type=str)
parser.add_argument('-out', '--outpath', help='the out_template_length_distribution__model_path', required=False, type=str)
args = parser.parse_args()
csv = args.csv
outpath = args.outpath

outdir = re.search(r'[^/]+$', csv).group(0)[:-4]
outdir_arr = outdir.split('.')
outdir = ""
for outdir_part in outdir_arr:
    outdir = outdir + "_" + outdir_part
outdir = "templateLengthDistribution" + outdir
if outpath != None:
        out = outpath
        out = out.strip()
        out = out.rstrip("\/")
        isExists=os.path.exists(out)
        if not isExists:
            os.makedirs(out)
else:
    out = outdir
    os.makedirs(out)
out = out + '/'

#read_csv = open("/mnt/X500/farmers/prod/runOut/OncoH/output/160002932BD/normal/5_recal_bam/160002932BD.insert.csv", 'r')
read_csv = open(csv, 'r')
try:
    csv_file = read_csv.read()
finally:
    read_csv.close()
file_array = csv_file.split('\n')
if not file_array[0].split(',')[0].isdigit() or not file_array[0].split(',')[1].isdigit():
    del file_array[0]
if file_array[-1] == '':
    file_array.pop()

x = []
y = []
xx = []
yy = []
xxx = []
yyy = []
first_one = int(file_array[0].split(',')[0])
last_one = int(file_array[-1].split(',')[0])
last_sub_first = float(last_one - first_one + 1)
for line in file_array:
    linearray = line.split(",")
    linearray0 = int(float(int(linearray[0]) - first_one + 1) / last_sub_first * 1000)
    x.append(linearray0)
    y.append(int(linearray[1]))
for i in range(len(x) - 3):
    x1 = x[i]
    y1 = y[i]
    x2 = x[i + 1]
    y2 = y[i + 1]
    x3 = x[i + 2]
    y3 = y[i + 2]
    x4 = x[i + 3]
    y4 = y[i + 3]
    A = y2 - y1
    B = x1 - x2
    C = x2 * y1 - x1 * y2
    flag1 = A * x3 + B * y3 + C
    flag2 = A * x4 + B * y4 + C
    flag = flag1 * flag2
    if flag < 0:
        xx.append(x3)
        yy.append(y3)

for i in range(len(xx) - 3):
    x1 = xx[i]
    y1 = yy[i]
    x2 = xx[i + 1]
    y2 = yy[i + 1]
    x3 = xx[i + 2]
    y3 = yy[i + 2]
    x4 = xx[i + 3]
    y4 = yy[i + 3]
    A = y2 - y1
    B = x1 - x2
    C = x2 * y1 - x1 * y2
    flag1 = A * x3 + B * y3 + C
    flag2 = A * x4 + B * y4 + C
    flag = flag1 * flag2
    if flag < 0:
        xxx.append(x3)
        yyy.append(y3)

inflectionPoint1 = xxx[0]
inflectionPoint2 = xxx[0]
for xxx_one in xxx:
    if abs(xxx_one - 166) < abs(inflectionPoint1 - 166):
        inflectionPoint1 = xxx_one
    if abs(xxx_one - 550) < abs(inflectionPoint2 - 550):
        inflectionPoint2 = xxx_one

write_parameter = open(out + "parameter.txt", "ab")
to_write_parameter = str(inflectionPoint1 + 1) + "\t" + str(inflectionPoint2 + 1) + "\n"
write_parameter.write(to_write_parameter)
write_parameter.close()

for i in range(len(x)):
    if x[i] == inflectionPoint1:
        inflectionPoint1_index = i
for i in range(len(x)):
    if x[i] == inflectionPoint2:
        inflectionPoint2_index = i

initial_X = []
initial_y = []
for line in file_array:
    linearray = line.split(",")
    X_one = []
    y_one = []
    linearray0 = int(float(int(linearray[0]) - first_one + 1) / last_sub_first * 1000)
    X_one.append(linearray0)
    y_one.append(int(linearray[1]))
    initial_X.append(X_one)
    initial_y.append(y_one)

X = initial_X[0:inflectionPoint1_index]
y = initial_y[0:inflectionPoint1_index]
X_test = X
y_test = y
xx = np.linspace(1, 1000, 1000)
perfect = 0
power = 1
for i in range(20):
    cubic_featurizer = PolynomialFeatures(degree=i + 1)
    X_train_cubic = cubic_featurizer.fit_transform(X)
    regressor_cubic = LinearRegression()
    regressor_cubic.fit(X_train_cubic, y)
    xx_cubic = cubic_featurizer.transform(xx.reshape(xx.shape[0], 1))

    X_test_cubic = cubic_featurizer.transform(X_test)
    if regressor_cubic.score(X_test_cubic, y_test) > perfect:
        perfect = regressor_cubic.score(X_test_cubic, y_test)
        power = i + 1
#print power
write_parameter = open(out + "parameter.txt", "ab")
to_write_parameter = str(power) + "\t"
write_parameter.write(to_write_parameter)
write_parameter.close()
cubic_featurizer = PolynomialFeatures(degree=power)
X_train_cubic = cubic_featurizer.fit_transform(X)
regressor_cubic = LinearRegression()
regressor_cubic.fit(X_train_cubic, y)
joblib.dump(regressor_cubic, out + "model_left.m")
xx_cubic = cubic_featurizer.transform(xx.reshape(xx.shape[0], 1))
X_test_cubic = cubic_featurizer.transform(X_test)
print regressor_cubic.score(X_test_cubic, y_test)

X = initial_X[inflectionPoint1_index:inflectionPoint2_index]
y = initial_y[inflectionPoint1_index:inflectionPoint2_index]
X_test = X
y_test = y
xx = np.linspace(1, 1000, 1000)
perfect = 0
power = 1
for i in range(20):
    cubic_featurizer = PolynomialFeatures(degree=i + 1)
    X_train_cubic = cubic_featurizer.fit_transform(X)
    regressor_cubic = LinearRegression()
    regressor_cubic.fit(X_train_cubic, y)
    xx_cubic = cubic_featurizer.transform(xx.reshape(xx.shape[0], 1))
    X_test_cubic = cubic_featurizer.transform(X_test)
    if regressor_cubic.score(X_test_cubic, y_test) > perfect:
        perfect = regressor_cubic.score(X_test_cubic, y_test)
        power = i + 1
#print power
write_parameter = open(out + "parameter.txt", "ab")
to_write_parameter = str(power) + "\t"
write_parameter.write(to_write_parameter)
write_parameter.close()
cubic_featurizer = PolynomialFeatures(degree=power)
X_train_cubic = cubic_featurizer.fit_transform(X)
regressor_cubic = LinearRegression()
regressor_cubic.fit(X_train_cubic, y)
joblib.dump(regressor_cubic, out + "model_middle.m")
xx_cubic = cubic_featurizer.transform(xx.reshape(xx.shape[0], 1))
X_test_cubic = cubic_featurizer.transform(X_test)
print regressor_cubic.score(X_test_cubic, y_test)

X = initial_X[inflectionPoint2_index:]
y = initial_y[inflectionPoint2_index:]
X_test = X
y_test = y
xx = np.linspace(1, 1000, 1000)
perfect = 0
power = 1
for i in range(20):
    cubic_featurizer = PolynomialFeatures(degree=i + 1)
    X_train_cubic = cubic_featurizer.fit_transform(X)
    regressor_cubic = LinearRegression()
    regressor_cubic.fit(X_train_cubic, y)
    xx_cubic = cubic_featurizer.transform(xx.reshape(xx.shape[0], 1))
    X_test_cubic = cubic_featurizer.transform(X_test)
    if regressor_cubic.score(X_test_cubic, y_test) > perfect:
        perfect = regressor_cubic.score(X_test_cubic, y_test)
        power = i + 1
#print power
write_parameter = open(out + "parameter.txt", "ab")
to_write_parameter = str(power)
write_parameter.write(to_write_parameter)
write_parameter.close()
cubic_featurizer = PolynomialFeatures(degree=power)
X_train_cubic = cubic_featurizer.fit_transform(X)
regressor_cubic = LinearRegression()
regressor_cubic.fit(X_train_cubic, y)
joblib.dump(regressor_cubic, out + "model_right.m")
xx_cubic = cubic_featurizer.transform(xx.reshape(xx.shape[0], 1))
X_test_cubic = cubic_featurizer.transform(X_test)
print regressor_cubic.score(X_test_cubic, y_test)
