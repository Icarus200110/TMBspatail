import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-in', '--input_tsv', help='the input_tsv', required=True, type=str)

args = parser.parse_args()
input_tsv = args.input_tsv

read_tsv = open(input_tsv, "rb")
try:
    tsv_file = read_tsv.read()
finally:
    read_tsv.close()
tsv_file_rows = tsv_file.split('\n')
while '' in tsv_file_rows:
    tsv_file_rows.remove('')


tsv_file_row_len = len(tsv_file_rows[0].split())
for tsv_file_row in tsv_file_rows:
    if len(tsv_file_row.split()) < tsv_file_row_len:
        print "Sorry, there is not enough data to verify. Please add the data! Or if you want to skip some of the samples, please set the value of if_calling to be N."
        exit()

if tsv_file_rows[0].split()[0].split('_')[-1] == "normal":
    if_calling = 0
    not_calling = 0
    for tsv_file_row in tsv_file_rows:
        every_tfrow_cols = tsv_file_row.split()
        if (every_tfrow_cols[-1] == '1') or (every_tfrow_cols[-1] == '0'):
            if_calling = if_calling + 1
        if every_tfrow_cols[-1] == '0':
            not_calling = not_calling + 1
    print not_calling, if_calling
    specificity = format(float(not_calling) / float(if_calling) * 100, '.4f')
    if float(specificity) >= 60:
        print "Congratulations! the specificity of your mutation calling software has reached " + specificity + "%"
    else:
        print "Sorry! the specificity of your mutation calling software is only " + specificity + "%"
else:
    if_calling = 0
    is_calling = 0
    for tsv_file_row in tsv_file_rows:
        every_tfrow_cols = tsv_file_row.split()
        if (every_tfrow_cols[-1] == '1') or (every_tfrow_cols[-1] == '0'):
            if_calling = if_calling + 1
        if every_tfrow_cols[-1] == '1':
            is_calling = is_calling + 1
    sensitivity = format(float(is_calling) / float(if_calling) * 100, '.4f')
    mutation_info = tsv_file_rows[0].split()[0]
    mutation_info_arr = mutation_info.split('_')[1:]
    gradient = mutation_info_arr.pop()
    if len(mutation_info_arr) > 1:
        mutation = "tandem_repeat"
    else:
        mutation = mutation_info_arr[0]
    if float(sensitivity) >= 60:
        print "Congratulations! In the gradient of " + gradient + ", the sensitivity of your " + mutation + " calling software has reached " + sensitivity + "%"
    else:
        print "Sorry! In the gradient of " + gradient + ", the sensitivity of your " + mutation + " calling software is only " + sensitivity + "%"

