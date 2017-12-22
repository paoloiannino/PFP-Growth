# encoding utf-8
import random

MAX_ID = 128;
MAX_ITEMS = 10;
N_TRASACTIONS = 10000

file_csv = open('./data/artificial_dataset_pfp.csv','w')
file_tsv = open('./data/artificial_dataset_pfp.tsv','w')

for i in range(0, N_TRASACTIONS):
    n_trans = random.sample(range(1, MAX_ITEMS), 1)[0]
    list_trans = random.sample(range(1, MAX_ID), n_trans)
    csv_string = str()
    tsv_string = str()
    for num in list_trans:
        csv_string = csv_string + ',' + str(num)
        tsv_string = tsv_string + '\t' + str(num)

    csv_string = csv_string + '\n'
    tsv_string = tsv_string + '\n'
    print i;

    file_csv.write(csv_string)
    file_tsv.write(tsv_string)

file_csv.close()
file_tsv.close()
