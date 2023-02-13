import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import collections

#Impiort data
df = pd.read_csv('rmsd_sh2_PTR_align.csv')

#get rows and columns
num_row = len(df)
num_col = len(df.columns)

#Array for each position
rmsd_pos = np.zeros(num_row)

#
for i in range(num_row):
    raw_data_i = df.iloc[i]
    data_i = []
    for n in raw_data_i:
        if n != -1 and n != '-1':
            if isinstance(n, (collections.Sequence, np.ndarray)):
                data_i.append(float(n.strip('[]'))*10)
            else:
                data_i.append(n*10)
    data_i.pop(0)
    rmsd_pos[i] = np.mean(data_i)

position = ['-4', '-3', '-2', '-1', '0', '+1', '+2', '+3', '+4', '+5', '+6', '+7', '+8']
plt.figure()
plt.bar(position, rmsd_pos)
plt.ylabel(r'Backbone RMSD ($\AA$)')
plt.xlabel('Position Relative to pTYR')
plt.title('Relative Flexibility of Peptide Residues')
plt.savefig('sh2_rmsd_cmpr_sh2_PTR_align.png')
