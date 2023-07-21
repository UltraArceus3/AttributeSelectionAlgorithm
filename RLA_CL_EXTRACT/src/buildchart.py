from cProfile import label
from matplotlib import pyplot as plt
import numpy as np
import csv
import os

x_single = [100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000]
y_single = [0.162222,0.373598,0.667444,1.00523,1.2732,1.51329,1.76971,2.03468,2.25925,2.50103]

x_complete = [100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000]
y_complete = [0.43321,1.26368,2.40165,3.84472,5.08712,6.65768,8.44603,10.6289,12.4928,14.631]

# filename = '/home/nachiket/rla_cl_exact/chart-data/output_stats-Qgram-K14.csv'
# with open(filename, newline='') as csvfile:
#     filereader = csv.reader(csvfile, delimiter=',', quotechar='|')
#     next(filereader)
#     for row in filereader:
#         if row[0] == 'Single Linkage':
#             x_single.append(row[5])
#             y_single.append(row[1])
#         else:
#             x_complete.append(row[5])
#             y_complete.append(row[1])

plt.plot(x_single, y_single, label = 'Without subtracting 97')

plt.plot(x_complete, y_complete,label = 'With subtracting 97')

plt.title ("ASCII- Manipulation Linking Time")
plt.legend(labels = ('Without subtracting 97','With subtracting 97'),loc='upper left')

plt.show()

