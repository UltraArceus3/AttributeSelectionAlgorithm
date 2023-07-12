# import modin.pandas as pd
# import modin.config as modin_cfg

# #df = pd.read_csv(r"C:\Users\Nachiket Deo\Documents\Rule-Based-ML\data\ds5.1.1")

# modin_cfg.Engine.put("ray")

import pandas as pd
import random


sample_ds_1 = r"C:\Users\Nachiket Deo\Documents\Rule-Based-ML\data\ds_sample.1.1"
sample_ds_2 = r"C:\Users\Nachiket Deo\Documents\Rule-Based-ML\data\ds_sample.1.2"

colnames=['SSN', 'Last_Name', 'First_Name', 'DOB','DOD']
df_1 = pd.read_csv(r"C:\Users\Nachiket Deo\Documents\Rule-Based-ML\data\ds5.1.1",sep='\t',header=None,names=colnames)
df_SSN = df_1['SSN'].tolist()

df_2 = pd.read_csv(r"C:\Users\Nachiket Deo\Documents\Rule-Based-ML\data\ds5.1.2",sep='\t',header=None,names=colnames)

# [Kamel,Rany,11,123] -- 120
# [Kamel,Rany,11,123] -- 1
# [Kamel,Ran,11,123]  -- 3



total_rate = 15
total_iterations_SSN = int(len(df_SSN) * (total_rate/100)) 

SSN_list_random = random.sample(df_SSN,total_iterations_SSN)

print(len(SSN_list_random))

output_record_sample = []

for key in SSN_list_random:
    tmp = df_1.loc[df_1['SSN'] == key].values.tolist()
    #output_record_sample.extend(iter(tmp))
    tmp.sort()

    l_1 = len(tmp)
    tmp_processed = []
    if l_1 > 1:
        if tmp[0] != tmp[1]:
            tmp_processed.append(tmp[0])
            if l_1 > 2:
                tmp_processed.extend(tmp[i] for i in range(1,3))
        elif tmp[l_1 - 1] != tmp[l_1 - 2]:
            tmp_processed.append(tmp[l_1-1])
            if l_1 > 2:
                tmp_processed.extend(tmp[i] for i in range(0,2))

    
    tmp_2 = df_2.loc[df_2['SSN'] == key].values.tolist()
    l_2 = len(tmp_2)
    tmp_2_processed = []

    if l_2 > 1:
        if tmp_2[0] != tmp_2[1]:
            tmp_2_processed.append(tmp_2[0])
            if l_2 > 2:
                tmp_2_processed.extend(tmp_2[i] for i in range(1,3))
        elif tmp_2[l_2 - 1] != tmp_2[l_2 - 2]:
            tmp_2_processed.append(tmp_2[l_2-1])
            if l_2 > 2:
                tmp_2_processed.extend(tmp_2[i] for i in range(0,2))

    output_record_sample.extend(iter(tmp_processed))

    for lst in tmp_2_processed:
        output_record_sample.append(lst)

print(output_record_sample[0])
print(output_record_sample[1])
print(output_record_sample[2])
print(output_record_sample[3])
print(output_record_sample[4])
print(output_record_sample[5])
print(output_record_sample[6])
print(output_record_sample[7])
print(output_record_sample[8])

df_output_record_sample = pd.DataFrame(output_record_sample,columns=colnames)

df_output_record_sample = df_output_record_sample.sample(frac=1).reset_index(drop=True)

total_len_half = len(df_output_record_sample) // 2


sample_df_1 = df_output_record_sample.iloc[:total_len_half,:]
sample_df_2 = df_output_record_sample.iloc[total_len_half:,:]
print(sample_df_1)

sample_df_1.to_csv(sample_ds_1,sep='\t',header=False,index=False)
sample_df_2.to_csv(sample_ds_2,sep='\t',header=False,index=False)







#print(df_SSN)