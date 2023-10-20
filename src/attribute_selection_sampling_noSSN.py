import polars as pl
import random
import Levenshtein
import yaml
from typing import Union
from edge_generation import edge_generation_c

import numpy as np

with open('../config.yaml', 'r') as f:
    config = yaml.safe_load(f)


BLOCKING_ATTRIBUTE = config["BLOCKING_ATTRIBUTE"]
# 3 - Last Name for Pseduopeople

SAMPLE_RATE = config["SAMPLE_RATE"]
THRESHOLD = config["THRESHOLD"]

dict_blocks = {}

output_records_sample = []

##
##
##

# DF -  [[1],[1],[1],[2],[2],[2]]
# DF_D - [[1],[2]]

# dict{0} --> 0
# dict{1} --> 3

"""

FUNCTION de_duplication
INPUT: dataset (Pandas DataFrame)
OUTPUT: de-duplicate Dataset (that contains unique records) and a dictionary with information about record 
in de-duplicate dataset and it's location in the original dataset.
information about the location is needed later to find exact matches.

    for e.g
    Original-  [Rany,Kamel], [Rany,Kamel]  [Nachiket,Deo] , [Nachiket Deo]
    De-duplicate-  [Rany,Kamel] [Nachiket,Deo]

"""


def de_duplication(dataset: pl.DataFrame):
    dataset_pand = dataset.to_pandas()  # Converts to pandas dataframe
    dataset_pand_lst = dataset_pand.values.tolist()  # Converts to list of list
    index_original_df = {}
    # adds first element to de duplicated data
    de_duplicate = [dataset_pand_lst[0]]
    index_original_df[0] = 0

    # ['Rany','Kamel'], ['Rany','Kaml']
    #

    for i in range(1, len(dataset_pand_lst)):
        if dataset_pand_lst[i] != dataset_pand_lst[i-1]:
            de_duplicate.append(dataset_pand_lst[i])
            len_d = len(de_duplicate)
            index_original_df[len_d-1] = i
        else:
            continue

    for i in range(len(de_duplicate)):
        de_duplicate[i] = [str(x).lower() for x in de_duplicate[i]]

    # pl_de_duplicate = pd.DataFrame(de_duplicate)
    # print("De-duplicate",de_duplicate[1])
    # print("Original",dataset_pand_lst[0:5])
    # print(index_original_df[1])
    # print(pl_de_duplicate)
    return de_duplicate, index_original_df


"""

FUNCTION: do_blocking
INPUT: de-duplicate dataset (list)
OUTPUT: DOES NOT RETURN ANYTHING but fills the information about records in the respective blocks.
Uses K-mer blocking.

REFERS: generate_k_mer()

"""


def do_blocking(de_duplicated_set: list):

    for i in range(len(de_duplicated_set)):
        tmp = de_duplicated_set[i]
        k_mer_list = generate_k_mer(tmp[BLOCKING_ATTRIBUTE], 3)

        for asc in k_mer_list:
            if asc.lower() in dict_blocks:
                dict_blocks[asc.lower()].append(i)
            else:
                dict_blocks[asc.lower()] = [i]


"""

FUNCTION generate_k_mer
INPUT: A single String, k

"""

def generate_k_mer(str_d, k:int):
    if len(str_d) <= k:
        return [str_d]

    return [str_d[i:i+k] for i in range(0, len(str_d)-(k-1))]


def generatePairs(df_sorted: pl.DataFrame, index_duplicate: dict, key: int):

    output_records_sample.append(df_sorted[index_duplicate[key]])
    if (key + 1 < len(index_duplicate)) and (index_duplicate[key] != index_duplicate[key+1] - 1):
        output_records_sample.append(df_sorted[index_duplicate[key] + 1])



"""

FUNCTION random_data_generation
INPUT:  df_sorted (dataframe containing original data in sorted order)
        de_duplicate(The list containing unique elements)
        index_duplicate(Dictionary containing information about location of records in de_duplicate inside original dataset)
        sample_rate(Sampling rate)

OUTPUT: DOES NOT RETURN ANYTHING but generates sample of the original dataset based on the sampling rate.

REFERS: generate_k_mer(),generatePairs()

"""

def random_data_generation(df_sorted: pl.DataFrame, de_duplicate: list, index_duplicate: dict, sample_rate: float = SAMPLE_RATE) -> None:
    """

        1.Generate Random records from de_duplicate
        2.Total number of records to be generated are total_number of records (de_duplicate) * SAMPLE_RATE
        3.For every randomly generated record, find it's k-mers and the blocks associated to it.
        4.Go through every other record in the block from previous step and find the levenshtein distance betweehn them
        5.If the distance between records is less than or equal to threshold then add them to output.
        6.Find the Exact associated records (Records with no errors) and push them to output (Call generatePairs function for this)

    """

    samp = random.sample(range(0, len(de_duplicate)),
                         int(len(de_duplicate) * sample_rate))

    # de_duplicate_np = np.array(de_duplicate,dtype=np.str)
    bytes_array = np.array([[s.encode() for s in row] for row in de_duplicate], dtype=bytes)

    samp_np = np.array(samp,dtype=int)
    dict_blocks_map = {key.encode(): np.array(value,dtype=int) for key, value in dict_blocks.items()}
    #print("Working till here!")

    record_pre_out = edge_generation_c(bytes_array,samp_np,BLOCKING_ATTRIBUTE,dict_blocks_map,THRESHOLD)
    
    for i in range(len(record_pre_out)):
        generatePairs(df_sorted, index_duplicate, record_pre_out[i])

    # for i in range(len(samp)):
    #     #print("Generating k-mers...", f"{i}/{len(samp)}", "\t" * 5, end='\r')
    #     kmer_index = set()
    #     record = de_duplicate[samp[i]]
    #     k_mer_list = generate_k_mer(record[BLOCKING_ATTRIBUTE], 3)

    #     for kmer in k_mer_list:
    #         if kmer in dict_blocks:
    #             # print(kmer)
    #             kmer_index = kmer_index.union(dict_blocks[kmer.lower()])

    #     # print("Generating pairs...")

    #     for index in kmer_index:
    #         # generatePairs(df_sorted, index_duplicate, index)
    #         tmp_record = de_duplicate[index]
    #         l_dist = 0

    #         for i in range(1, len(tmp_record)):
    #             a = tmp_record[i]
    #             b = record[i]

    #             l_dist += Levenshtein.distance(a, b)

    #             if l_dist > THRESHOLD:
    #                 l_dist = THRESHOLD + 1
    #                 break
            
    #         # if l_dist <= THRESHOLD:
    #         #     record_pre_out.append(index)

    #         if l_dist <= THRESHOLD:
    #             generatePairs(df_sorted, index_duplicate, index)

    


"""

FUNCTION write_output
INPUT: sample_ds_1(String,file_name for output_1),sample_ds_2(String,file_name for output_2)
OUTPUT:Writes the record of the Sample to output
Two output paths are needed as dataset is split into two files (Necessary condition for running the RLA code)

"""


def write_output(cols, _sep="\t", sample_ds_1="sample_ds_1.txt", sample_ds_2="sample_ds_2.txt"):
    _flip = 0  # 0 for ds_1, 1 for ds_2
    def flip_cond(): return _flip % 2

    with open(sample_ds_1, 'w', encoding="utf-8") as f1, open(sample_ds_2, 'w', encoding="utf-8") as f2:

        for record in output_records_sample:
            def rec(ind): return str(record[ind].item())
            # record_csv = _sep.join([rec("SSN") , rec("Last_Name"), rec("First_Name"), rec("DOB"), rec("DOD")]) + "\n"
            record_csv = _sep.join([rec(x) for x in cols]) + "\n"

            # print(record_csv)
            if flip_cond():
                f1.write(record_csv)
            else:
                f2.write(record_csv)
            _flip += 1


def time_code(df: pl.DataFrame, columns: list, tr_samps: list, out_file_1="./data/pse_sample.1.1", out_file_2="./data/pse_sample.1.2") -> list:
    import time

    times = []

    for tr in tr_samps:
        print(f"Total Rate: {tr}")
        output_records_sample.clear()

        start = time.time()

        de_duplicate, index_duplicate = de_duplication(df)
        do_blocking(de_duplicated_set=de_duplicate)

        ''' Random Data Generation '''
        random_data_generation(df, de_duplicate, index_duplicate, tr)
        end = time.time()
        times.append(end - start)

        print(f"Time taken: {times[-1]} seconds")

        write_output(columns, sample_ds_1=out_file_1, sample_ds_2=out_file_2)

    return times


def attribute_sampling(df: pl.DataFrame, columns: list, output_files: Union[str, list]  = "./out.1"):

    tr_samp = [SAMPLE_RATE]

    if type(output_files) is str:
        out1 = f"{output_files}.1"
        out2 = f"{output_files}.2"
    elif type(output_files) is list:
        out1 = output_files[0]
        out2 = output_files[1]
    else:
        raise TypeError("output_files must be a string or a list of strings")

    times = time_code(df, columns, tr_samp, out_file_1=out1, out_file_2=out2)

    return times


"""

NOT USED
WAS NEEDED TO GENERATE PLOTS FOR RESEARCH PAPER

"""


def _plot(times, tr_samples, fig_name="./fig.png", xticks=None, yticks=None):
    import matplotlib.pyplot as plt
    plt.plot(times, [int(x*100)
             for x in tr_samples], marker='o', linestyle='--')
    plt.xlabel("Run Time (s)")
    plt.ylabel("Total Rate (%)")

    if xticks:
        plt.xticks(xticks)

    if yticks:
        plt.yticks(yticks)

    plt.grid()
    plt.savefig(fig_name)
    plt.show()


if __name__ == "__main__":
    _plot(*attribute_sampling(),
          xticks=[5500, 6000, 6500, 7000, 7500, 8000, 8500], yticks=[2, 3, 4, 5, 6])


# print(df_SSN)
