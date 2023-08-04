import polars as pl
import random
import Levenshtein




BLOCKING_ATTRIBUTE = 3
## 3 - Last Name for Pseduopeople

TOTAL_RATE = 8/100
THRESHOLD = 1

dict_blocks = {}

output_records_sample = []

##
##
##

# DF -  [[1],[1],[1],[2],[2],[2]]
# DF_D - [[1],[2]]

# dict{0} --> 0
# dict{1} --> 3


def de_duplication(dataset: pl.DataFrame):
    dataset_pand = dataset.to_pandas() ## Converts to pandas dataframe
    dataset_pand_lst = dataset_pand.values.tolist() ## Converts to list of list
    index_original_df = {}
    de_duplicate = [dataset_pand_lst[0]] ## adds first element to de duplicated data
    index_original_df[0] = 0

    # ['Rany','Kamel'], ['Rany','Kaml']
    #

    for i in range(1,len(dataset_pand_lst)):
        if dataset_pand_lst[i] != dataset_pand_lst[i-1]:
            de_duplicate.append(dataset_pand_lst[i])
            len_d = len(de_duplicate)
            index_original_df[len_d-1] = i
        else:
            continue
    
    for i in range(len(de_duplicate)):
        de_duplicate[i] = [str(x).lower() for x in de_duplicate[i]]


    #pl_de_duplicate = pd.DataFrame(de_duplicate)
    # print("De-duplicate",de_duplicate[1])
    # print("Original",dataset_pand_lst[0:5])
    # print(index_original_df[1])
    # print(pl_de_duplicate)
    return de_duplicate,index_original_df

def do_blocking(de_duplicated_set: list):

    for i in range(len(de_duplicated_set)):
        tmp = de_duplicated_set[i]
        k_mer_list = generate_k_mer(tmp[BLOCKING_ATTRIBUTE],3)
        
        for asc in k_mer_list:
           if asc.lower() in dict_blocks:
                dict_blocks[asc.lower()].append(i)
           else:
                dict_blocks[asc.lower()] = [i]

def generate_k_mer(str_d, k:int) -> list:
    if len(str_d) <= k:
        return [str_d]

    return [str_d[i:i+k] for i in range(0,len(str_d)-(k-1))]



def generatePairs(df_sorted: pl.DataFrame, index_duplicate:dict, key:int):
    
    output_records_sample.append(df_sorted[index_duplicate[key]])
    if (key + 1 < len(index_duplicate)) and (index_duplicate[key] != index_duplicate[key+1] - 1):
        output_records_sample.append(df_sorted[index_duplicate[key] + 1])
    
    
##  Original [Rany,Kamel], [Rany,Kamel]  [Nachiket,Deo] , [Nachiket Deo]   
##  De-duplicate [Rany,Kamel] [Nachiket,Deo]
    # key = 0 [Rany Kamel]
    # Key + 1 = 2 --> [Nachiket Deo]

##
## TO BE IMPLEMENTED BY RANY
##


def random_data_generation(df_sorted: pl.DataFrame, de_duplicate:list, index_duplicate: dict, total_rate: float = TOTAL_RATE) -> None:
    
    ## Generate Random records from de_duplicate
    ## Total number of records to be generated are total_number of records (de_duplicate) * TOTAL_RATE
    ## Find the Exact associated records and push them to output (Call generatePairs function for this)

    ## Random record --> take it's index
    ## Generate K-mer of blocking attribute
    ## Query in dict_blocks and get indexes
    ## Pass the indexes to generatePairs()
    ## [Rany,Kamel] (1) ... [Rane,Kamel] (10)
    ## Kam,Ame,Mel
    # dict_blocks[Kam] = [1,10..]
    # dict_blocks[Ame] = [1,10..]
    # dict_blocks[Mel] = [1,10..]

    # get random sample of indexes
    samp = random.sample(range(0, len(de_duplicate)), int(len(de_duplicate) * total_rate))



    for i in range(len(samp)):
        print("Generating k-mers...", f"{i}/{len(samp)}", "\t" * 5, end='\r')
        kmer_index = set()
        record = de_duplicate[samp[i]]
        k_mer_list = generate_k_mer(record[BLOCKING_ATTRIBUTE], 3)

        for kmer in k_mer_list:
            if kmer in dict_blocks:
                #print(kmer)
                kmer_index = kmer_index.union(dict_blocks[kmer.lower()])

        #print(kmer_index)

        # print("Generating pairs...")
        
        for index in kmer_index:
            #generatePairs(df_sorted, index_duplicate, index)
            tmp_record = de_duplicate[index]
            l_dist = 0
            
            for i in range(1,len(tmp_record)):
                a = tmp_record[i]
                b = record[i]

                l_dist += Levenshtein.distance(a,b)
                
                if l_dist > THRESHOLD:
                    l_dist = THRESHOLD + 1
                    break

            if l_dist <= THRESHOLD:
                generatePairs(df_sorted, index_duplicate, index)
            

def write_output(_sep = "\t", sample_ds_1 = "sample_ds_1.txt", sample_ds_2 = "sample_ds_2.txt", headers = []):
    _flip = 0 # 0 for ds_1, 1 for ds_2
    flip_cond = lambda : _flip % 2


    header = _sep.join(headers) + "\n"

    with open(sample_ds_1, 'w',encoding = "utf-8") as f1, open(sample_ds_2, 'w',encoding = "utf-8") as f2:
        f1.write(header)
        f2.write(header)

        for record in output_records_sample:
            rec = lambda ind : str(record[ind].item())
            #record_csv = _sep.join([rec("SSN") , rec("Last_Name"), rec("First_Name"), rec("DOB"), rec("DOD")]) + "\n"
            record_csv = _sep.join([rec("simulant_id") , rec("first_name"), rec("middle_initial"), rec("last_name"), rec("age"), rec("date_of_birth"), rec("street_number"), rec("street_name"), rec("unit_number"), rec("city"), rec("state"), rec("zipcode"), rec("relation_to_reference_person"), rec("sex"), rec("race_ethnicity")]) + "\n"

            #print(record_csv)
            if flip_cond():
                f1.write(record_csv)
            else:
                f2.write(record_csv)
            _flip += 1


# AALTO     ┆ ROBERT     ┆ 6222003 ┆ 5091963
# aaltop    ┆ robert     ┆ 6222003 ┆ 5091963

def time_code(df: pl.DataFrame, tr_samps: list, out_file_1 = "./data/pse_sample.1.1", out_file_2 = "./data/pse_sample.1.2") -> list:
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

        write_output(sample_ds_1 = out_file_1, sample_ds_2 = out_file_2)

    return times

def attribute_sampling(df: pl.DataFrame, output_files: str | list = "./out.1"):

    tr_samp = [TOTAL_RATE]

    if type(output_files) is str:
        out1 = f"{output_files}.1"
        out2 = f"{output_files}.2"
    elif type(output_files) is list:
        out1 = output_files[0]
        out2 = output_files[1]
    else:
        raise TypeError("output_files must be a string or a list of strings")

    times = time_code(df, tr_samp, out_file_1 = out1, out_file_2 = out2)

    #print(times)
    if type(times) is list:
        times_n = [l[0] for l in times]
    else:
        times_n = times

    return times_n, tr_samp
    


def _plot(times, tr_samples, fig_name = "./fig.png", xticks = None, yticks = None):
    import matplotlib.pyplot as plt
    plt.plot(times, [int(x*100) for x in tr_samples], marker='o', linestyle='--')
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
    _plot(*attribute_sampling(), xticks = [5500,6000,6500,7000,7500,8000,8500], yticks = [2,3,4,5,6])











#print(df_SSN)