
import Levenshtein
import numpy as np
from cython.parallel cimport prange


def generate_k_mer(str_d, k:int):
    if len(str_d) <= k:
        return [str_d]

    return [str_d[i:i+k] for i in range(0, len(str_d)-(k-1))]


def edge_generation(str[:,:] de_duplicate, long[:] samp,blocking_attribute,dict_blocks,threshold):

    cdef list record_pre_out = []

    for i in range(len(samp)):
        #print("Generating k-mers...", f"{i}/{len(samp)}", "\t" * 5, end='\r')
        kmer_index = set()
        record = de_duplicate[<int>samp[i]]
        k_mer_list = generate_k_mer(record[blocking_attribute], 3)

        for kmer in k_mer_list:
            if kmer in dict_blocks:
                # print(kmer)
                kmer_index = kmer_index.union(dict_blocks[kmer.lower()])

        # print("Generating pairs...")

        for index in kmer_index:
            # generatePairs(df_sorted, index_duplicate, index)
            tmp_record = de_duplicate[index]
            l_dist = 0

            for i in range(1, len(tmp_record)):
                a = tmp_record[i]
                b = record[i]

                l_dist += Levenshtein.distance(a, b)

                if l_dist > threshold:
                    l_dist = threshold + 1
                    break
            
            if l_dist <= threshold:
                record_pre_out.append(index)
    return record_pre_out