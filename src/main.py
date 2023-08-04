import yaml
import polars as pl
from attribute_selection_sampling_noSSN import attribute_sampling
import os, subprocess
import pandas as pd
from parse_xml import parse_xml

with open('../config.yaml', 'r') as f:
    config = yaml.safe_load(f)


#attribute_selection_nosnn -> 
#C++ for generating processed dataset -> 
#feature_selection_data_prune -> 
#apriori.py

RLA_DIR = "../RLA_CL_EXTRACT/"


def call_RLA(out_1, out_2):
    parse_xml()


    sample_data = pd.read_csv(out_1, sep = '\t', header = None)

    generate_processed_data = subprocess.check_call([os.path.join(RLA_DIR, "bin/rlacl")] + ["config_test.xml", str(len(sample_data)), "0","0","0","0"])

    subprocess.run(["mv", os.path.join(RLA_DIR, "feature_selection_processed_data_file.csv"), "../data"])


if __name__ == "__main__":
    #call_RLA(*list(config["sample_output"].values()))
    #exit()

    id = config['id_column']

    col_names = config['header']

    cols_wo_id = [x for x in col_names if x != id]

    df1_file = config["input_files"]["input_file_1"]
    df2_file = config["input_files"]["input_file_2"]

    df1 = pl.read_csv(df1_file, separator = '\t', has_header = False, new_columns = col_names,
                   dtypes={x: pl.Utf8 for x in col_names})
    df2 = pl.read_csv(df2_file, separator = '\t', has_header = False, new_columns = col_names,
                   dtypes={x: pl.Utf8 for x in col_names})
    
    ''' Formatting Data '''
    df_merged = pl.concat([df1, df2])
        #df_sorted = df_merged.sort(by=['Last_Name','First_Name','DOB','DOD'],descending=False)
    df_sorted = df_merged.sort(by = cols_wo_id, descending = False)


    attribute_sampling(df_sorted, output_files = list(config["sample_output"].values()))

    call_RLA()