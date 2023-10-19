import yaml
import polars as pl
from attribute_selection_sampling_noSSN import attribute_sampling
import os
import subprocess
import pandas as pd
from parse_xml import parse_xml
from feature_selection_data_prune import prune_feature_data
from record_linkage_feature_selection_apriori import run_pipeline as run_apriori
import time

with open('../config.yaml', 'r') as f:
    config = yaml.safe_load(f)


# attribute_selection_nosnn ->
# C++ for generating processed dataset ->
# feature_selection_data_prune ->
# apriori.py

RLA_DIR = "../RLA_CL_EXTRACT/"


def sample_generation():
    id = config['id_column']

    col_names = config['header']

    cols_wo_id = [x for x in col_names if x != id]

    df1_file = config["input_files"]["input_file_1"]
    df2_file = config["input_files"]["input_file_2"]

    df1 = pl.read_csv(df1_file, separator='\t', has_header=False, new_columns=col_names,
                      dtypes={x: pl.Utf8 for x in col_names})
    df2 = pl.read_csv(df2_file, separator='\t', has_header=False, new_columns=col_names,
                      dtypes={x: pl.Utf8 for x in col_names})

    ''' Formatting Data '''
    df_merged = pl.concat([df1, df2])
    # df_sorted = df_merged.sort(by=['Last_Name','First_Name','DOB','DOD'],descending=False)
    df_sorted = df_merged.sort(by=cols_wo_id, descending=False)

    attribute_sampling(df_sorted, col_names, output_files=list(
        config["sample_output"].values()))


def call_RLA():
    parse_xml()

    sample_data = pd.read_csv(
        config["input_files"]["input_file_1"], sep='\t', header=None)

    generate_processed_data = subprocess.check_call([os.path.join(
        RLA_DIR, "bin/rlacl")] + [config["rla_xml_file"], str(len(sample_data)), "0", "0", "0", "0"])

    # subprocess.run(["mv", os.path.join(
    #     RLA_DIR, "feature_selection_processed_data_file.csv"), "../data"])
    subprocess.run(["mv","feature_selection_processed_data_file.csv", "../data"])

def feature_prune():
    prune_feature_data(path="../data/feature_selection_processed_data_file.csv",
                       output=config["pruned_output"], keep=config["FEATURE_PRUNE_KEEP"])


def apriori():
    run_apriori(config["pruned_output"], output_file=config["rules_output"])


if __name__ == "__main__":
    # call_RLA(*list(config["sample_output"].values()))
    # exit()

    if config["run"]["sample_generation"]:
        
        print("Running attribute selection sampling...\n")
        start = time.time()
        sample_generation()
        end = time.time()
        print("Time Taken for Sampling:",(end-start))
        
    if config["run"]["processed_data_generation"]:
        print("Running RLA_CL_EXTRACT...\n")
        call_RLA()

    if config["run"]["feature_prune"]:
        print("Running feature pruning...\n")
        feature_prune()

    start = time.time()
    if config["run"]["attribute_selection"]:
        print("Running apriori...\n")
        apriori()
    end = time.time()
    print("Time Taken for Frequent Itemset",(end-start))

    print("Done!")
