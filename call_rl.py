import subprocess
import os
import pandas as pd

wd = os.getcwd()

os.chdir("RLA_CL_EXTRACT")

git_stash = subprocess.check_call(['git'] + ['stash'])

git_checkout_status = subprocess.check_call(['git'] + ["checkout","feature_selection"])

make_call = subprocess.run(['make'])

sample_data = pd.read_csv(r'/home/nachiket/AttributeSelectionAlgorithm/RLA_CL_EXTRACT/data/pse/pse_sample3.1.1',sep='\t',header=None)

generate_processed_data = subprocess.check_call(['bin/rlacl'] + ["/home/nachiket/AttributeSelectionAlgorithm/config_pse_sample3.xml",str(len(sample_data)),"0","0","0","0"])

subprocess.run(["mv", "feature_selection_processed_data_file.csv", "../data"])


