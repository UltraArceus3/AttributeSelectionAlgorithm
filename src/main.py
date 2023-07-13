import yaml
import polars as pl

with open('config.yaml', 'r') as f:
    config = yaml.safe_load(f)

col_names_psedu = ['simulant_id', 'first_name', 'middle_initial', 'last_name', 'age',     
       'date_of_birth', 'street_number', 'street_name', 'unit_number', 'city',
       'state', 'zipcode', 'relation_to_reference_person', 'sex',
       'race_ethnicity']

col_names = config['header']
"""
    'simulant_id':pl.Utf8, 
    'first_name':pl.Utf8, 
    'middle_initial':pl.Utf8, 
    'last_name':pl.Utf8, 
    'age':pl.Utf8,     
    'date_of_birth':pl.Utf8, 
    'street_number':pl.Utf8, 
    'street_name':pl.Utf8, 
    'unit_number':pl.Utf8, 
    'city':pl.Utf8,
    'state':pl.Utf8, 
    'zipcode':pl.Utf8, 
    'relation_to_reference_person':pl.Utf8, 
    'sex':pl.Utf8,
    'race_ethnicity':pl.Utf8
"""

df_1 = pl.read_csv(r"../data/pse_dec.1.1", separator = '\t', has_header = False, new_columns = col_names,
                   dtypes={x: pl.Utf8 for x in col_names})
df_2 = pl.read_csv(r"../data/pse_dec.1.2", separator = '\t', has_header = False, new_columns = col_names,
                   dtypes={x: pl.Utf8 for x in col_names})

#attribute_selection_nosnn -> 
#C++ for generating processed dataset -> 
#feature_selection_data_prune -> 
#apriori.py



        ''' Formatting Data '''
        df_merged = pl.concat([df1, df2])
        #df_sorted = df_merged.sort(by=['Last_Name','First_Name','DOB','DOD'],descending=False)
        df_sorted = df_merged.sort(by=['first_name', 'middle_initial', 'last_name', 'age',     
       'date_of_birth', 'street_number', 'street_name', 'unit_number', 'city',
       'state', 'zipcode', 'relation_to_reference_person', 'sex',
       'race_ethnicity'], descending=False)

if __name__ == "__main__":
