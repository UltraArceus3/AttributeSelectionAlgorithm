import pseudopeople as psp
import time
from record_linkage_feature_selection_apriori import generate_rules,save_to_csv  
from mlxtend.preprocessing import TransactionEncoder
from mlxtend.frequent_patterns import apriori, association_rules,fpgrowth
from anomaly_detection import interest_measures
import pandas as pd

config_dec = {
    'decennial_census': {
        'column_noise': {
            'first_name': {
                'make_typos': {
                    'cell_probability': 1.0,
                },
            },
            'street_name' : {
                'make_typos' : {
                    'cell_probability': 0.2,
                },
            },
        },
    },
}



def perform_association_analysis():
    
    start = time.Time()
    val = 1
    #df = df[df['match'] != val].drop_duplicates(subset=['Attr_1', 'Attr_2', 'Attr_3', 'Attr_4'])
    input_dataset = census.values.tolist()
    tr = TransactionEncoder()
    tr_arr = tr.fit(input_dataset).transform(input_dataset)
    census = pd.DataFrame(tr_arr, columns=tr.columns_)
    #print(df)
    frequent_itemsets = fpgrowth(census, min_support = 0.001, use_colnames = True)
    #print(frequent_itemsets)
    rules = generate_rules(frequent_itemsets)
    N = 15
    final_rules = dict(sorted(rules.items(), key=lambda item: item[1][0], reverse=True)[:N])

    print(final_rules)
    # print("Saving rules to CSV...")
    #save_to_csv(rules)
    # print("Saved to CSV!")
    end = time.time()

    print("Time Taken", (end - start))


def main():

    sample_ds_1 = r"C:\Users\Nachiket Deo\Documents\Rule-Based-ML\data\pse_dec.1.1"
    sample_ds_2 = r"C:\Users\Nachiket Deo\Documents\Rule-Based-ML\data\pse_dec.1.2"

    census = psp.generate_decennial_census(source=r"C:\Users\Nachiket Deo\Documents\pseudopeople_input_data_ri_1_0_0\pseudopeople_input_data_ri_1.0.0")    
    census = census.astype(str)
    
    census.apply(lambda x: x.lower())

    print(census.columns)

    total_len_half = len(census) // 2

    sample_df_1 = census.iloc[:total_len_half,:]
    sample_df_2 = census.iloc[total_len_half:,:]
    print(len(sample_df_1))

    ## First Copy
    sample_df_1.to_csv(sample_ds_1,sep='\t',header=False,index=False)
    sample_df_2.to_csv(sample_ds_2,sep='\t',header=False,index=False)
    
    ## Second Copy

    sample_df_1.to_csv(sample_ds_1,sep='\t',mode='a',header=False,index=False)
    sample_df_2.to_csv(sample_ds_2,sep='\t',mode='a',header=False,index=False)

    ## Third Copy

    sample_df_1.to_csv(sample_ds_1,sep='\t',mode='a',header=False,index=False)
    sample_df_2.to_csv(sample_ds_2,sep='\t',mode='a',header=False,index=False)


    #################################################################################################################################################################################
    census_2 = psp.generate_decennial_census(source=r"C:\Users\Nachiket Deo\Documents\pseudopeople_input_data_ri_1_0_0\pseudopeople_input_data_ri_1.0.0",config=config_dec)    
    census_2 = census_2.astype(str)
    
    print(census_2.columns)

    census_2.apply(lambda x: x.lower())

    total_len_half = len(census_2) // 2

    sample_df_1 = census_2.iloc[:total_len_half,:]
    sample_df_2 = census_2.iloc[total_len_half:,:]
    print(len(sample_df_1))

    ## Error induced copy
    sample_df_1.to_csv(sample_ds_1,sep='\t',mode='a',header=False,index=False)
    sample_df_2.to_csv(sample_ds_2,sep='\t',mode='a',header=False,index=False)

    
if __name__ == "__main__":
    main()