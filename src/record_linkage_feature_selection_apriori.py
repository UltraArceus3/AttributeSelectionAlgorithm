import pandas as pd
from mlxtend.preprocessing import TransactionEncoder
from mlxtend.frequent_patterns import fpgrowth
import time
from generateChart import plotChart
from decimal import *
import os
from dataclasses import dataclass
#from arulespy.arules import Transactions, apriori, parameters



"""
    
    FUNCTION save_to_csv
    INPUT: Association rules (list), _file (Output file path), _not_have(string)
    OUTPUT: Output file in the data folder containing good attributes.

"""


def save_to_csv(rules: list, _file="../.output/rules.csv", _not_have="-"):
    vals = {}
    m = []
    lifts = []
    levs = []
    convs = []

    for rule in rules:
        rfs = list(rule[0])  # Breaks down rule frozenset
        rule_added = []

        for r in rfs:
            lr = r.split("_")

            if lr[0] not in vals.keys():
                vals[lr[0]] = [_not_have] * len(m)

            vals[lr[0]].append(lr[1])
            rule_added.append(lr[0])

        for k in vals.keys():
            if k not in rule_added:
                vals[k].append(_not_have)

        m.append(list(rule[1])[0])

        lifts.append(rule[2])
        levs.append(rule[3])
        convs.append(rule[4])

    # Sorts Rule Collection (so it is 1, 2, 3, 4)
    vals = dict(sorted(vals.items(), key = lambda x: int(x[0])))
    vals.keys()

    # Create CSV
    with open(_file, "w") as f:
        vkeys = ",".join(map(str, vals.keys()))

        f.write(f"{vkeys},M,lift,leverage,convergence\n")
        for i in range(len(m)):
            rule_str = ""
            for k in vals.keys():
                rule_str += str(vals[k][i]) + ","

            f.write(f"{rule_str}{m[i]},{lifts[i]},{levs[i]},{convs[i]}\n")


"""
    FUNCTION generate_rules
    INPUT: frequent_itemsets (pandas DataFrame)
    OUTPUT: Association Rules with highest value for lift,conviction and leverage

"""


def generate_rules(frequent_itemsets):
    """
        1. Create a dict to store [frequent_itemset --> support_count]
        2. Find Support Count for 'match' == 1
        3. Iterate through the frequent itemsets and find ones that contains 1.
        4. Split the frequent itemset from previous step into format (X --> Y) where Y is 1 and X is [frequent_itemset - Y]
        5. Find the value for lift,conviction and leverage and store it along with assciation rule [X --> Y]
        6. Sort the stored association rules based on lift and conviction in descending order.

    """

    final_rules = []

    dict_freq_item = {}

    for ind in frequent_itemsets.index:
        dict_freq_item[frequent_itemsets['itemsets']
                       [ind]] = frequent_itemsets['support'][ind]

    # print("Frequent_dict",dict_freq_item[frozenset(['1'])])

    supp_y = dict_freq_item[frozenset(['1'])]

    for ind in frequent_itemsets.index:
        lst_itemsets = list(frequent_itemsets['itemsets'][ind])
        if '1' in lst_itemsets and len(lst_itemsets) > 1:
            tmp = [x for x in lst_itemsets if x != '1']
            fz = frozenset(tmp)
            sp_x = dict_freq_item[fz]
            sp_x_y = frequent_itemsets['support'][ind]
            conf = sp_x_y / sp_x
            lift = conf / supp_y
            conviction = (1 - supp_y) / (1 - conf)
            leverage = sp_x_y - (sp_x * supp_y)
            final_rules.append([fz, '1', lift, leverage, conviction])

        # print(type(frequent_itemsets['itemsets'][ind]))

        # print(frequent_itemsets['support'][ind], frequent_itemsets['itemsets'][ind])

    return sorted(final_rules, key=lambda x: (-x[2], -x[4]))


"""
    FUNCTION run_pipeline
    INPUT: src(path to the input file), output (path to output file)
    OUTPUT: Generates the attributes with good lift,leverage and conviction in form of csv file.
    This function contains the attribute selection algorithm
    REFERS: fpgrowth() [from mlxtend library], generate_rules()

"""


def run_pipeline(src, output_file="../.output/rules.csv"):
    """
        1. Read the data from input file and transofrm it.
        2. Pass it through fp-growth function to generate frequent itemsets.
        3. Pass the frequent itemset to generate_rules() function so that association rules (X --> Y) where Y is match are generated
        4. Push the good rules to output.

    """
    df = pd.read_csv(src, header=None)

    # names=['Attr_1', 'Attr_2', 'Attr_3', 'Attr_4','Attr_5','Attr_6','Attr_7','Attr_8','Attr_9','Attr_10','Attr_11','Attr_12','Attr_13','Attr_14','match']

    cols = [f"Attr_{i}" for i in range(1, len(df.columns))] + ['match']

    df.columns = cols

    df = df.astype(str)
    val = 1
    # df = df[df['match'] != val].drop_duplicates(subset=['Attr_1', 'Attr_2', 'Attr_3', 'Attr_4'])
    input_dataset = df.values.tolist()
    tr = TransactionEncoder()
    tr_arr = tr.fit(input_dataset).transform(input_dataset)
    df = pd.DataFrame(tr_arr, columns=tr.columns_)
    # print(df)
    frequent_itemsets = fpgrowth(df, min_support= 0.0001, use_colnames=True)
    print("Frequent itemsets generated")
    # print(frequent_itemsets)
    rules = generate_rules(frequent_itemsets)
    print(rules)
    # print(rules[0:15])

    print("Saving rules to CSV...")
    save_to_csv(rules, _file=output_file)
    print("Saved to CSV!")



def main():

    # lst = [r'../data/feature_selection_processed_data_file5.csv',r'../data/feature_selection_processed_data_file8.csv',r'../data/feature_selection_processed_data_file10.csv',r'../data/feature_selection_processed_data_file12.csv',r'../data/feature_selection_processed_data_file15.csv']
    lst = [r'../data/fs_pruned_data_pse_3.csv']
    time_data = []
    record_data = [3, 4, 5]

    # df = pd.read_csv(lst[0],names=['Attr_1', 'Attr_2', 'Attr_3', 'Attr_4','match'])
    # df = df.astype(str)
    # val = 1
    # df = df[df['match'] != val].drop_duplicates(subset=['Attr_1', 'Attr_2', 'Attr_3', 'Attr_4'])
    # print(df)

    for src in lst:
        start = time.time()
        run_pipeline(src)
        end = time.time()
        time_data.append(round((end-start), 3))
    print("Time-Data", time_data)
    # plotChart(time_data,record_data,'Time(s)','Sample (%)')


if __name__ == "__main__":
    main()
