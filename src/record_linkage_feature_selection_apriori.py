import pandas as pd
from mlxtend.preprocessing import TransactionEncoder
from mlxtend.frequent_patterns import fpgrowth
import time
from generateChart import plotChart
from decimal import *
import os
from dataclasses import dataclass

@dataclass
class Rule:
    antecedent: str
    precendent: str
    lift: float

def truncate(n, decimals=0):
    multiplier = 10 ** decimals
    return int(n * multiplier) / multiplier

class interest_measures:
    lift = lambda p_xIy, supp_x, supp_y, **_ : p_xIy / (supp_x * supp_y) # Range [0, inf), independence = 1
    jaccard = lambda p_xIy, supp_x, supp_y, **_ : p_xIy / (supp_x + supp_y - p_xIy) # Range [0, 1]
    leverage = lambda p_xIy, supp_x, supp_y, **_ : p_xIy - (supp_x * supp_y) # Range [-1, 1], independence = 0
    added_value = lambda p_xIy, supp_x, supp_y, **_ : (truncate(p_xIy, 6) / supp_x) - supp_y # Range [-.5, 1]
    conviction = lambda supp_x, conf_xy, **_ : (1 - supp_x) / (1 - conf_xy) # Range [0, inf), independence = 1
    certainty_factor = lambda conf_xy, supp_y, **_ : (conf_xy - supp_y) / (1 - supp_y) # Range [-1, 1], independence = 0

def certainty_factor(conf_xy,supp_y):
    return (conf_xy - supp_y) / (1 - supp_y)

def added_value(p_xIy,supp_x,supp_y):
    return (p_xIy / supp_x) - supp_y


def save_to_csv(rules: list, _file = "../.output/rules.csv", _not_have = "-"):
    vals = {}
    m = []
    lifts = []
    levs = []
    convs = []

    for rule in rules:
        rfs = list(rule[0]) # Breaks down rule frozenset
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


    vals = dict(sorted(vals.items())) # Sorts Rule Collection (so it is 1, 2, 3, 4)
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


def generate_other_measures(frequent_itemsets):

    for ind in frequent_itemsets.index:
        fz_freq = frequent_itemsets['itemsets'][ind]
        fz_freq_ante = []
        sc_y = 0
        y = ''
        if len(fz_freq) > 1 and ( '0' in fz_freq):
            for item in fz_freq:
                if item != '0':
                    fz_freq_ante.append(item)    
            frequent_itemsets_y = frequent_itemsets[frequent_itemsets['itemsets'] == frozenset('0')]        
            sc_y = frequent_itemsets_y['support'].values[0]
            y = '0'
        elif len(fz_freq) > 1 and ( '1' in fz_freq):
            
            for item in fz_freq:
                if item != '1':
                    fz_freq_ante.append(item)
            frequent_itemsets_y = frequent_itemsets[frequent_itemsets['itemsets'] == frozenset('1')]        
            sc_y = frequent_itemsets_y['support'].values[0]
            y = '1'

        if y != '':
            fz_freq_ante = frozenset(fz_freq_ante)
            frequent_itemsets_ante = frequent_itemsets[frequent_itemsets['itemsets'] == fz_freq_ante]
            sc_x = frequent_itemsets_ante['support'].values[0]
            sc_x_y = frequent_itemsets['support'][ind]
            sc_x = float(sc_x)
            sc_x_y = float(sc_x_y)
            confidence =  sc_x_y / sc_x
            #print("conf",confidence)
            cf = certainty_factor(confidence,sc_y)
            fz_freq_ante_f = frozenset([fz_freq_ante,frozenset(y)])
            if fz_freq_ante_f in final_rules:
                final_rules[fz_freq_ante_f].append(cf)

def generate_top_10_rules_frozen(Q_itemsets: tuple, itemsets_normal: tuple, target, k: int, total_records: int, interest_measure = interest_measures.certainty_factor):

    frequent_item_dict = {}
    itemsets = Q_itemsets[0]
    itemsets_n = itemsets_normal[0]
    rules = []
    for i, freq_itemset in enumerate(itemsets_n):
        frequent_item_dict[frozenset(freq_itemset)] = itemsets_normal[1][i]
        #print(freq_itemset,itemsets_normal[1][i])
    
    for i in range(len(itemsets)):

        support_count__target = 0
        support_count_x = 0

        ##
        ## F-q --> f`, t
        ##

        if type(itemsets[i]) == list and itemsets[i] != target:
            temp = [x for x in itemsets[i] if x not in target]
            trg = [x for x in itemsets[i] if x in target]

            if  frozenset(temp) in frequent_item_dict:
                support_count_x = Decimal(frequent_item_dict[frozenset(temp)]) / Decimal(total_records)
                sup_x = frequent_item_dict[frozenset(temp)]
                                
            if trg:
                sup_y = frequent_item_dict[frozenset(trg)]
                support_count__target = Decimal(frequent_item_dict[frozenset(trg)]) / Decimal(total_records)
            
            if support_count_x == 0.0 or support_count__target == 0.0:
                continue
            
            
            p_xIy = Decimal(Q_itemsets[1][i]) / Decimal(total_records)

            # TODO: Confidence is NOT in range (it's sometimes > 1 ).
            # TODO: This could have to do with how P(X Inter Y) is calculated.
            conf_xy = Decimal(p_xIy) / Decimal(support_count_x)
            #print("I,CONF,SUPP_X",Q_itemsets[1][i],conf_xy,sup_x,support_count_x,Q_itemsets[0][i])
            if conf_xy > 1:
                print("I,CONF,SUPP_X",Q_itemsets[1][i],conf_xy,sup_x,support_count_x,Q_itemsets[0][i])
                continue

            supp_xUy = (sup_x + sup_y - (Q_itemsets[1][i]))/float(total_records)
            # lift = supp(x,y) / (supp(x) * supp(y))
            #xlift = ( (Q_itemsets[1][i]) / float(total_records)) / float((support_count_x) * (support_count__target))

            # added value = (supp(x,y) / supp(x)) - supp(y)
            # tm = ((Q_itemsets[1][i]) / total_records)
            # tm = truncate(tm,6)
            # lift = (tm / support_count_x) - support_count__target

            # jaccard = supp(x U y) / (supp(x) + supp(y) - supp(x U y))
            #jacc = (p_xIy) / (support_count_x + support_count__target - p_xIy)
            
            # Leverage = P(xIy) - P(x)P(y)
            #lift = p_xIy - (support_count_x * support_count__target)

            lift = interest_measure(p_xIy = p_xIy, 
                                    supp_x = Decimal(support_count_x), 
                                    supp_y = Decimal(support_count__target),
                                    conf_xy = Decimal(conf_xy),
                                    supp_xUy = supp_xUy)

            #print("CONF\n", conf_xy, "\nCONF\n")

            if lift >= 0:
                rule = Rule(trg, temp, lift)
                #print(Q_itemsets[0][i],tm,support_count_x,temp,trg,support_count__target,frequent_item_dict[frozenset(trg)],lift)
                rules.append(rule)


    ## Sorting the data based on the target and then lift in decreasing order
    rules.sort(key=lambda x: (x.antecedent,x.lift), reverse=True)
    #print(rules)

    ##
    ## Algorithm for Rule selection.
    ## For each group (associated with each target) we take the rules with top 5 lifts.
    ## The sorting in previous line helps in doing the selection.
    ##

    out_rules = []
    cnt = 0
    current_ante = 0
    for rule_ in rules:

        if cnt < 10 and current_ante == rule_.antecedent:
            current_ante = rule_.antecedent
            out_rules.append(rule_)
            cnt += 1
        elif cnt < 10:
            current_ante = rule_.antecedent
            out_rules.append(rule_)
            cnt = 1
        elif cnt == 10:
            if current_ante != rule_.antecedent:
                cnt = 0

    #print(out_rules)

    os.makedirs('../.output', exist_ok = True)

    #print("\nRules:", rules, "\n")

    print("Rules generated. Saving...")
    with open(r'../.output/rules_3.txt', 'w') as fp:
        max_rule = -1
        for item in out_rules:
        # write each item on a new line
            if item.lift >= -5:
                fp.write("%s\n" % item)
                print(item.precendent,item.antecedent,item.lift)
            if item.lift > max_rule:
                max_rule = item.lift
        print("Max-Rule",max_rule)
    print('Done')




def generate_rules(frequent_itemsets):
    
    final_rules = []

    dict_freq_item = {}

    for ind in frequent_itemsets.index:
        dict_freq_item[frequent_itemsets['itemsets'][ind]] = frequent_itemsets['support'][ind]

    print("Frequent_dict",dict_freq_item[frozenset(['1'])])

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
            leverage =  sp_x_y - (sp_x * supp_y)
            final_rules.append([fz,'1',lift,leverage,conviction])


        #print(type(frequent_itemsets['itemsets'][ind]))

        #print(frequent_itemsets['support'][ind], frequent_itemsets['itemsets'][ind])


    return sorted(final_rules, key=lambda x: (-x[2], -x[4]))
    


def run_pipeline(src):

    df = pd.read_csv(src,names=['Attr_1', 'Attr_2', 'Attr_3', 'Attr_4','Attr_5','Attr_6','Attr_7','Attr_8','Attr_9','Attr_10','Attr_11','Attr_12','Attr_13','Attr_14','match'])
    df = df.astype(str)
    val = 1
    #df = df[df['match'] != val].drop_duplicates(subset=['Attr_1', 'Attr_2', 'Attr_3', 'Attr_4'])
    input_dataset = df.values.tolist()
    tr = TransactionEncoder()
    tr_arr = tr.fit(input_dataset).transform(input_dataset)
    df = pd.DataFrame(tr_arr, columns=tr.columns_)
    #print(df)
    frequent_itemsets = fpgrowth(df, min_support = 0.001, use_colnames = True)
    print("Frequent itemsets generated")
    #print(frequent_itemsets)
    rules = generate_rules(frequent_itemsets)
    print(rules[0:15])

    print("Saving rules to CSV...")
    save_to_csv(rules)
    print("Saved to CSV!")

def main():
    
    #lst = [r'../data/feature_selection_processed_data_file5.csv',r'../data/feature_selection_processed_data_file8.csv',r'../data/feature_selection_processed_data_file10.csv',r'../data/feature_selection_processed_data_file12.csv',r'../data/feature_selection_processed_data_file15.csv']
    lst = [r'../data/fs_pruned_data_pse_3.csv']
    time_data = []
    record_data = [3,4,5]
    
    # df = pd.read_csv(lst[0],names=['Attr_1', 'Attr_2', 'Attr_3', 'Attr_4','match'])
    # df = df.astype(str)
    # val = 1
    # df = df[df['match'] != val].drop_duplicates(subset=['Attr_1', 'Attr_2', 'Attr_3', 'Attr_4'])
    # print(df)


    for src in lst:
        start = time.time()
        run_pipeline(src)    
        end = time.time()
        time_data.append(round((end-start),3))
    print("Time-Data",time_data)
    #plotChart(time_data,record_data,'Time(s)','Sample (%)')

if __name__ == "__main__":
    main()


# [[frozenset({'7_0', '1_1', '5_0'}), '1', 132.7594472117149, 0.006138230034880164, inf], 
#  [frozenset({'7_0', '1_1', '3_0', '5_0'}), '1', 132.7594472117149, 0.006138230034880164, inf], 
# [frozenset({'7_0', '4_0', '1_1', '5_0'}), '1', 132.7594472117149, 0.006138230034880164, inf], 
# [frozenset({'2_0', '7_0', '1_1', '5_0'}), '1', 132.7594472117149, 0.006138230034880164, inf], 
# [frozenset({'7_0', '1_1', '14_0', '5_0'}), '1', 132.7594472117149, 0.006138230034880164, inf], 
# [frozenset({'7_0', '1_1', '13_0', '5_0'}), '1', 132.7594472117149, 0.006138230034880164, inf], 
# [frozenset({'7_0', '1_1', '6_0', '5_0'}), '1', 132.7594472117149, 0.006138230034880164, inf], 
# [frozenset({'11_0', '1_1', '7_0', '5_0'}), '1', 132.7594472117149, 0.006138230034880164, inf], 
# [frozenset({'7_0', '1_1', '9_0', '5_0'}), '1', 132.7594472117149, 0.006138230034880164, inf], 
# [frozenset({'7_0', '1_1', '5_0', '12_0'}), '1', 132.7594472117149, 0.006138230034880164, inf], 
# [frozenset({'7_0', '1_1', '8_0', '5_0'}), '1', 132.7594472117149, 0.006138230034880164, inf], 
# [frozenset({'7_0', '1_1', '5_0', '10_0'}), '1', 132.7594472117149, 0.006138230034880164, inf], 
# [frozenset({'7_0', '4_0', '1_1', '3_0'}), '1', 132.7594472117149, 0.006138230034880164, inf], 
# [frozenset({'2_0', '7_0', '1_1', '3_0'}), '1', 132.7594472117149, 0.006138230034880164, inf], 
# [frozenset({'7_0', '1_1', '3_0', '13_0'}), '1', 132.7594472117149, 0.006138230034880164, inf]]