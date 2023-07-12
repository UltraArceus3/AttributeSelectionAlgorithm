from math import ceil

def prune_feature_data(path: str, output: str, keep: float = 0.20):
    with open(path, "r") as f:
        data = f.readlines()

        features_prune = {}
        features_keep = {}

        for d in data:
            if int(d[-2]) == 0: # Gets last character, if mismatch
                if d in features_prune:
                    features_prune[d] += 1
                else:
                    features_prune[d] = 1
            else:
                if d in features_keep:
                    features_keep[d] += 1
                else:
                    features_keep[d] = 1
                
    
    with open(output, "w") as f:
        for d, count in features_keep.items():
            for _ in range(count):
                f.write(d)
        
        for d, count in features_prune.items():
            prune_perc = ceil(count * keep)
            #print(prune_perc, count)
            for _ in range(prune_perc):
                f.write(d)



if __name__ == "__main__":
    prune_feature_data(path = "../data/feature_selection_processed_data_file_pse_3.csv", output = "../data/fs_pruned_data_pse_3.csv", keep = .2)
    prune_feature_data(path = "../data/feature_selection_processed_data_file_pse_4.csv", output = "../data/fs_pruned_data_pse_4.csv", keep = .2)
    prune_feature_data(path = "../data/feature_selection_processed_data_file_pse_5.csv", output = "../data/fs_pruned_data_pse_5.csv", keep = .2)
