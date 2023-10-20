# AttributeSelectionAlgorithm
1. Information about the dataset / sampling rate / stages of pipeline etc are to be entered in **config.yaml**

2. Run the following code to install the dependencies associated with Python packages
  `python3 -m pip install -r requirements.txt`

## Information about the pipeline
1. We generate the random sample from the dataset. (Sampling rate is set in config.yaml)
2. We use RLA_CL code written in C++ to generate the processed data. (RLA_CL is code for record linkage algorithm that was published in 2016)
3. The processed data file contains lot of extraneous data so we prune in this step.
4. Attribute selection algorithm using association rule mining is run.

## Steps to run the program

### Input dataset
1. Place the dataset files (atleast two are required) into the **data** subdirectory.

### Update configuration file

2. Update the config.yaml file by changing the attribute **input_files** with the names of the input dataset files.
3. Update the **header** attribute  with the names of the attributes for your dataset. (Make sure the headers are not present in your dataset)
4. Change the **id_column** attribute to denote the unique identifier in your dataset.
5. Update the **THRESHOLD** with value suitable for your dataset (This would help program to correct link link the records even when there are errors such as typos are present). The value is usually 1 or 2 (Please note this value should be a positive integer)
6. If you are going to create sample then update **TOTAL_RATE** attribute with the percentage of total data should be considered for random sampling.
7. Since we are using k-mer based blocking while generation of sample, **BLOCKING_ATTRIBUTE** is needed to identify the number of the attribute to be used for blocking. (note blocking speeds up the computation of pairs so choosing right attribute would impact speed)
8. Populate the **COMPARISON_ATTRIBUTE** based on number of attributes in your dataset. for e.g if your dataset has 5 attributes and 0th attribute is unique_identifier then `[1 2 3 4]` should be correct value for comparision attribute. (Note we are using all attributes so that all attributes are included in attribute selection algorithm).

### No Sample run
1. If your dataset is small (total records less than 50,000) then you might not want to generate samples. In the case disable **sample_generation** inside **run** attribute by setting it to false.
2. Update **sample_output** attribute with the path to the input dataset (Same value as **input_files**).
3. DO NOT touch this setting if you are using sampling.

### Running the program
1. Run `python3 setup.py build_ext --inplace` inside **src** folder. This compiles the cython code.
1. In order to run the program, run the file `main.py` inside **src** folder.
2. It is advisable to create sample (Stage 1) by setting the value to True in **run** attribute while keeping every other attribute inside **run** to False. This will create sample and store it in appropriate location.
3. Other stages of the pipeline can run togther by changing the sub attributes in **run** to True and setting **sample_generation** to False. This is done so that program doesn't crash.
4. Once all stages are completed, in order to find the attributes please take a look at **rules_out.csv** inside **data** folder.

### Interpret the Output
1. The headers of the output are attributes from the input dataset set as integer value [0,1,2,3] and M,lift,leverage,convergence where M stands what value was there in Y for association rule (X --> Y)
2. The contents may have values like 0,1,2 which stands for the edit distance for the respective attribute.
3. The output is sorted with best rules at the top (Having high value for lift,conviction and leverage)
   

