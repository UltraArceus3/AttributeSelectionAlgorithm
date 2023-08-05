# AttributeSelectionAlgorithm
1. Information about the dataset / sampling rate / stages of pipeline etc are to be entered in **config.yaml**

2. Run the following code to install the dependencies associated with Python pcakages
  `python3 -m pip install -r requirements.txt`

## Infomation about the pipeline
1. We generate the random sample from the dataset. (Sampling rate is set in config.yaml)
2. We use RLA_CL code written in C++ to generate the processed data. (RLA_CL is code for record linkage algorithm that was published in 2016)
3. The processed data file contains lot of extraneous data so we prune in this step.
4. Attribute selection algorithm using association rule mining is run.

