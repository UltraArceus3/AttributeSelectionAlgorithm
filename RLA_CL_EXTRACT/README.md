# Record Linkage Complete Linkage implementation

## Installing on Linux

### Required packages

Use the following command to install all required packages for Debian, Ubuntu Linux

```
sudo apt install libxml2-dev libboost-regex-dev libmpich-dev libboost-log-dev
```

might need also

```
sudo apt install libcr-dev mpich mpich-doc openmpi-bin
```

### How to build

```shell
make clean
make
```

### How to run

```shell
$ bin/rlacl <config_path> <n_records>
```

This program takes two arguments as detailed below:

1. The first argument `<config_path>` is the path for the configuration file. The configuration file has all the parameters required for running the linkage process.

2. The second argument `<n_records>` is the total number of records to read from the dataset(s).


#### Example

To run the linkage job detailed by the parameters of the configuration file `in/config.xml` and limit the number of records to 400 000 records, invoke the following command.

```shell
$ bin/rlacl in/config.xml 400000
```

### Setting the dataset 

```
cd rla_cl_exact/in/

```
1.Open the **config.xml file** and navigate to the tag named **'dataset'**
The value component of dataset has a path to the file where dataset is stored. You need to have at least two dataset files for the program to run correctly.

2.Inside **ComparisonGroup**tag of the config.xml file there's a property called <dist_calc_method> which specifies what type of distance measure is to be used. Please find the following value which program would accept
	**1 - Edit Distance
	2 - REVERSAL Distance
	3 - Truncate Distance
	5 - Q-Gram Distance
	6 - Hausdorff Distance **
3.The **comparing_attribute_indices** property mentions the position of the column from the input dataset that is being used by the algorithm. (Incorrect values would lead to un-defined behaviour)

4.**threshold property** is used by the algorithm to identify when should we consider two strings are not equal and should be in different clusters

5.**priority property** like the comparing attribute indices specify the location of the column in input data that would be used for clustering and which needs to be given priority

### Logging results

The results of the linkage job is saved in `results.csv` file by default in the working directory.
However, one can set a different name in the configuration file.
See XML node: `results_logging/filename/value` and the associated comments.
Also, this can be overridden by specifying the name for logging the results via one of the command line switches `-L` or `--log`.

