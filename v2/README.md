# bayesian_hdx Version 2
Calcualting the distribution of potential protection factors at each residue from MS-HDX deuterium incorporation data.

By explicitly modeling residue-resolved protection factors and accounting for variable noise in the experimental data, we can generate sets of protection factors that fit to MS-HDX peptide deuterium incorporation observations. The result is a set of 

### Requirements
* python 2.7 or 3.X. 
* numpy
* scikit-learn
* matplotlib
* pylab

## Usage

#### As python scripts
The code runs as a python script that can be modified by the user.  Look at the modeling_new.py script in `examples/simulated_system/` for further explanation.  Simply run the script as: 
```
./path/to/v2/setup.sh python modeling_new.py
```
The resulting output folder specified in the modeling script will contain, for each state, an output file and all of the experimental data for that state, consolidated into a single JSON file (.hxd). 

## Input Data
HDX-MS data can read directly from the output .csv file from HDXWorkbench and MassSpec Studio. Support for other HDX data formats will be added as requested.

Data can also be read from a .csv file with the following format:
```
# peptide_seq, start_res, time, D_inc, Score [optional] 
AAMNST, 1, 10, 3.212346
AAMNST, 1, 30, 8.5279405
AAMNST, 1, 90, 20.9023379
```
Currently, the "Score" column is not used, but may be utilized to weight data in future updates.

See the `examples/simulated_system/data` for examples of this file format. 

See `examples/HDXWorkbench_example/modeling_new.py` for a detailed workflow using an MS Workbench file.

### Running Time
Analysis for a system of 50-75 peptides of 10,000 steps takes between one to four hours, depending on the amount of overlap and speed of the processor.
```
run_type = "benchmark"
```
To get a rough estimate for how long your run will take, edit the modeling script such that run_type="benchmark". The script will calculate the approximate time to run 1000 frames.

## Output
Data output is delivered in into a user-defined output directory. A single output file is created containing links to the original data, individual models/scores and sigmas and a conversion table for converting model values to protection factors. 

## Analysis

Multiple independent runs using the same protshould be performed. These runs can be analyzed to extract the best scoring solutions, test for sampling convergence and return clusters of good scoring models. In addition, results from different runs can be compared to test for similarity, or identify the mean and significance of DHDX (similar to v1.0). Examples of these functions are described in:

`examples/CytC/analysis_convergence.py`

A number of plots are used to view HDX data. 
* HDX Plots - Plots the ensemble of HDX protection factors from the best scoring models at each residue for one or more states.

* Fragment Chi Plots - Plots for each state showing the fragment overlap and colored by the fit of that fragment data to the model. Chi values over 20 are suspect.

* Fragment fit-to-data - Plots of time vs. %D incorporation showing the fit of the model to the data. Poor fits can indicate insufficient sampling and/or inconsistent data. 

These are described in 
`examples/CytC/analysis_plotting.py`
