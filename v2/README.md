# bayesian_hdx Version 2
Calcualting the distribution of potential protection factors at each residue from MS-HDX deuterium incorporation data.

### Requirements
* python 2.7 or 3.X. 
* numpy
* scikit-learn
* matplotlib
* pylab

## Usage

#### As python script
The code runs as a python script that can be modified by the user.  Look at the modeling_new.py script in `examples/simulated_system/` for further explanation.  Simply run the script as: 
```
./path/to/v2/setup.sh python modeling_new.py
```

## Input Data
HDX-MS data can read directly from the output .csv file from HDXWorkbench and MassSpec Studio

Data can also be read from a .csv file with the following format:
```
# peptide_seq, start_res, time, D_inc, Score [optional] 
AAMNST, 1, 10, 3.212346
AAMNST, 1, 30, 8.5279405
AAMNST, 1, 90, 20.9023379
```
Currently, the "Score" column is not used, but will be in upcoming versions

See the `examples/simulated_system/data` for examples.  Support for other HDX data formats will be added as requested.

### Running Time
Analysis for a system of 50-75 peptides takes between 3-6 hours, depending on the amount of overlap and speed of the processor.
```
run_type = "benchmark"
```
To get a rough estimate for how long your run will take, edit the modeling script such that run_type="benchmark". The script will calculate the approximate time to run 1000 frames.


## Output
Data output is delivered in into a user-defined output directory, default = `output`


## Analysis
A number of plots are used to view HDX data. 
* HDX Plots - Plots the ensemble of HDX protection factors from the best scoring models at each residue for one or more states.

* Fragment Chi Plots - Plots for each state showing the fragment overlap and colored by the fit of that fragment data to the model. Chi values over 20 are suspect.

* Fragment fit-to-data - Plots of time vs. %D incorporation showing the fit of the model to the data. Poor fits can indicate insufficient sampling and/or inconsistent data. 
