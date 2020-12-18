[![Build Status](https://github.com/salilab/bayesian_hdx/workflows/build/badge.svg?branch=master)](https://github.com/salilab/bayesian_hdx/actions?query=workflow%3Abuild)
[![codecov](https://codecov.io/gh/salilab/bayesian_hdx/branch/master/graph/badge.svg)](https://codecov.io/gh/salilab/bayesian_hdx)

# bayesian_hdx
Bayesian Analysis of HDX-MS Data. 

----------------------
## Version 2.0 is live
The major change in V2.0 is the explicit calculation of protection factors for each residue, rather than the difference between two states. Estimations of the magnitude and significance of difference, the focus of v1.0, is essentially unchanged, however, I do see some small, but insignificant numerical differences between the two methods. 

The documentation for the v2.0 workflow is found in the v2 directory.  Eventually, v1 will be deprecated.  

### Notable changes/additions
* 2-10 fold increase in speed
* Calculation of intrinsic unfolded exchange rates based on Bai, Englander 1993.
* Consolidated output file to allow restarting of modeling runs
* Bug fixes
----------------------

## Version 1.0 doc

Calculate the magnitude and significance of a perturbation in HDX at quasi-residue resolution.

### Requirements
* python 2.7 or 3.X. 
* numpy
* scikit-learn
* scipy
* matplotlib

## Usage

#### As python script
The code runs as a python script that can be modified by the user.  Look at the modeling.py script in `examples/simulated_system/` for further explanation.  Simply run the script as: 
```
python modeling.py
```
#### From the command line with a workbech file
Alternatively, if you have an HDX Workbench file, you can input it directly to the `workbench_executable.py` script with the following format:
```
python workbench_executable.py "path/to/workbench/file.csv" "output_directory"
```
Additional command line arguments can be added. Run `python workbench_executable.py -h` to see all options.  Ensure that the ./pyext/src/ folder is in your PYTHONPATH or add it at invocation with the flag `--path "path/to/code/pyext/src"`.

## Input Data
HDX-MS data can read directly from the output .csv file from HDXWorkbench.

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
Data output is delivered in into the `output` (or user-defined) directory. 

* HDX Plots - Plots the ensemble of calculated HDX protection factors at each residue for one or more states. Residue number is on the X-axis (top) 

![Violin plot](img/violins.png)

* Fragment Chi Plots - Plots for each state showing the fragment overlap and colored by the fit of that fragment data to the model.

* Fragment fit-to-data - Plots of time vs. %D incorporation showing the fit of the model to the data.
