# bayesian_hdx
Bayesian Analysis of DHDX-MS Data. Calculate the magnitude and significance of a perturbation in HDX at quasi-residue resolution.

## Usage
The code runs as a python script referencing libraries contained in the repository.  See the two examples for sample scripts.

## Input
HDX-MS data can read directly from the output .csv file from HDXWorkbench.

Data can also be read from a .csv file with the following format:
```
# peptide_seq, start_res, end_res, time, D_inc 
AAMNST, 1, 6, 10, 3.212346
AAMNST, 1, 6, 30, 8.5279405
AAMNST, 1, 6, 90, 20.9023379
```

See the `examples/simulated_system/data` for examples.  Support for other HDX data formats will be added as requested.

### Running Time
Analysis for a system of 50-75 peptides takes between 3-6 hours, depending on the amount of overlap and speed of the processor.

## Output
Data output is delivered in into the `output` (or user-defined) directory. 

* DHDX Plots - Plots the DHDX magnitude and significance between the apo (or reference) state and perturbed state.

* Fragment Chi Plots - Plots for each state showing the fragment overlap and colored by the fit of that fragment data to the model.

* Fragment fit-to-data - Plots of time vs. %D incorporation showing the fit of the model to the data.
