# Normalization-Paper
This is a supplementary repository that houses code for experiments and figures for the paper: 

Normalizing Basis Functions: Approximate Stationary Models for Large Spatial Data
Authors: Antony Sikorski, Daniel McKenzie, and Douglas Nychka

A preprint of the paper is currently available at: LINKHERE

---

Description of Files: 

Folder "LatticeKrigRPackage": 
This folder contains the experimental, research version of the LatticeKrig package. All scripts in this repo should urge you to download this in order to reproduce the experiments. All one needs to do is change the working directory when prompted. A local install of the "LatticeKrig" package will then be made on your machine.

Folder "latticekrig_functions":
- This folder highlights which functions in the LatticeKrig package were modified. It is not necessary to do anything with these, as they will be installed on your machine through "LatticeKrigRPackage". Regardless, if one is interested, they may look here. 
- Also contains the "LKrig.sim.conditional.foreach.R" function for conditional simulation in parallel (optional) and the fillGrid.R  function (convenient grid checking).

Folder "timing_error": 
This folder contains the majority of the timing and error experiments, along with the bulk of the figure generation for the paper. 
- timing_error_script.R: Performs all of the experiments. This takes multiple separate nights to run, and thus I would avoid running it. The end product of this script is all of the ".rda" files that appear in the "timing_error/dataframes" folder.
- figure_generation.Rmd: Creates all visualizations for Sections 1-4 of the paper, along with a bit of basic hyperparameter exploration at the end. This file should not take long to run. This script creates almost everything in the "timing_error/figures" folder. The process visualization for FFT normalization was created using an external image editing application.


