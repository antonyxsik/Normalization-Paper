# Normalizing Basis Functions

This is a supplementary repository that houses experiments and figures for the paper: 

*Normalizing Basis Functions: Approximate Stationary Models for Large Spatial Data*

**Authors**: Antony Sikorski, Daniel McKenzie, and Douglas Nychka

The article can be found [here](https://onlinelibrary.wiley.com/doi/10.1002/sta4.70015), and the revised preprint is available on [arXiv](https://arxiv.org/abs/2405.13821).

## Cite this as: 

If you use this code in your work or simply need to cite this paper, please do so with the following BibTeX: 

```{bibtex}
@article{sikorski2024normalizing,
  title={Normalizing Basis Functions: Approximate Stationary Models for Large Spatial Data},
  author={Sikorski, Antony and McKenzie, Daniel and Nychka, Douglas},
  journal={Stat},
  volume={13},
  number={4},
  pages={e70015},
  year={2024},
  publisher={Wiley Online Library}
}
```

## Abstract

> In geostatistics, traditional spatial models often rely on the Gaussian Process (GP) to fit stationary covariances to data. It is well known that this approach becomes computationally infeasible when dealing with large data volumes, necessitating the use of approximate methods. A powerful class of methods approximate the GP as a sum of basis functions with random coefficients. Although this technique offers computational efficiency, it does not inherently guarantee a stationary covariance. To mitigate this issue, the basis functions can be “normalized” to maintain a constant marginal variance, avoiding unwanted artifacts and edge effects. This allows for the fitting of nearly stationary models to large, potentially non-stationary datasets, providing a rigorous base to extend to more complex problems. Unfortunately, the process of normalizing these basis functions is computationally demanding. To address this, we introduce two fast and accurate algorithms to the normalization step, allowing for efficient prediction on fine grids. The practical value of these algorithms is showcased in the context of a full spatial workflow on a large dataset, where significant computational speedups are achieved. While implementation and testing is done specifically within the LatticeKrig framework, these algorithms can be adapted to other basis function methods operating on regular grids. 

---

## LatticeKrig

The **`LatticeKrig`** R package is required to reproduce the scripts, visualizations, and experiments in this repository. Previously, this repository housed a development version of the package. The package has since been published to [CRAN](https://cran.r-project.org/web/packages/LatticeKrig/) and is on version 9.3.0. In order to install **`LatticeKrig`**, simply type the following in R. 
```R
install.packages("LatticeKrig") 
```
This line is at the top of all scripts and RMarkdown notebooks in this repository. 

## Contents

### `timing_error/` 

This folder contains the majority of the timing and error experiments, along with the bulk of the figure generation for the paper. 
- `timing_error_script.R`:  Performs all of the experiments. This takes multiple separate nights to run, and thus I would avoid running it. The end product of this script is all of the `.rda` files that appear in the `timing_error/dataframes/` folder.
- `figure_generation.Rmd` :  Creates all visualizations for Sections 1-4 of the paper, along with a bit of basic hyperparameter exploration at the end. This file should not take long to run. This script creates almost everything in the `timing_error/figures/` folder. The process visualization for FFT normalization was created using an external image editing application.

### `big_data_analysis/`:
This folder contains the experiments for Section 5, where the entire workflow of a large, simulated climate dataset is timed, and accuracy is investigated.

- `simulated_fit_predict.Rmd` : This creates the full data, train and test split, fits the model, and generates prediction surfaces. The data is saved in the `big_data_analysis/dataframes/` folder. The resulting fits and predictions, along with the times, would typically be saved in the `big_data_analysis/results/` folder. The scripts are still set up to save them there, but due to the fits being such large files, we instead choose to store them in [this Google Drive folder](https://drive.google.com/drive/folders/1zPAspblrd8kHy2XLfZfRspxae5X4__HI?usp=sharing). The `big_data_analysis/results/` folder on Github should also contain a text file with the link to this folder.
- `simulated_figures.Rmd` : This script loads in all data, fits, and predictions. It assess accuracy using various metrics, and creates tables summarizing the results. It also creates tons of figures of data, predicted surfaces, and artifacts. The two figures in Section 5 of the paper are created here, and are saved in the `big_data_analysis/figures/` folder. 

### `latticekrig_functions/`:
This folder contains the functions that were modified or added to the **`LatticeKrig `** v9.3.0 release, along with some helper functions used in the experiments. It is not necessary to do anything with these, as all necessary functions will be installed on your machine when the package is downloaded. Regardless, if one is interested, they may look here. A few noteworthy functions are: 

- The `LKrigNormalizeBasisFFTInterpolate_scale.R` function contains the FFT based algorithm for specific grid sizes. 
- The `LKrigNormalizeBasisFFTInterpolate.R` function is a slight modification that allows for much more flexibility with coarse and fine grid sizes, but results in content bleeding. A shift at the end of the algorithms attempts to account for content bleeding. 
- The `LKRectangleFastNormalization.R` function calls `findNormNEW.f` to perform the exact normalization using the Kronecker based algorithm. 
- Also contains: `LKrig.sim.conditional.foreach.R` function for conditional simulation in parallel (optional), the `fillGrid.R` function (convenient grid checking), and the `interp.surface.FFT.R` function, which can help someone build their own fft upsampling function for any application. The fft interp function currently uses a scale factor, but can be easily adjusted to take any coarse and fine grid sizes by following the method in `LKrigNormalizeBasisFFTInterpolate.R`. 


