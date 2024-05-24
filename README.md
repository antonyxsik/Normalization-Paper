# Normalization-Paper

This is a supplementary repository that houses experiments and figures for the paper: 

*Normalizing Basis Functions: Approximate Stationary Models for Large Spatial Data*

**Authors**: Antony Sikorski, Daniel McKenzie, and Douglas Nychka

A preprint of the paper is currently available at: [arXiv](https://arxiv.org/abs/2405.13821)

## Cite this as: 

If you wish to cite this paper, please do so with the following BibTeX: 

```{bibtex}
@misc{sikorski2024normalizing,
      title={Normalizing Basis Functions: Approximate Stationary Models for Large Spatial Data}, 
      author={Antony Sikorski and Daniel McKenzie and Douglas Nychka},
      year={2024},
      eprint={2405.13821},
      archivePrefix={arXiv},
      primaryClass={stat.CO}
}
```

## Abstract

> In geostatistics, traditional spatial models often rely on the Gaussian Process (GP) to fit stationary covariances to data. It is well known that this approach becomes computationally infeasible when dealing with large data volumes, necessitating the use of approximate methods. A powerful class of methods approximate the GP as a sum of basis functions with random coefficients. Although this technique offers computational efficiency, it does not inherently guarantee a stationary covariance. To mitigate this issue, the basis functions can be “normalized” to maintain a constant marginal variance, avoiding unwanted artifacts and edge effects. This allows for the fitting of nearly stationary models to large, potentially non-stationary datasets, providing a rigorous base to extend to more complex problems. Unfortunately, the process of normalizing these basis functions is computationally demanding. To address this, we introduce two fast and accurate algorithms to the normalization step, allowing for efficient prediction on fine grids. The practical value of these algorithms is showcased in the context of a full spatial workflow on a large dataset, where significant computational speedups are achieved. While implementation and testing is done specifically within the LatticeKrig framework, these algorithms can be adapted to other basis function methods operating on regular grids. 

---

## Contents

### `LatticeKrigRPackage/`

This folder contains the development (research) version of the **`LatticeKrig`** package. 
- All scripts in this repo should prompt you to download this in order to reproduce the experiments. The first code block will look like this:
- The user simply needs to change the working directory when prompted. A local install of the **`LatticeKrig`** package will then be made on your machine. The line should look like this:
  ```R
  install.packages("LatticeKrigRPackage/LatticeKrig", repos = NULL, type="source")
  ```
- If one wishes to uninstall the development version and re-download the default version of **`LatticeKrig`** from CRAN, they can do so with these two lines (should also be available to uncomment in the `figure_generation.Rmd` script): 
  ```R
  #removes development version
  remove.packages("LatticeKrig")
  #downloads default from CRAN
  install.packages("LatticeKrig")  
  ```

### `timing_error/` 

This folder contains the majority of the timing and error experiments, along with the bulk of the figure generation for the paper. 
- `timing_error_script.R`:  Performs all of the experiments. This takes multiple separate nights to run, and thus I would avoid running it. The end product of this script is all of the `.rda` files that appear in the `timing_error/dataframes/` folder.
- `figure_generation.Rmd` :  Creates all visualizations for Sections 1-4 of the paper, along with a bit of basic hyperparameter exploration at the end. This file should not take long to run. This script creates almost everything in the `timing_error/figures/` folder. The process visualization for FFT normalization was created using an external image editing application.

### `big_data_analysis/`:
This folder contains the experiments for Section 5, where the entire workflow of a large, simulated climate dataset is timed, and accuracy is investigated.

- `simulated_fit_predict.Rmd` : This creates the full data, train and test split, fits the model, and generates prediction surfaces. The data is saved in the `big_data_analysis/dataframes/` folder. The resulting fits and predictions, along with the times, would typically be saved in the `big_data_analysis/results/` folder. The scripts are still set up to save them there, but due to the fits being such large files, we instead choose to store them in [this Google Drive folder](https://drive.google.com/drive/folders/1zPAspblrd8kHy2XLfZfRspxae5X4__HI?usp=sharing). The `big_data_analysis/results/` folder on Github should also contain a text file with the link to this folder.
- `simulated_figures.Rmd` : This script loads in all data, fits, and predictions. It assess accuracy using various metrics, and creates tables summarizing the results. It also creates tons of figures of data, predicted surfaces, and artifacts. The two figures in Section 5 of the paper are created here, and are saved in the `big_data_analysis/figures/` folder. 

### `latticekrig_functions/`:
This folder highlights which functions in the LatticeKrig package were modified. It is not necessary to do anything with these, as they will be installed on your machine through `LatticeKrigRPackage/`. Regardless, if one is interested, they may look here. 

- The `LKrigNormalizeBasisFFTInterpolate_scale.R` function contains the FFT based algorithm for specific grid sizes. 
- The `LKrigNormalizeBasisFFTInterpolate.R` function is a slight modification that allows for much more flexibility with coarse and fine grid sizes, but results in content bleeding. A shift at the end of the algorithms attempts to account for content bleeding. 
- The `LKRectangleFastNormalization.R` function calls `findNormNEW.f` to perform the exact normalization using the Kronecker based algorithm. 
- Also contains: `LKrig.sim.conditional.foreach.R` function for conditional simulation in parallel (optional), the `fillGrid.R` function (convenient grid checking), and the `interp.surface.FFT.R` function, which can help someone build their own fft upsampling function for any application. The fft interp function currently uses a scale factor, but can be easily adjusted to take any coarse and fine grid sizes by following the method in `LKrigNormalizeBasisFFTInterpolate.R`. 


