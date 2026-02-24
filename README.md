This repository contains the R implementation of the methodology presented in:

**A Poisson Factor Mixture Model for the Analysis of Linguistic Competence in Italian University Students’ Writing**

The proposed method introduces a Poisson Factor Mixture Model for the analysis of multivariate count data. The model combines a finite mixture structure to account for unobserved heterogeneity with a latent layer that captures possible dependence among variables.

The included files are the following:

- `pfmm_functions.R`  
  It contains the core functions implementing the Poisson factor mixture model. 
  The main function requires the following inputs:
  
  - `y`: data matrix  
  - `k`: number of groups  
  - `r`: number of latent variables  
  - `x`: covariate matrix  
  
  The script includes the routines for model estimation, likelihood evaluation, latent structure modeling, and clustering based on posterior probabilities.

- `example_data_generation_and_model_run.R`  
  It is an example file illustrating how to generate synthetic count data and how to run the proposed model using the functions provided in `pfmm_functions.R`. The script shows the full procedure, from data generation to model fitting and extraction of results.

The code is written in R and reproduces the methodological framework and empirical analyses described in the paper.