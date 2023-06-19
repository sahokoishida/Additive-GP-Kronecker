# Additive-GP-Kronecker

This is a repository for the [paper](https://arxiv.org/abs/2305.07073) : Efficient additive Gaussian process models for large scale data and application to NO2 concentration in London

### Repository structure
The structure of the repository is something like ...

```
├── Code
│   ├── Stan
│   │   ├── GP_helpers.stan *Collection of user-specified functions for estimation*
│   │   ├── GP_helpers_pred.stan *Collection of user-specified functions for prediction*
│   │   ├── GP_1d.stan *An example file for a simple 1d GP regression*
|   |   ├── GP_3d_kron_alltwoway_est.stan *For estimating model parameters in 3d grid structure data, but with 3d grid and all two-way interaction model*
|   |   └── GP_3d_kron_alltwoway_pred.stan *Same as above, but for prediction*
│   ├── R
|   |   ├── NO2_London_example.R *data analysis example with London NO2 data (with stan files, MCMC), to be updated*
|   |   ├── GP_helpers.R *Helper functions, including user specified kernels, to come*
|   |   ├── GP_2d_kron.R *Model parameter estimation by optimising marginal likelihood, to come*
|   |   └──
│   └──Python
│       ├── NO2_London_data_scraping.ipynb *data scraping, to come*
|       └── NO2_London_example.ipynb *data analysis example with London to come*
└── Data
    ├── London_air_station.csv *information on monitoring stations*v
    ├── NO2_STadjust.csv *NO2 concentration from monitoring stations*
    └── NO2_imputation_STadjust.csv *NO2 concentration from monitoring stations after imputation*

```
The main file for data analysis is `NO2_London_example.R`. 

### Resource
The raw data is collected from [London Air](https://www.londonair.org.uk/)
