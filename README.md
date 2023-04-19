# Additive-GP-Kronecker

This is a repository for the paper: Efficient additive Gaussian process models for large scale data and application to NO2 concentration in London

### Repository structure
The structure of the repository are something like

├── Code
│   ├── Stan
│   │   ├── GP-helpers.stan *Collection of user-specified functions*
│   │   ├── GP-1d.stan *An example file for estimating model parameters in a simple 1d GP regression*
|   |   ├── GP-2d-kron-main.stan *For estimating model parameters in 2d grid structure data, main effect model*
|   |   ├── GP-3d-kron-all-twoway.stan *Same as above, but with 3d grid and all two-way interaction model*
|   |   └── GP-3d-kron-all-twoway-pred.stan *Same as above, but for prediction*
│   ├── R
|   |   ├── NO2-London-example.R *data analysis example with London NO2 data (with stan files, MCMC)*
|   |   ├── NO2-London-data-imutation.R  *missing value imputation*
|   |   ├── GP-helpers.R *Helper functions, including user specified kernels, to come*
|   |   ├── GP-2d-kron.R *Model parameter estimation by optimising marginal likelihood, to come*
|   |   └──
│   └──Python
│       ├── NO2-London-data-scraping.ipynb *data scraping*
|       ├──
|       └──
└── Data
    ├── London-NO2-raw.csv *NO2 concentration from monitoring stations*
    ├── London-air-station.csv *information on monitoring stations*
    └── NO2-imputation.csv *NO2 concentration from monitoring stations after imputation*

### Resource
The raw data is collcted from [London Air](https://www.londonair.org.uk/)
