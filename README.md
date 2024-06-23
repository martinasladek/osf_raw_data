# RMT: Data and Models in Psychology

This repository contains data and supplementary material for Chapter 5 of the thesis Role of Robust Statistical Methods in the Credibility Movement of Psychological Science. 

The R folder contains data and all the analysis scripts. The folder is organised as follows: 

- `data` folder: 
  - `helper_data/raw_data_database_osf_public.csv` contains the list of included and excluded studies, reasons for exclusions, links to raw data and dois/links to papers. Each paper has been assigned an ID. This ID is used to for the matching dataset files and analysis files specified below. 
  - `processed_data/mod_df.RDS` contains the processed data. Use this file along with `r_docs/data_processing` to get additional information about the metrics of interest
  - `raw_data` contains all raw data from included studies. 
- `r_docs` folder: 
  - Data processing (`data_processing.qmd`) and analyses (`data_analysis.qmd`) ran on processed data to obtain typical estimates for metrics (skewness, kurtosis, heteroscedasticity, etc)
  - Individual analysis files for each paper labelled with appropriate ID. 
- `objects` folder: 
  - `lm` folder has the exports of all linear models from individual papers 
  - `models` contains the results and MCMC draws from all Bayesian models fitted to obtain typical estimates. 
- `scripts` folder contains helper functions used for re-analyses of linear models and for extracting information about residual metrics. 


