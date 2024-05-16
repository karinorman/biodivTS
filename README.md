# Functional Change through Time
This repository contains the code necessary to reproduce analyses and figures associated with the manuscript "No general trend for multiple dimensions of functional diversity in bird and mammal communities".

It is archived on Zenodo at: https://doi.org/10.5281/zenodo.5514335

## Data

All scripts for accessing and processing raw data are found in the `data-raw` folder with scripts of the same name as the data source. For data sources with more than one script they, execution order follows the number prefix. 

Raw data and interim and final data products are archived on Zendo at https://doi.org/10.5281/zenodo.6499442. 

## Workflow

All scripts for the analysis are stored in `analysis/` and should be executed following the number prefix. Brief descriptions of the contents of each script as well as input and output data products are given below. More detailed comments are given in the body of the scripts.

File Name | Description | Input | Output
--------- | ----------- | ----- | ------
01_filter_by_traits.Rmd | Merge trait and timeseries data. | `biotime_data.rda`, `elton_mamm.rda`, `elton_bird.rda`, `amphibio.rda`, `meta.rda` | `bt_traitfiltered.rda`, `trait_ref.rda`
02_rarefy_timeseries| Calculate rarefied metrics. | `bt_traitfiltered.Rmd` | samples and metrics in `rarefied_metrics/` and `rarefied_samples/`, (file for each sample)
03_rarefy_null_models.Rmd | Get null model samples for each rarefied sample, calculated metrics. | `rarefied_samples/`, `bt_traitfiltered.Rmd` | `null_table.Rmd`
04_collate_rarefied_resamps_median.Rmd | Combine rarefied metrics with null model stats to get final dataframe of metrics. | `rarefied_metrics/`, `null_table.rda`, | `rarefied_metrics.rda`
05_model_metrics.Rmd | Format data for modeling and run all models. | `meta.rda`, `rarefied_metrics.rda` | `model_data.rda`, `meta_clean.rda`, `metric_model_table.rda`, `indv_mod_table.rda`, `model_output.rda`, `sensitivity_output.rda`
06_study_change.Rmd | Examine and group types of study-level change. | `model_output.rda`, `model_data.rda`| `change_group_ses_desc.rda`
