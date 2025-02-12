# PNAS-WBGT-fertility-Nigeria
Replication code for "Humid-heat shifts conceptions and increases pregnancy termination risks in Nigeria" published in PNAS 2025.

The materials in this repository includes public data and the scripts used to produce the figures and tables appearing in the main text and supplementary information of the paper.

If you find errors in the code or have questions or suggestions, please contact Nina Brooks at nrbrooks@bu.edu.

## Organization of repository
To replicate the analysis your directory should be structured as followed:

data/raw: stores raw DHS data (see note below)

data/clean: stores clean data files for analysis

scripts/: scripts for cleaning, merging data, and analyzing data.

figures/: stores pdf figures output by analysis script.

tables/: stores latex tables produced by analysis script.

## Data
### DHS
We used Nigeria DHS data extracted from [IPUMS DHS](https://www.idhsdata.org/idhs/). DHS data is publicly available, with a registered user account. To replicate this analysis, you must first register for a DHS account and if you wish to use IPUMS DHS to extract the data, also register for an IPUMS DHS account (both free). Both IPUMS and the DHS have extensive guides and tutorials on working with the data. The [IPUMS Climate Change and Health Research Hub](https://tech.popdata.org/dhs-research-hub/) also has resources on how to import IPUMS  data into R using the `ipumsr` package. Once you have accounts, this analysis utilizes the 2008 and 2018 surveys with the fertility calendar module and the complete births history from the 1990, 2003, 2008, 2018 waves. In addition to demographic, fertility, and conceptrceptive use variables, you must also obtain access to the GPS data (which requires specific permission). 

### Heat
We estimate daily maximum shaded wet bulb globe temperature (WBGTmax) for Nigeria DHS clusters from 1983 - 2016 using CHIRTS-daily. The raw CHIRTSdaily heat record is publicly available from the University of California, Santa Barbara Climate Hazards Center at https://data.chc.ucsb.edu/products/CHIRTSdaily/. The approach to produce the WBGTmax from CHIRTSdaily is described in C Tuholske, et al., Global urban population exposure to extreme heat. *Proc. Natl. Acad. Sci.* 118, e2024792118 (2021) and the output is available at https://data.chc.ucsb.edu/people/cascade/UHE-daily/wbgtmax/. This repository contains code to extract the WBGTmax data for the Nigeria DHS clusters.

## Scripts
There are 8 scripts to reproduce the data cleaning and analysis for this paper. Note that some are python scripts and some are Rmd files.

`setup.R`: this script loads necessary R packages and defines custom functions that will be used in the subsequent scripts. 

`01_extractWBGT.py`: this code extracts the daily WBGTmax record in 5km buffer zones around each Nigeria DHS cluster.

`02_annual_avg_wbgtmax.ipynb`: this code produces the annual average WBGTmax across Nigeria used to create Figure 1.

`03_prepare_WBGTdata.Rmd`: creates monthly and annual WBGTmax (from daily record) for each DHS cluster.

`04_import_clean_DHS.Rmd`: imports IPUMS DHS data (from `data/raw`) and cleans variables for analysis. Creates woman-level data, woman-month (from calendar), and the unbalanced panel of births and saves them to `data/clean'.

`05_merge_heat_DHS.Rmd`: merges WBGTmax with all dhs datasets. specifically, merges heat and woman-month calendar data to create survival structured data for discrete-time analysis, merges heat with woman-month calendar to create annual pregnancy-level data, merges heat with birth panel to create annual births datasets.

`06_makeFig1_Nigeria_AvgHeat83-2016.ipynb`: this script creates the map of annual average WBGT and change across Nigeria from 1983-2016.

`07_run_analysis_heat_fertility.Rmd`: runs all analysis, including supplementary analysis, and outputs figures and tables from main paper and SI.
