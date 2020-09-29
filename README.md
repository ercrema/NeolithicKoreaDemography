# Title title title: source code, data, and scripts

This repository contains all data and scripts required to fully reproduce all analyses presented in the following paper: 

_Paper Title_

The repository is organised in the following main directories:
  - ./runscript ... contains R scripts for executing core analyses
  - ./Rimages ... contains R image files with results
  - ./data ... contains all datasets required for the analysis
  - ./src ... contains additional bespoke R functions for Approximate Bayesian Computation
  - ./figures_and_results ... contains all figures and results summaries as well as the R script required to generate them (`figurelog.R` and `results_summary.R`)
  
## Data Sets and Data Preparation

### Data Sets
The core dataset consists of a collection of radiocarbon dates gathgered from excatation reports and journal article by Habeom Kim (`./data/Neolithic_C14_dates.csv`), marine sediment data and associated radiocarbon dates originally published by [Kim et al 2000](https://doi.org/10.1016/j.quascirev.2004.08.010) (`./data/kim2004/SSDP_102.Kim.2004.csv` and `./data/kim2004/SSDP_102.Kim.2004-chron.csv`) and retrieved from [LiPDverse Global Holocene Repository](http://lipdverse.org/globalHolocene/current_version/SSDP_102.Kim.2004.html), and shapefiles for the Korean Peninsula from the [GADM](https://gadm.org/index.html) (`/data/shp/*`).

### Pre-processing radiocarbon dates
The original set of radiocarbon dates have been processed using the R scipts `./runscript/data_clean.R` and `./runscript/data_prep.R`. The former processes `./data/Neolithic_C14_dates.csv` by:

1. Removing all instances of dates without LabCode or a clear source
2. Removing all dates outside the time window 6200-2800 14C Age
3. Subsetting only to terrestial dates associated with charcoal or wood
4. Removing sites without geographic coordinates
5. Assigning a unique site identifier (`SiteID`)
6. Creating a new field `milletAsso` for all cases where the radiocarbon dates are associated to millets.
7. Keeping a subset of fields for the analysis and storing the output in the file `./data/Neolithic_C14_dates_cleaned.csv`

The `./runscript/data_prep.R` function processes `./data/Neolithic_C14_dates_cleaned.csv` by grouping site into clusters using the DBSCAN algorithm (to control for CRM-based arbitrary definitions of archaeological sites) and subsequently assigning each cluster to either the _coastal_ or _inland_ region based on distance from coast. The function also calibrates all dates using IntCal20 and generates radiocarbon bins for each cluster to control for inter-site variation in sampling intensity. The results are store in the R image file `/R_image_files/koreanC14.RData`. 

## SPD Analyses
The file `./runscript/statistical_tests.R` contains R scripts for executing the following set of analysis via the [rcarbon](https://cran.r-project.org/web/packages/rcarbon/index.html) R package:
 - Composite Kernel Density Estimates with Bootstrapping for Coastal and Inland regions.
 - Mark Permutation test of Coastal vs Inland regions.
 - Theorethical Growth Model (Exponential and Logistic) for the whole study area and per region.

The results are stored in the R image file `/R_image_files/spd_test_results.RData`

## Age Depth Model
The file `./runscript/age_depth_model.R` contains R scripts for an age-depth model of Kim et al 2004 sedimentary data (`./data/kim2004/SSDP_102.Kim.2004.csv` and `./data/kim2004/SSDP_102.Kim.2004-chron.csv`) using the Compound Poisson-Gamma chronology model provided by the [Bchron](https://cran.r-project.org/web/packages/Bchron/index.html) R package. Results are stored in the R image file `./R_image_files/kim2004_agedepthmodel.RData`.

## Approximate Bayesian Computation
The files `./runscript/abc_rej_*.R` includes the R scripts for run the Approximate Bayesian Computation (ABC) analysis on the observed SPD of each region (`abc_rej_coastal.R` and `abc_rej_inland.R`) as well as the entire dataset (`abc_rej_general.R`). The script requires the [rcarbon](https://cran.r-project.org/web/packages/rcarbon/index.html) R package as well as additional bespoke functions contained in the `./src/` directory which includes a code for steamlined faster calibration (`fastCalibrate.R`) and the core simulation model that is fitted to the data (`sim_model.R`). Results of the ABC, as well as posterior predictive checks (executed using the scripts `./runscript/ppcheck_*.R`) are stored in individual images files (`./R_image_files/resABC_laplace_*.RData` and `./R_image_files/predcheck_results_*.RData`).

# File Structure

.
├── data
│   ├── kim2004
│   │   ├── SSDP_102.Kim.2004-chron.csv
│   │   └── SSDP_102.Kim.2004.csv
│   ├── koreanC14.RData
│   ├── Neolithic_C14_dates_cleaned.csv
│   ├── Neolithic_C14_dates.csv
│   └── shp
│       ├── polyline_korea.cpg
│       ├── polyline_korea.dbf
│       ├── polyline_korea.prj
│       ├── polyline_korea.sbn
│       ├── polyline_korea.sbx
│       ├── polyline_korea.shp
│       ├── polyline_korea.shp.xml
│       └── polyline_korea.shx
├── figures_and_results
│   ├── abc_result_hpd90.csv
│   ├── age_depth_hpd90.csv
│   ├── figure_ckde.pdf
│   ├── figure_climate_vs_changepoint.pdf
│   ├── figure_event_comparisons.pdf
│   ├── figure_kim2004_reanalysis.pdf
│   ├── figurelog.R
│   ├── figure_millet_multiplot.pdf
│   ├── figure_millet_spd.pdf
│   ├── figure_modelTests.pdf
│   ├── figure_permtest_millet.pdf
│   ├── figure_permtest.pdf
│   ├── figure_posterior_coastal_vs_inland.pdf
│   ├── figure_posterior_general.pdf
│   ├── figure_ppcheck_coastal_vs_inland.pdf
│   ├── figure_ppcheck_general.pdf
│   ├── figure_site_map.pdf
│   ├── figure_stacked_spd.pdf
│   ├── nhst_result.csv
│   └── results_summary.R
├── README.md
├── R_image_files
│   ├── kim2004_agedepthmodel.RData
│   ├── koreanC14.RData
│   ├── predcheck_results_coastal.RData
│   ├── predcheck_results_general.RData
│   ├── predcheck_results_inland.RData
│   ├── predcheck_results.RData
│   ├── resABC_laplace_coastal.RData
│   ├── resABC_laplace_general.RData
│   ├── resABC_laplace_inland.RData
│   ├── resABC_laplace.RData
│   ├── spd_test_results.RData
│   └── test_results.RData
├── runscript
│   ├── abc_rej_coastal.R
│   ├── abc_rej_general.R
│   ├── abc_rej_inland.R
│   ├── age_depth_model.R
│   ├── data_clean.R
│   ├── data_prep.R
│   ├── ppcheck_coastal.R
│   ├── ppcheck_general.R
│   ├── ppcheck_inland.R
│   └── statistical_tests.R
└── src
    ├── fastCalibrate.R
    └── sim_model.R
```
```


# R Settings

```
attached base packages:
[1] parallel  stats     graphics  grDevices utils    
[6] methods   base     

other attached packages:
 [1] doSNOW_1.0.18     snow_0.4-3       
 [3] doParallel_1.0.15 iterators_1.0.12 
 [5] foreach_1.5.0     Bchron_4.7.1     
 [7] dplyr_1.0.2       raster_3.3-13    
 [9] rgeos_0.5-3       rgdal_1.5-12     
[11] maptools_1.0-1    rworldmap_1.3-6  
[13] dbscan_1.1-5      sp_1.4-2         
[15] rcarbon_1.4.1    

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.5            lattice_0.20-41      
 [3] deldir_0.1-28         digest_0.6.25        
 [5] ggforce_0.3.2         plyr_1.8.6           
 [7] R6_2.4.1              ggridges_0.5.2       
 [9] evaluate_0.14         spam_2.5-1           
[11] ggplot2_3.3.2         tensor_1.5           
[13] pillar_1.4.6          rlang_0.4.7          
[15] rstudioapi_0.11       rpart_4.1-15         
[17] Matrix_1.2-18         goftest_1.2-2        
[19] startup_0.14.1        rmarkdown_2.3        
[21] splines_4.0.2         stringr_1.4.0        
[23] foreign_0.8-80        polyclip_1.10-0      
[25] munsell_0.5.0         spatstat.data_1.4-3  
[27] compiler_4.0.2        xfun_0.17            
[29] pkgconfig_2.0.3       mgcv_1.8-31          
[31] htmltools_0.5.0       tidyselect_1.1.0     
[33] tibble_3.0.3          codetools_0.2-16     
[35] crayon_1.3.4          MASS_7.3-51.6        
[37] grid_4.0.2            nlme_3.1-148         
[39] gtable_0.3.0          lifecycle_0.2.0      
[41] magrittr_1.5          scales_1.1.1         
[43] stringi_1.5.3         farver_2.0.3         
[45] spatstat_1.64-1       ellipsis_0.3.1       
[47] generics_0.0.2        vctrs_0.3.4          
[49] spatstat.utils_1.17-0 RColorBrewer_1.1-2   
[51] tools_4.0.2           glue_1.4.2           
[53] tweenr_1.0.1          purrr_0.3.4          
[55] maps_3.3.0            fields_10.3          
[57] abind_1.4-5           yaml_2.2.1           
[59] colorspace_1.4-1      dotCall64_1.0-0      
[61] knitr_1.29   
```

# Funding
E.R.Crema was funded by a Philip Leverhulme Prizes.

# Licence
CC-BY 3.0
