# Bayesian analyses of the temporal relationship between climatic and demographic events: a case-study from the Chulmun period, Korea: source code, data, and scripts

This repository contains all data and scripts required to fully reproduce all analyses presented in the following paper: 

Kim,H., Lee,G., Crema, E.R. Bayesian analyses of the temporal relationship between climatic and demographic events: a case-study from the Chulmun period, Korea.

The repository is organised as follows:
  -  `./runscript` ... contains R scripts for executing core analyses
  - `./R_image_files` ... contains R image files with results
  - `./data` ... contains all datasets required for the analysis
  - `./figures_and_results` ... contains all figures and results summaries and the R scripts required for their creation.
  
## Data Sets and Data Preparation

### Data Sets
The core dataset consists of:

  - a collection of radiocarbon dates gathered from various excavation reports and journal articles (`./data/Neolithic_C14_dates.csv`);
  - sediment data and associated radiocarbon dates from the SSDP-102 core (`./data/Kim_etal_2004/*`) originally published by [Kim et al 2004](https://doi.org/10.1016/j.quascirev.2004.08.010) And retrieved from [LiPDverse Global Holocene Repository](http://lipdverse.org/globalHolocene/current_version/SSDP_102.Kim.2004.html).
  - arboreal to total pollen ratio and associated radiocarbon dates from the Pomaeho sediment core (`./data/Constantine_etal_2020/*`) originally published by [Constantine et al 2020](https://doi.org/10.1017/qua.2018.132).
  - arboreal to total pollen ratio and associated radiocarbon dates from the GY-1 sediment core (`./data/Park_etal_2019/*`) originally published by [Park et al 2019](https://doi.org/10.1038/s41598-019-47264-8).
  - ESRI Shapefiles of the Korean Peninsula from [GADM](https://gadm.org/index.html) (`/data/shp/*`).

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

## SPD Permutation Tests
The file `./runscript/statistical_tests.R` contains R scripts for executing mark permutation tests for the Coastal vs Inland regions and for Millet vs non-Millet dates. Results are stored in the R image file `/R_image_files/spd_test_results.RData`

## Age Depth Models
The file `./runscript/age_depth_model.R` contains R scripts for fitting Compound Poisson-Gamma age-depth models for the SSDP-102, Pomaeho, and the GY-1 sedimentary data using the
[Bchron](https://cran.r-project.org/web/packages/Bchron/index.html) R package. Results are stored in the R image file `./R_image_files/kim2004_agedepthmodel.RData`.

## Bayesian Analysis
The files `./runscript/mcmc_coastal.R`, `./runscript/mcmc_inland.R`, and `./runscript/mcmc_all.R` contains the scripts for fitting a bounded double exponential model via the [nimble](https://r-nimble.org/) and [nimbleCarbon](https://cran.r-project.org/web/packages/nimbleCarbon/index.html) R packages, whilst `mcmc_diagnostic_ppc.R` contains the R code for MCMC diagnostics and posterior predictive checks. R image files containing the MCMC outputs (`mcmc_samples_inland.RData`,`mcmc_samples_coatsal.RData`, and `mcmc_samples_all.RData`) and the posterior predictive checks (`mcmcdiagnostic_postpredcheck.RData`) are loacted in the `R_image_files` directory.


# File Structure

```
├── data
│   ├── Constantine_etal_2020
│   │   ├── pomaeho_constantine_etal_2020_aptp.csv
│   │   └── pomaeho_constantine_etal_2020_c14.csv
│   ├── Kim_etal_2004
│   │   ├── SSDP102_Kim_etlal_2004_c14.csv
│   │   └── SSDP102_Kim_etlal_2004_temp.csv
│   ├── koreanC14.RData
│   ├── Neolithic_C14_dates_cleaned.csv
│   ├── Neolithic_C14_dates.csv
│   ├── Park_etal_2019
│   │   ├── GY_Park_etal_2019_aptp.csv
│   │   └── GY_Park_etal_2019_c14.csv
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
│   ├── age_depth_hpd90.csv
│   ├── figure1.pdf
│   ├── figure2.pdf
│   ├── figure3.pdf
│   ├── figure4.pdf
│   ├── figureS10.pdf
│   ├── figureS11.pdf
│   ├── figureS12.pdf
│   ├── figureS13.pdf
│   ├── figureS1.pdf
│   ├── figureS2.pdf
│   ├── figureS3.pdf
│   ├── figureS4.pdf
│   ├── figureS5.pdf
│   ├── figureS6.pdf
│   ├── figureS7.pdf
│   ├── figureS8.pdf
│   ├── figureS9.pdf
│   ├── main_text_figures.R
│   ├── mcmc_summary.csv
│   ├── results_summary.R
│   └── supplementary_figures.R
├── NeolithicKoreaDemography.Rproj
├── README.md
├── R_image_files
│   ├── agedepthmodels.RData
│   ├── koreanC14.RData
│   ├── mcmcdiagnostic_postpredcheck.RData
│   ├── mcmc_samples_all.RData
│   ├── mcmc_samples_coastal.RData
│   ├── mcmc_samples_inland.RData
│   ├── spd_test_results.RData
│   └── test_results.RData
└── runscript
    ├── age_depth_model.R
    ├── data_clean.R
    ├── data_prep.R
    ├── mcmc_all.R
    ├── mcmc_coastal.R
    ├── mcmc_diagnostic_ppc.R
    ├── mcmc_inland.R
    └── statistical_tests.R

```


# R SessionInfo

```
attached base packages:
[1] stats     graphics  grDevices utils     methods   base     

other attached packages:
 [1] latex2exp_0.4.0    Bchron_4.7.3       coda_0.19-4        truncnorm_1.0-8   
 [5] nimbleCarbon_0.1.2 nimble_0.10.1      raster_3.3-13      rgeos_0.5-3       
 [9] rgdal_1.5-18       maptools_1.0-1     rworldmap_1.3-6    dbscan_1.1-5      
[13] sp_1.4-5           rcarbon_1.4.1      here_0.1           dplyr_1.0.2       

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.5            lattice_0.20-41       deldir_0.2-9         
 [4] snow_0.4-3            rprojroot_2.0.2       foreach_1.5.1        
 [7] ggforce_0.3.2         plyr_1.8.6            R6_2.5.0             
[10] ggridges_0.5.2        spam_2.5-1            ggplot2_3.3.2        
[13] tensor_1.5            pillar_1.4.7          rlang_0.4.8          
[16] rstudioapi_0.13       rpart_4.1-15          Matrix_1.2-18        
[19] goftest_1.2-2         startup_0.14.1        splines_4.0.3        
[22] stringr_1.4.0         foreign_0.8-80        igraph_1.2.6         
[25] polyclip_1.10-0       munsell_0.5.0         spatstat.data_1.7-0  
[28] compiler_4.0.3        xfun_0.20             pkgconfig_2.0.3      
[31] mgcv_1.8-31           doSNOW_1.0.19         tidyselect_1.1.0     
[34] tibble_3.0.4          codetools_0.2-16      crayon_1.3.4         
[37] MASS_7.3-51.6         grid_4.0.3            nlme_3.1-148         
[40] gtable_0.3.0          lifecycle_0.2.0       magrittr_2.0.1       
[43] scales_1.1.1          stringi_1.5.3         farver_2.0.3         
[46] spatstat_1.64-1       ellipsis_0.3.1        generics_0.1.0       
[49] vctrs_0.3.5           spatstat.utils_1.20-2 iterators_1.0.13     
[52] tools_4.0.3           glue_1.4.2            tweenr_1.0.1         
[55] purrr_0.3.4           maps_3.3.0            fields_10.3          
[58] abind_1.4-5           parallel_4.0.3        yaml_2.2.1           
[61] colorspace_2.0-0      dotCall64_1.0-0       knitr_1.30     
```

# Funding
This project was supported by the Korean Studies Promotion Service/Academy of Korean Studies (Laboratory Program for Korean Studies, AKS-2015-Lab-2250001)

# Licence
CC-BY 3.0
