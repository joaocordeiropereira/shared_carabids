Code author: João M. Cordeiro Pereira
Affiliation: Chair of Wildlife Ecology & Management, University of Freiburg, Germany

Publication title: Specialist carabid beetles in mixed montane forests show positive links to roe deer and to biodiversity-oriented forestry
Authors: João M. Cordeiro Pereira, Sebastian Schwegmann, Clàudia Massó Estaje, Martin Denter, Grzegorz Mikusinski, Ilse Storch

Last updated: 1. Nov 2023

This README file explains which data files are available, which are necessary to run the analysis R script, and provides metadata about all variables contained in those data files

------ Files available --------------------------------

- R script (carabidae_script.r)
- two files read by R script: all plot-level data (alldata_carabids_df.csv) and carabid traits and classifications (carabid_traits.csv)
- one additional file: trap-level carabid counts (raw_trap_counts.csv)

------ Metadata all_data_carabids_df.csv: -------------

- plot_id: Unique plot ID (CFB followed by three digits)
- X: UTM easting coordinates of plot (CRS: EPSG:31467, DHDN 3-degree Gauss-Kruger zone 3)
- Y: UTM northing coordinates of plot (CRS: EPSG:31467, DHDN 3-degree Gauss-Kruger zone 3)
- n_traps: number of useable traps per plot (from which activity-densities were calculated)
- date_set: date of trap placement (format yyyy-mm-dd)
- duration: number of days traps were active
- trapdays: n_traps*duration, as a measure of sampling effort
- carabidae_n: total carabid activity-density (pooled across traps)
- abax_ovalis to licinus_hoffmanseggii: number of individuals captured for each species (pooled across traps)
- avg_alt: average plot altitude above sea level (m)
- DBHMean: plot-level mean of tree Diameter at Breast Height (DBH), for all trees > 7 cm DBH (in mm)
- lying_dw_volume: volume of coarse woody debris (cubic meters)
- pc_broadleaf: fraction of plot basal area belonging to broadleaved tree species
- roe_deer: daily rate of roe deer detection events
- canopycover: fraction of plot area containing vegetation above 5 m height
- for_cover_100ha: percentage of 100-ha (1 square km) circular buffer around plot that is covered by forest
- for_cover_2500ha: percentage of 2500-ha (25 square km) circular buffer around plot that is covered by forest
- sd_slope: standard deviation of slope values (in degrees) across plot
- northness: northness index (-1 to 1), calculated from plot aspect, -1 being south-facing slopes and 1 north-facing slopes

------ Metadata carabid_traits.csv: ------------------

- species: species name, according to taxonomy in [Ref 1]
- forest_specialist: "yes" or "no", according to [Ref 2]
- montane_specialist: "yes" or "no", according to [Ref 2]
- flight: brachypterous or flight-able, according to [Ref 1]
- body_length: average body length (in mm), according to [Ref 3]
- redlisted_bw: Red list status of species in the state of Baden-Württemberg, as of 2005 [Ref 4]
- redlisted_de: Red list status of species in the Federal Republic of Germany, as of 2016 [Ref 5]
- redlisted: "R" when species is in any category other than * (ungefährdet) in any of the previous two columns, otherwise NA

[Ref 1] Trautner, J. (Ed.), 2017. Die Laufkäfer Baden-Württembergs. Ulmer Eugen Verlag, Stuttgart (Germany).

[Ref 2] Gesellschaft für Angewandte Carabidologie, e.V.G. Lebens- Raumpräferenzen Der Laufkäfer Deutschlands – Wissensbasierter Katalog, Pp. 1–45. Angewandte Carabidologie, Supplement V., 2009.

[Ref 3] Müller-Motzfeld, Gerd. Die Käfer Mitteleuropas 2. Adephaga 1: Carabidae (Laufkäfer). Edited by H. Freude, K.W. Harde, and G.A. Lohse. 2nd ed. Vol. 2. 11 vols. Krefeld: Goecke u. Evers, 2004.

[Ref 4] Trautner, Jürgen, Michael Bräunicke, Josef Kiechle, Mathias Kramer, Jörg Rietze, Arno Schanowski, and Karin Wolf-Schwenninger. Rote Liste Und Artenverzeichnis Der Laufkäfer Baden-Württembergs (Coleoptera: Carabidae). 3rd ed. Naturschutz-Praxis, Artenschutz 9, 2005.

[Ref 5] Schmidt, J., Jürgen Trautner, and Gerd Müller-Motzfeld. ‘Rote Liste Und Gesamtartenliste Der Laufkäfer (Coleoptera: Carabidae) Deutschlands’. In Rote Liste Der Gefährdeten Tiere, Pflanzen Und Pilze Deutschlands, edited by H. Gruttke, S. Balzer, M. Binot-Hafke, H. Haupt, N. Hofbauer, G. Ludwig, G. Matzke-Hajek, and M. Ries, Vol. Band 4: Wirbellose Tiere (Teil 2). Naturschutz Und Biologische Vielfalt, 70 (4). Bonn: Bundesamt für Naturschutz, 2016.

------ Metadata raw_trap_counts.csv ------------------

- plot_id, date_set, duration: see above on all_data_carabids_df.csv
- trap: trap location - CE (plot center), NW (northwestern quadrant) or SE (southeastern quadrant)
- carabidae_n: total number of carabid individuals captured in each trap
- abax_ovalis to licinus_hoffmanseggii: number of individuals captured for each species, per trap
- comments: explanations of why several traps were not used in further analyses

------ R session information -------------------------

R version 4.2.2 (2022-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19045)

Matrix products: default

locale:
[1] LC_COLLATE=Portuguese_Portugal.utf8  LC_CTYPE=Portuguese_Portugal.utf8    LC_MONETARY=Portuguese_Portugal.utf8
[4] LC_NUMERIC=C                         LC_TIME=Portuguese_Portugal.utf8    

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] MuMIn_1.47.5    ggeffects_1.2.3 DHARMa_0.4.6    MASS_7.3-58.2   lme4_1.1-31     Matrix_1.5-3    car_3.1-2      
 [8] carData_3.0-5   vegan_2.6-4     lattice_0.20-45 permute_0.9-7   iNEXT_3.0.0     corrplot_0.92   gridExtra_2.3  
[15] metR_0.14.0     ggrepel_0.9.3   scales_1.2.1    cowplot_1.1.1   lubridate_1.9.2 forcats_1.0.0   stringr_1.5.0  
[22] dplyr_1.1.2     purrr_1.0.1     readr_2.1.4     tidyr_1.3.0     tibble_3.2.1    ggplot2_3.4.2   tidyverse_2.0.0
[29] plyr_1.8.8     

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.10       utf8_1.2.3        R6_2.5.1          backports_1.4.1   stats4_4.2.2      pillar_1.9.0      rlang_1.1.1   
 [8] rstudioapi_0.14   minqa_1.2.5       data.table_1.14.8 nloptr_2.0.3      checkmate_2.2.0   splines_4.2.2     munsell_0.5.0
[15] compiler_4.2.2    pkgconfig_2.0.3   mgcv_1.8-41       tidyselect_1.2.0  fansi_1.0.4       tzdb_0.4.0        withr_2.5.0  
[22] nlme_3.1-162      gtable_0.3.3      lifecycle_1.0.3   magrittr_2.0.3    cli_3.6.1         stringi_1.7.12    cachem_1.0.8 
[29] reshape2_1.4.4    generics_0.1.3    vctrs_0.6.3       boot_1.3-28.1     tools_4.2.2       glue_1.6.2        hms_1.1.3    
[36] abind_1.4-5       parallel_4.2.2    fastmap_1.1.1     timechange_0.2.0  colorspace_2.1-0  cluster_2.1.4     memoise_2.0.1  


