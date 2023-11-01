Code author: João M. Cordeiro Pereira
Affiliation: Chair of Wildlife Ecology & Management, University of Freiburg, Germany

Publication title: Specialist carabid beetles in mixed montane forests show positive links to roe deer and to biodiversity-oriented forestry
Authors: João M. Cordeiro Pereira, Sebastian Schwegmann, Clàudia Massó Estaje, Martin Denter, Grzegorz Mikusinski, Ilse Storch

Last updated: 1. Nov 2023

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


