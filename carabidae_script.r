# Paper title: Specialist carabid beetles in mixed montane forests show positive links to roe deer and to biodiversity-oriented forestry
# Authors: João M. Cordeiro Pereira, Sebastian Schwegmann, Clàudia Massó, Martin Denter, Grzegorz Mikusinski, Ilse Storch
# code written by: João M. Cordeiro Pereira

#-------------Session info:--------------------#
# R version 4.2.2 (2022-10-31 ucrt)            #
# Platform: x86_64-w64-mingw32/x64 (64-bit)    #
# Running under: Windows 10 x64 (build 19045)  #
#----------------------------------------------#

# all file paths are relative to directory of project "shared_carabids.Rproj"

#------1. load packages--------

require(plyr) # needed for functions of package s3cR (and needs to be loaded before dplyr and tidyverse, to avoid conflicts)
require(tidyverse) # includes dplyr and ggplot2, for data wrangling and plotting
require(cowplot) # for a nice ggplot2 theme
require(scales) # for manipulation of ggplot2 scales
require(forcats) # for manipulating factors
require(ggrepel) # to avoid overlaps in text labels on ggplot
require(metR) # to add labels to ordisurf contours in ggplot in an easier way
require(grid) # for arranging plots on a grid
require(gridExtra) # add-on to package "grid"
require(corrplot) # to examine multi-collinearities
require(iNEXT) # for species accummulation curves and richness estimation
require(vegan) # for ordination analyses
require(car) # an calculate VIF for predictors in a model
require(lme4) # required for GLMM
require(MASS) # for negative binomial model (glm.nb)
require(DHARMa) # for model diagnostics based on simulations
require(ggeffects) # to extract model predictions, which can be plotted with ggplot2
require(MuMIn) # to calculate pseudo-R-squared of GLMM

#------2. functions from GitHub package "s3cR" (Gaüzére 2018, https://github.com/pgauzere/s3cR)-------------------------------

# only functions cwmean() and cwvar(), used for calculating community-weighted mean and variance of carabid body sizes
# Gaüzère, Pierre, Guilhem Doulcier, Vincent Devictor, and Sonia Kéfi. ‘A Framework for Estimating Species-Specific Contributions to Community Indicators’. Ecological Indicators 99 (April 2019): 74–82. https://doi.org/10.1016/j.ecolind.2018.11.069.

cwmean <-
  function(df, trait_val_col="trait_val"){
    
    #   Compute the community weighted mean of the dataframe.
    # 
    #     Args:
    #         df (dataframe)       : Containing a column "n" (number of individuals) and a column referencing the trait value.
    #         trai_val_col (text)  : names of the column referencing trait value 
    # 
    #     Return: 
    #         Community weighted mean.
    
    return(  sum( df[,trait_val_col] *  df$n / sum(df$n))) }

cwvar <-
  function(df, trait_val_col="trait_val", bessel=TRUE, cwm=NA){
    #   Compute the community weighted mean and community weighted variance of the dataframe.
    # 
    #     Args:
    #         df (dataframe)       : Containing a column "n" (number of individuals) and a column referencing the trait value.
    #         trai_val_col (text)  : names of the column referencing trait value 
    #         bessel (bool)        : If TRUE, use bessel correction for an unbiased variance estimator (N/(N-1)).
    #         cwm                  : If Community weighted means has already been computed (it is only to speed up computations)
    #     Return: 
    #         Community weighted mean and community weighted variance.
    
    df$trait_val<-df[[trait_val_col]]
    
    if (is.na(cwm) == T) {cwm <- cwmean(df)}
    
    if (bessel==TRUE){
      corrective_term = sum(df$n)  / (sum(df$n) -1)
    } else {
      corrective_term = 1
    }
    
    n2  <- sum ( df$trait_val^2 * df$n / sum(df$n))
    return(cwv=corrective_term * (n2 - cwm^2))}

#----- 3. Import and process data -------------------------------------------------

all_data_carabids <- read.table("alldata_carabids_df.csv", sep = ",", header = TRUE)
# check data structure
str(all_data_carabids)
head(all_data_carabids)

sum(all_data_carabids$carabidae_n) # total number of Carabidae specimens

# add plot-level species richness (carabid_sr)
names(all_data_carabids) # first check column indices
all_data_carabids <- all_data_carabids %>% dplyr::mutate(carabid_sr = rowSums(.[9:49] > 0))

# import data with functional group classification, and check data structure
carabid_traits <- read.table("carabid_traits.csv", sep = ",", header = TRUE)
head(carabid_traits)

# create vectors of species names for each functional group of interest
sp_forspec <- carabid_traits$species[carabid_traits$forest_specialist == "yes"]
sp_monspec <- carabid_traits$species[carabid_traits$montane_specialist == "yes"]
sp_brachy <- carabid_traits$species[carabid_traits$flight == "brachypterous"]
sp_large <- carabid_traits$species[carabid_traits$body_length > 22]
sp_nonfor <- carabid_traits$species[carabid_traits$forest_specialist == "no"] # non-forest species (for section 5)

# for forest species, forest specialists, montane specialists and brachypterous species, calculate abundance and species richness
all_data_carabids <- all_data_carabids %>%
  dplyr::mutate(carabid_forspec_n = rowSums(.[names(all_data_carabids) %in% sp_forspec]),
         carabid_forspec_sr = rowSums(.[names(all_data_carabids) %in% sp_forspec] > 0),
         carabid_monspec_n = rowSums(.[names(all_data_carabids) %in% sp_monspec]),
         carabid_monspec_sr = rowSums(.[names(all_data_carabids) %in% sp_monspec] > 0),
         carabid_brachy_n = rowSums(.[names(all_data_carabids) %in% sp_brachy]),
         carabid_brachy_sr = rowSums(.[names(all_data_carabids) %in% sp_brachy] > 0),
         carabid_large_n = rowSums(.[names(all_data_carabids) %in% sp_large]),
         carabid_large_sr = rowSums(.[names(all_data_carabids) %in% sp_large] > 0),
         carabid_non_for_n = rowSums(.[names(all_data_carabids) %in% sp_nonfor])
  )

sum(all_data_carabids$carabid_forspec_n) # total forest specialist individuals
sum(all_data_carabids$carabid_monspec_n) # total montane specialist individuals
sum(all_data_carabids$carabid_brachy_n) # total brachypterous individuals
sum(all_data_carabids$carabid_large_n) # total large individuals
# calculate corresponding percentages (of total number of carabid specimens)
(sum(all_data_carabids$carabid_forspec_n)/sum(all_data_carabids$carabidae_n))*100
(sum(all_data_carabids$carabid_monspec_n)/sum(all_data_carabids$carabidae_n))*100
(sum(all_data_carabids$carabid_brachy_n)/sum(all_data_carabids$carabidae_n))*100
(sum(all_data_carabids$carabid_large_n)/sum(all_data_carabids$carabidae_n))*100

# create table with total abundances per species
species_totals <- all_data_carabids %>% dplyr::summarise(across(abax_ovalis:licinus_hoffmanseggii, sum)) %>%
  pivot_longer(cols = everything(), names_to = "species", values_to = "n")
# add percentages of dominance
species_totals <- species_totals %>% dplyr::mutate(dominance = (n/sum(species_totals$n))*100)
# add body size to this table
species_totals <- species_totals %>% dplyr::left_join(carabid_traits[c(1,5)], by = "species")

# convert body lengths to long format, so you can plot a body size distribution over all individuals
species_totals_replicated <- rep(species_totals$body_length, species_totals$n)
# plot body size distribution, obtain mean and median
hist(species_totals_replicated)
mean(species_totals_replicated)
median (species_totals_replicated)

# create table with frequencies of each species across all plots
species_freqs <- all_data_carabids %>% dplyr::summarise(across(abax_ovalis:licinus_hoffmanseggii, ~ sum(. > 0))) %>%
  pivot_longer(cols = everything(), names_to = "species", values_to = "freq")

# calculate community weighted means, variance and skewness, using functions from package s3cR (code section 1)
# first we need to extract a dataframe containing only plot_id and the species' abundances
community_df_cwv <- all_data_carabids %>% dplyr::select(plot_id, abax_ovalis:licinus_hoffmanseggii)
# add empty columns for community-weighted means (cwm) and community-weighted variances (cwv)
community_df_cwv$cwm <- NA
community_df_cwv$cwv <- NA

# loop across each row, to apply the cwmean() and cwvar() functions to populate the empty columns
# these functions require a single-site dataframe with columns for species, their counts and corresponding trait value
# so, on each iteration of the loop, we need to reshape the data to long format, join in the body length value from the carabid_traits dataframe, retrieve cwm and cwv values, and store them on the i-th row of our initial dataframe
for (i in 1:nrow(community_df_cwv)){
  
  community_df_cwv$cwm[i] <- pivot_longer(community_df_cwv[i,], cols = abax_ovalis:licinus_hoffmanseggii, names_to = "species", values_to = "n") %>%
  as.data.frame() %>%
  left_join(carabid_traits[c(1,5)], by = "species") %>%
  cwmean(trait_val_col = "body_length")
  
  community_df_cwv$cwv[i] <- pivot_longer(community_df_cwv[i,], cols = abax_ovalis:licinus_hoffmanseggii, names_to = "species", values_to = "n") %>%
    as.data.frame() %>%
    left_join(carabid_traits[c(1,5)], by = "species") %>%
    cwvar(trait_val_col = "body_length")
}

head(community_df_cwv) # check if community-weighted means and variances were correctly added (last two columns of dataframe)

# add these variables back to the main data frame
all_data_carabids <- all_data_carabids %>% left_join(community_df_cwv[c(1,43,44)], by = "plot_id")

# inspect distribution of response variables
hist(all_data_carabids$carabidae_n, breaks = 15) # Activity-density of all carabids
hist(all_data_carabids$carabid_sr, breaks = 15) # Species richness of all carabids
hist(all_data_carabids$carabid_brachy_n, breaks = 15) # Activity-density of brachypterous carabids
hist(all_data_carabids$carabid_brachy_sr, breaks = 15) # Species richness of brachypterous carabids
hist(all_data_carabids$carabid_large_n, breaks = 15) # Activity-density of large carabids
hist(all_data_carabids$carabid_monspec_n, breaks = 15) # Activity-density of montane specialist carabids
hist(all_data_carabids$carabid_monspec_sr, breaks = 15) # Species richness of montane specialist carabids
hist(all_data_carabids$cwm) # Community-weighted mean of body size
hist(all_data_carabids$cwv) # Community-weighted variance of body size

# check whether cwm and cwv are correlated
cor.test(all_data_carabids$cwm, all_data_carabids$cwv, method = "pearson")

# inspect distribution of predictors
hist(all_data_carabids$for_cover_100ha, breaks = 15)
hist(all_data_carabids$for_cover_2500ha, breaks = 15)
hist(all_data_carabids$avg_alt, breaks = 15)
hist(all_data_carabids$canopycover, breaks = 15)
hist(all_data_carabids$roe_deer, breaks = 15)
hist(all_data_carabids$DBHMean, breaks = 15)
hist(all_data_carabids$lying_dw_volume, breaks = 15)
hist(all_data_carabids$pc_broadleaf, breaks = 15)
hist(all_data_carabids$sd_slope, breaks = 15)
hist(all_data_carabids$northness, breaks = 15)

# check for collinearities among predictors
# first check column indices of predictors
names(all_data_carabids) # indices 50 to 59
cor_matrix <- cor(all_data_carabids[c(50:59)], use = "pairwise.complete.obs", method = "spearman")
# visualize correlation matrix, using function corrplot.mixed() from package corrplot
corrplot.mixed(cor_matrix, tl.pos='lt', tl.cex=0.6, number.cex=0.5, addCoefasPercent=T)

# check log-transformations for highly-skewed predictors (lying deadwood and roe deer)
hist(log(all_data_carabids$lying_dw_volume), breaks = 15)
hist(log(all_data_carabids$roe_deer), breaks = 15)
# add log-transformed predictors (with suffix "_tf") to our dataset
all_data_carabids <- all_data_carabids %>% dplyr::mutate(
  lying_dw_volume_tf = log(lying_dw_volume),
  roe_deer_tf = log(roe_deer)
)

# scale all predictors to mean of 0 and standard deviation of 1
# (remember to use log-transformed roe deer and lying deadwood data)
# scaled predictors receive suffix "_stdr"
all_data_carabids <- all_data_carabids %>% dplyr::mutate(
  avg_alt_stdr = (avg_alt - mean(avg_alt))/sd(avg_alt),
  DBHMean_stdr = (DBHMean - mean(DBHMean))/sd(DBHMean),
  lying_dw_volume_stdr = (lying_dw_volume_tf - mean(lying_dw_volume_tf))/sd(lying_dw_volume_tf),
  pc_broadleaf_stdr = (pc_broadleaf - mean(pc_broadleaf))/sd(pc_broadleaf),
  roe_deer_stdr = (roe_deer_tf - mean(roe_deer_tf))/sd(roe_deer_tf),
  canopycover_stdr = (canopycover - mean(canopycover))/sd(canopycover),
  for_cover_100ha_stdr = (for_cover_100ha - mean(for_cover_100ha))/sd(for_cover_100ha),
  for_cover_2500ha_stdr = (for_cover_2500ha - mean(for_cover_2500ha))/sd(for_cover_2500ha),
  sd_slope_stdr = (sd_slope - mean(sd_slope))/sd(sd_slope),
  northness_stdr = (northness - mean(northness))/sd(northness)
)

#----- 4. Species accumulation curve using iNext ---------------------------------

# extract community data frame containing only species data (and the number of traps, for now)
community_df <- all_data_carabids %>% dplyr::select(plot_id, n_traps, abax_ovalis:licinus_hoffmanseggii) %>% as.data.frame()
rownames(community_df) <- community_df$plot_id # set plot IDs as row names
community_df$plot_id <- NULL # remove plot ID column

# for iNext, species should be on rows and sites on columns, so transpose the community dataframe
chao_samples <- as.data.frame(t(community_df))
# remove plots for which number of traps is below 3, so variance of estimates is not inflated
chao_samples <- chao_samples %>% dplyr::select(where(~ first(.) == 3))
ncol(chao_samples) # 59 sites kept

# then remove first row (number of traps)
chao_samples <- chao_samples[-1,]
# remove n_traps from community_df as well, since we won't need it anymore
community_df <- community_df %>% dplyr::select(-n_traps)

#for incidence-based (sample-based) analysis, convert to presence/absence
chao_samples[chao_samples > 1] <- 1

#calculate asymptotic Chao richness
ChaoRichness(chao_samples, datatype = "incidence_raw", conf = 0.95)

# calculate incidence-based (sample-based) rarefaction/extrapolation curves
# first convert data to an incidence_freq datatype (incidence_raw is not working well)
chao_samples_freq <- chao_samples %>% mutate(freq = rowSums(.)) %>%
  dplyr::select(freq)
# for incidence_freq data type, the first row must be the total number of samples
samples_n <- data.frame(freq = ncol(chao_samples)) # create that value
chao_samples_freq <- rbind(samples_n, chao_samples_freq) # add it on top of the dataframe
# run iNEXT analysis
inext_samples <- iNEXT(chao_samples_freq, datatype = "incidence_freq", conf = 0.95, nboot = 200)
inext_samples

#save Chao richness estimates and confidence intervals on their own dataframe for plotting
chao_samples_asymp <- inext_samples$AsyEst[1,]

# plot curve for incidence-based (sample-based) richness, using ggiNEXT(), a ggplot2 extension
# gnore warning message
sac2 <- ggiNEXT(inext_samples, type = 1, se = TRUE) +
  labs(x = "Pooled sampling plots", y = "Pooled species richness") +
  geom_pointrange(data = chao_samples_asymp, aes(x = 130, y = Estimator, ymin = LCL, ymax = UCL), inherit.aes = FALSE, size = 1, colour = "green4", linewidth = 1) + #add asymptotic richness and its confidence interval in green
  geom_hline(yintercept = chao_samples_asymp[1,4], linetype = "dashed", color = "green4") + #add reference line for estimated asymptotic richness
  ylim(0,85) +
  theme_cowplot() +
  theme(axis.line = element_line("gray25"), text = element_text(size = 12), legend.position = "none")
sac2

#----- 5. GLMM model for non-forest species ----------------------------------------------------

# run GLMM, with observation-level random effect
m_nonfor <- glmer(cbind(carabid_non_for_n, carabid_forspec_n) ~ (1|plot_id) + canopycover_stdr + for_cover_2500ha_stdr + for_cover_100ha_stdr, data = all_data_carabids, family = binomial(link = "logit"))
summary(m_nonfor)
vif(m_nonfor)

# use likelihood ratio test, to assess significance of fixed effects (more reliable than Wald-z tests included in glmer summary)
# first create nested models removing each of the predictors
m_nonfor_1 <- update(m_nonfor, . ~ . - canopycover_stdr)
m_nonfor_2 <- update(m_nonfor, . ~ . - for_cover_2500ha_stdr)
m_nonfor_3 <- update(m_nonfor, . ~ . - for_cover_100ha_stdr)
# then compare each nested model with the full model (first nested model, then full model)
anova(m_nonfor_1, m_nonfor)
anova(m_nonfor_2, m_nonfor)
anova(m_nonfor_3, m_nonfor)

# check model diagnostics with package DHARMa
simulation_m_nonfor <- simulateResiduals(m_nonfor, plot = TRUE)
plotResiduals(simulation_m_nonfor, form = all_data_carabids$avg_alt_stdr)
plotResiduals(simulation_m_nonfor, form = all_data_carabids$for_cover_100ha_stdr)
plotResiduals(simulation_m_nonfor, form = all_data_carabids$for_cover_100ha_stdr)
testDispersion(simulation_m_nonfor)
testOutliers(simulation_m_nonfor)
testZeroInflation(simulation_m_nonfor)

# get Nakagawa & Schielzeth pseudo-R-squared (look only at marginal one)
MuMIn::r.squaredGLMM(m_nonfor)

mean(all_data_carabids$carabidae_n) # 84.45455 (mean no. carabids per plot)

# generate transformed values for 50%, 70% and 90% canopy cover
# first create function to transform values, then apply it to 0.5, 0.7 and 0.9
canopycover_transf <- function(x){(x - mean(all_data_carabids$canopycover))/sd(all_data_carabids$canopycover)}
canopy_cover_refvalues <- sapply(c(0.5,0.7,0.9), canopycover_transf)

# generate predictions for those values, using ggpredict() from package ggeffects
preds_canopy_cover <- ggpredict(m_nonfor, terms = "canopycover_stdr [canopy_cover_refvalues]")
# predictions must then be multiplied by 84 species (rounded average of no. specimens per plot)
preds_canopy_cover$predicted * 84

#----- 6. Run GLM models and inspect outputs -------------------------------------------------

names(all_data_carabids) # check out variable names

#run quasipoisson model for carabid species richness (full community) (model m1)
m1 <- glm(carabid_sr ~ avg_alt_stdr + roe_deer_stdr + canopycover_stdr + pc_broadleaf_stdr + DBHMean_stdr + lying_dw_volume_stdr + trapdays, data = all_data_carabids, family = quasipoisson(link = "log"))
summary(m1) # model summary and explained deviance
vif(m1) # check variance inflation factors
# model diagnostics
par(mfrow = c(2,2))
plot(m1)
par(mfrow = c(1,1))

# create alternative model without roe deer
m1.1 <- update(m1, . ~ . - roe_deer_stdr)
# test comparison between models with and without roe deer
anova(m1 , m1.1, test="F")

# run negative binomial model for carabid activity-density (full community) (model m2)
m2 <- glm.nb(carabidae_n ~ avg_alt_stdr + roe_deer_stdr + canopycover_stdr + pc_broadleaf_stdr + DBHMean_stdr +  lying_dw_volume_stdr + offset(log(trapdays)), data = all_data_carabids)
summary(m2)
vif(m2)
# check model diagnostics with DHARMa package
simulation_m2 <- simulateResiduals(m2, n = 1000, plot = TRUE) # QQ plot and uniformity of residuals (residuals vs. fitted)
# then residuals vs. each of the predictors
plotResiduals(simulation_m2, form = all_data_carabids$avg_alt_stdr)
plotResiduals(simulation_m2, form = all_data_carabids$roe_deer_stdr)
plotResiduals(simulation_m2, form = all_data_carabids$canopycover_stdr)
plotResiduals(simulation_m2, form = all_data_carabids$pc_broadleaf_stdr)
plotResiduals(simulation_m2, form = all_data_carabids$DBHMean_stdr)
plotResiduals(simulation_m2, form = all_data_carabids$lying_dw_volume_stdr)
testDispersion(simulation_m2) # test for remaining overdispersion
testOutliers(simulation_m2) # test for outliers
testZeroInflation(simulation_m2) # test for zero-inflation

#create alternative model removing roe deer, and conduct model comparison with likelihood ratio test
m2.1 <- update(m2, . ~ . - roe_deer_stdr)
anova(m2.1, m2) # likelihood ratio test, first the nested model and then the complex one

# run equivalent single-species models for all most common species:

# Abax parallelepipedus
m2.ap <- glm.nb(abax_parallelepipedus ~ avg_alt_stdr + roe_deer_stdr + canopycover_stdr + pc_broadleaf_stdr + DBHMean_stdr + lying_dw_volume_stdr + offset(log(trapdays)), data = all_data_carabids)
summary(m2.ap)
simulation_m2.ap <- simulateResiduals(m2.ap, plot = TRUE) # quick model diagnostic with DHARMa

# Pterostichus burmeisteri
m2.pb <- glm.nb(pterostichus_burmeisteri ~ avg_alt_stdr + roe_deer_stdr + canopycover_stdr + pc_broadleaf_stdr + DBHMean_stdr + lying_dw_volume_stdr + offset(log(trapdays)), data = all_data_carabids)
summary(m2.pb)
simulation_m2.pb <- simulateResiduals(m2.pb, plot = TRUE)

# Pterostichus burmeisteri (2nd model, with same selection of predictors as model m6 for montane specialists: without pc_broadleaf, but adding sd_slope and northness)
m6.pb <- glm.nb(pterostichus_burmeisteri ~ avg_alt_stdr + roe_deer_stdr + canopycover_stdr + DBHMean_stdr + lying_dw_volume_stdr + sd_slope_stdr + northness_stdr + offset(log(trapdays)), data = all_data_carabids)
summary(m6.pb)
simulation_m6.pb <- simulateResiduals(m6.pb, plot = TRUE)

# Abax ovalis
m2.ao <- glm.nb(abax_ovalis ~ avg_alt_stdr + roe_deer_stdr + canopycover_stdr + pc_broadleaf_stdr + DBHMean_stdr + lying_dw_volume_stdr + offset(log(trapdays)), data = all_data_carabids)
summary(m2.ao)
simulation_m2.ao <- simulateResiduals(m2.ao, plot = TRUE)

# Carabus auronitens
m2.ca <- glm.nb(carabus_auronitens ~ avg_alt_stdr + roe_deer_stdr + canopycover_stdr + pc_broadleaf_stdr + DBHMean_stdr + lying_dw_volume_stdr + offset(log(trapdays)), data = all_data_carabids)
summary(m2.ca)
simulation_m2.ca <- simulateResiduals(m2.ca, plot = TRUE)
# Carabus auronitens (2nd model, analogous to m6)
m6.ca <- glm.nb(carabus_auronitens ~ avg_alt_stdr + roe_deer_stdr + canopycover_stdr + DBHMean_stdr + lying_dw_volume_stdr + sd_slope_stdr + northness_stdr + offset(log(trapdays)), data = all_data_carabids)
summary(m6.ca)
simulation_m6.ca <- simulateResiduals(m6.ca, plot = TRUE)

# Carabus nemoralis
m2.cn <- glm.nb(carabus_nemoralis ~ avg_alt_stdr + roe_deer_stdr + canopycover_stdr + pc_broadleaf_stdr + DBHMean_stdr + lying_dw_volume_stdr + offset(log(trapdays)), data = all_data_carabids)
summary(m2.cn)
simulation_m2.cn <- simulateResiduals(m2.cn, plot = TRUE)

# Molops piceus
m2.mp <- glm.nb(molops_piceus ~ avg_alt_stdr + roe_deer_stdr + canopycover_stdr + pc_broadleaf_stdr + DBHMean_stdr + lying_dw_volume_stdr + offset(log(trapdays)), data = all_data_carabids)
summary(m2.mp)
simulation_m2.mp <- simulateResiduals(m2.mp, plot = TRUE)

# run quasipoisson model for species richness of brachypterous carabids (model m3)
m3 <- glm(carabid_brachy_sr ~ avg_alt_stdr + roe_deer_stdr + canopycover_stdr + pc_broadleaf_stdr + DBHMean_stdr + lying_dw_volume_stdr + trapdays, data = all_data_carabids, family = quasipoisson(link = "log"))
summary(m3)
par(mfrow = c(2,2))
plot(m3)
par(mfrow = c(1,1))

# run negative binomial model for activity-density of brachypterous carabids, and run DHARMa diagnostics (model m4)
m4 <- glm.nb(carabid_brachy_n ~ avg_alt_stdr + roe_deer_stdr + canopycover_stdr + pc_broadleaf_stdr + DBHMean_stdr + lying_dw_volume_stdr + offset(log(trapdays)), data = all_data_carabids)
summary(m4)
simulation_m4 <- simulateResiduals(m4, plot = TRUE)

# run quasipoisson model for species richness of montane specialist carabids (model m5)
m5 <- glm(carabid_monspec_sr ~ avg_alt_stdr + roe_deer_stdr + canopycover_stdr + DBHMean_stdr + lying_dw_volume_stdr + sd_slope_stdr + northness_stdr + trapdays, data = all_data_carabids, family = quasipoisson(link = "log"))
summary(m5)
par(mfrow = c(2,2))
plot(m5)
par(mfrow = c(1,1))

# run negative binomial model for activity-density of montane specialist carabids, and run DHARMa diagnostics (model m6)
m6 <- glm.nb(carabid_monspec_n ~ avg_alt_stdr + roe_deer_stdr + canopycover_stdr + DBHMean_stdr + lying_dw_volume_stdr + sd_slope_stdr + northness_stdr + offset(log(trapdays)), data = all_data_carabids)
summary(m6)
simulation_m6 <- simulateResiduals(m6, plot = TRUE)

# run negative binomial model for activity-density of large carabids, and run DHARMa diagnostics (model m7)
m7 <- glm.nb(carabid_large_n ~ avg_alt_stdr + roe_deer_stdr + canopycover_stdr + pc_broadleaf_stdr + DBHMean_stdr + lying_dw_volume_stdr + offset(log(trapdays)), data = all_data_carabids)
summary(m7)
simulation_m7 <- simulateResiduals(m7, plot = TRUE)

# run linear model for CWM body size
m8 <- lm(cwm ~  avg_alt_stdr + roe_deer_stdr + canopycover_stdr + pc_broadleaf_stdr + DBHMean_stdr + lying_dw_volume_stdr + trapdays, data = all_data_carabids)
summary(m8)
par(mfrow = c(2,2))
plot(m8)

# run linear model for CWV body size
m9 <- lm(cwv ~  avg_alt_stdr + roe_deer_stdr + canopycover_stdr + pc_broadleaf_stdr + DBHMean_stdr + lying_dw_volume_stdr + trapdays, data = all_data_carabids, na.action = "na.fail")
summary(m9)
plot(m9)
par(mfrow = c(1,1))

#----- 7. Plot significant modelled relationships  (with ggplot2) -------------------------------------------------

power_e <- function(x) round(exp(x), digits = 3) # define function to transform x axis log values back to original values (used for roe deer and lying deadwood)

# extract mean and sd of all relevant predictors, so they can be used to rescale x axis values to the original predictor values
mean_alt <- mean(all_data_carabids$avg_alt)
sd_alt <- sd(all_data_carabids$avg_alt)
mean_roe <- mean(all_data_carabids$roe_deer_tf)
sd_roe <- sd(all_data_carabids$roe_deer_tf)
mean_dbh <- mean(all_data_carabids$DBHMean)
sd_dbh <- sd(all_data_carabids$DBHMean)
mean_broad <- mean(all_data_carabids$pc_broadleaf)
sd_broad <- sd(all_data_carabids$pc_broadleaf)
mean_cover <- mean(all_data_carabids$canopycover)
sd_cover <- sd(all_data_carabids$canopycover)
mean_dw <- mean(all_data_carabids$lying_dw_volume_tf)
sd_dw <- sd(all_data_carabids$lying_dw_volume_tf)

# for model 1: effect of altitude on total species richness

# using ggpredict() from package ggeffects, generate model predictions from 500 values of altitude
predicted_df1_alt <- ggpredict(m1, terms = "avg_alt_stdr [n = 500]")

ggsr_alt <- ggplot(all_data_carabids, aes(x = avg_alt, y = carabid_sr)) +
  geom_point(cex = 2, alpha = 0.7) + # add data points (with observed values)
  labs(x = "Altitude (m a.s.l.)", y = "Total SR") + #axis labels
  geom_line(aes(x = x*sd_alt + mean_alt, y = predicted), data = predicted_df1_alt, colour = "darkgreen", linewidth = 0.8) + # add line between all predicted values
  geom_ribbon(aes(x = x*sd_alt + mean_alt, y = predicted, ymin = conf.low, ymax = conf.high), data = predicted_df1_alt, alpha = 0.1) + # add ribbon with confidence intervals for predicted values
  theme_cowplot() # use nice-looking cowplot theme
ggsr_alt

# for model 1: effect of log roe deer on total species richness

predicted_df1_roe <- ggpredict(m1, terms = "roe_deer_stdr [n = 500]")

ggsr_deer <- ggplot(all_data_carabids, aes(x = roe_deer_tf, y = carabid_sr)) +
  geom_point(cex = 2, alpha = 0.7) +
  labs(x = "Relative abundance of roe deer (log scale)", y = "Total SR") +
  geom_line(aes(x = x*sd_roe + mean_roe, y = predicted), data = predicted_df1_roe, colour = "darkgreen", linewidth = 0.8) +
  geom_ribbon(aes(x = x*sd_roe + mean_roe, y = predicted, ymin = conf.low, ymax = conf.high), data = predicted_df1_roe, alpha = 0.1) +
  scale_x_continuous(labels = power_e) + # reverse log transformation on X axis labels
  theme_cowplot()
ggsr_deer

# for model 1: effect of DBH mean on total species richness

predicted_df1_dbh <- ggpredict(m1, terms = "DBHMean_stdr [n = 500]")

ggsr_dbh <- ggplot(all_data_carabids, aes(x = DBHMean/10, y = carabid_sr)) + # convert DBH from mm to cm
  geom_point(cex = 2, alpha = 0.7) +
  labs(x = "Mean tree DBH (cm)", y = "Total SR") +
  geom_line(aes(x = (x*sd_dbh + mean_dbh)/10, y = predicted), data = predicted_df1_dbh, colour = "darkgreen", linewidth = 0.8) +
  geom_ribbon(aes(x = (x*sd_dbh + mean_dbh)/10, y = predicted, ymin = conf.low, ymax = conf.high), data = predicted_df1_dbh, alpha = 0.1) +
  theme_cowplot()
ggsr_dbh

# for model 2: effect of altitude on total carabid abundance

predicted_df2_alt <- ggpredict(m2, terms = "avg_alt_stdr [n = 500]")

ggab_alt <- ggplot(all_data_carabids, aes(x = avg_alt, y = (carabidae_n/trapdays)*mean(trapdays))) + # standardize sampling effort by mean trap-days
  geom_point(cex = 2, alpha = 0.7) +
  labs(x = "Altitude (m a.s.l.)", y = "Total AD") +
  geom_line(aes(x = x*sd_alt + mean_alt, y = predicted), data = predicted_df2_alt, colour = "darkgreen", linewidth = 0.8) +
  geom_ribbon(aes(x = x*sd_alt + mean_alt, y = predicted, ymin = conf.low, ymax = conf.high), data = predicted_df2_alt, alpha = 0.1) +
  theme_cowplot()
ggab_alt

# for model 2: effect of log roe deer on total carabid abundance

predicted_df2_roe <- ggpredict(m2, terms = "roe_deer_stdr [n = 500]")

ggab_deer <- ggplot(all_data_carabids, aes(x = roe_deer_tf, y = (carabidae_n/trapdays)*mean(trapdays))) +
  geom_point(cex = 2, alpha = 0.7) +
  labs(x = "Relative abundance of roe deer (log scale)", y = "Total AD") +
  geom_line(aes(x = x*sd_roe + mean_roe, y = predicted), data = predicted_df2_roe, colour = "darkgreen", linewidth = 0.8) +
  geom_ribbon(aes(x = x*sd_roe + mean_roe, y = predicted, ymin = conf.low, ymax = conf.high), data = predicted_df2_roe, alpha = 0.1) +
  scale_x_continuous(labels = power_e) +
  theme_cowplot()
ggab_deer

# for model 2: effect of pc_broadleaf on total carabid abundance

predicted_df2_broad <- ggpredict(m2, terms = "pc_broadleaf_stdr [n = 500]")

ggab_broad <- ggplot(all_data_carabids, aes(x = pc_broadleaf*100, y = (carabidae_n/trapdays)*mean(trapdays))) + # convert pc_broadleaf from fraction to percentage
  geom_point(cex = 2, alpha = 0.7) +
  labs(x = "% Broadleaf", y = "Total AD") +
  geom_line(aes(x = (x*sd_broad + mean_broad)*100, y = predicted), data = predicted_df2_broad, colour = "darkgreen", linewidth = 0.8) +
  geom_ribbon(aes(x = (x*sd_broad + mean_broad)*100, y = predicted, ymin = conf.low, ymax = conf.high), data = predicted_df2_broad, alpha = 0.1) +
  theme_cowplot()
ggab_broad

# for model 3: effect of mean DBH on brachypterous beetle richness

predicted_df3_dbh <- ggpredict(m3, terms = "DBHMean_stdr [n = 500]")

ggbrachysr_dbh <- ggplot(all_data_carabids, aes(x = DBHMean/10, y = carabid_brachy_sr)) +
  geom_point(cex = 2, alpha = 0.7) +
  labs(x = "Mean tree DBH (cm)", y = "Brach. spp. SR") +
  geom_line(aes(x = (x*sd_dbh + mean_dbh)/10, y = predicted), data = predicted_df3_dbh, colour = "darkgreen", linewidth = 0.8) +
  geom_ribbon(aes(x = (x*sd_dbh + mean_dbh)/10, y = predicted, ymin = conf.low, ymax = conf.high), data = predicted_df3_dbh, alpha = 0.1) +
  scale_y_continuous(breaks = c(3, 6, 9, 12)) +
  theme_cowplot()
ggbrachysr_dbh

# for model 4: effect of altitude on brachypterous beetle abundance

predicted_df4_alt <- ggpredict(m4, terms = "avg_alt_stdr [n = 500]")

ggbrachy_alt <- ggplot(all_data_carabids, aes(x = avg_alt, y = (carabid_brachy_n/trapdays)*mean(trapdays))) +
  geom_point(cex = 2, alpha = 0.7) +
  labs(x = "Altitude (m a.s.l.)", y = "Brach. spp. AD") +
  geom_line(aes(x = x*sd_alt + mean_alt, y = predicted), data = predicted_df4_alt, colour = "darkgreen", linewidth = 0.8) +
  geom_ribbon(aes(x = x*sd_alt + mean_alt, y = predicted, ymin = conf.low, ymax = conf.high), data = predicted_df4_alt, alpha = 0.1) +
  theme_cowplot()
ggbrachy_alt

# for model 4: effect of roe deer on brachypterous beetle abundance

predicted_df4_roe <- ggpredict(m4, terms = "roe_deer_stdr [n = 500]")

ggbrachy_deer <- ggplot(all_data_carabids, aes(x = roe_deer_tf, y = (carabid_brachy_n/trapdays)*mean(trapdays))) +
  geom_point(cex = 2, alpha = 0.7) +
  labs(x = "Relative abundance of roe deer (log scale)", y = "Brach. spp. AD") +
  geom_line(aes(x = x*sd_roe + mean_roe, y = predicted), data = predicted_df4_roe, colour = "darkgreen", linewidth = 0.8) +
  geom_ribbon(aes(x = x*sd_roe + mean_roe, y = predicted, ymin = conf.low, ymax = conf.high), data = predicted_df4_roe, alpha = 0.1) +
  scale_x_continuous(labels = power_e) +
  theme_cowplot()
ggbrachy_deer

# for model 4: effect of pc_broadleaf on brachypterous beetle abundance

predicted_df4_broad <- ggpredict(m4, terms = "pc_broadleaf_stdr [n = 500]")

ggbrachy_broad <- ggplot(all_data_carabids, aes(x = pc_broadleaf*100, y = (carabid_brachy_n/trapdays)*mean(trapdays))) +
  geom_point(cex = 2, alpha = 0.7) +
  labs(x = "% Broadleaf", y = "Brach. spp. AD") +
  geom_line(aes(x = (x*sd_broad + mean_broad)*100, y = predicted), data = predicted_df4_broad, colour = "darkgreen", linewidth = 0.8) +
  geom_ribbon(aes(x = (x*sd_broad + mean_broad)*100, y = predicted, ymin = conf.low, ymax = conf.high), data = predicted_df4_broad, alpha = 0.1) +
  theme_cowplot()
ggbrachy_broad

#for model 5: effect of altitude on montane specialist richness

predicted_df5_alt <- ggpredict(m5, terms = "avg_alt_stdr [n = 500]")

ggmonsr_alt <- ggplot(all_data_carabids, aes(x = avg_alt, y = carabid_monspec_sr)) +
  geom_point(cex = 2, alpha = 0.7) +
  labs(x = "Altitude (m a.s.l.)", y = "Montane SR") +
  geom_line(aes(x = x*sd_alt + mean_alt, y = predicted), data = predicted_df5_alt, colour = "darkgreen", linewidth = 0.8) +
  geom_ribbon(aes(x = x*sd_alt + mean_alt, y = predicted, ymin = conf.low, ymax = conf.high), data = predicted_df5_alt, alpha = 0.1) +
  theme_cowplot()
ggmonsr_alt

#for model 5: effect of canopy cover on montane specialist richness

predicted_df5_cover <- ggpredict(m5, terms = "canopycover_stdr [n = 500]")

ggmonsr_can <- ggplot(all_data_carabids, aes(x = canopycover*100, y = carabid_monspec_sr)) + # convert canopy cover from fraction to percentage
  geom_point(cex = 2, alpha = 0.7) +
  labs(x = "% Canopy Cover", y = "Montane SR") +
  geom_line(aes(x = (x*sd_cover + mean_cover)*100, y = predicted), data = predicted_df5_cover, colour = "darkgreen", linewidth = 0.8) +
  geom_ribbon(aes(x = (x*sd_cover + mean_cover)*100, y = predicted, ymin = conf.low, ymax = conf.high), data = predicted_df5_cover, alpha = 0.1) +
  theme_cowplot()
ggmonsr_can

#for model 5: effect of lying dead wood on montane specialist richness

predicted_df5_dw <- ggpredict(m5, terms = "lying_dw_volume_stdr [n = 500]")

ggmonsr_dw <- ggplot(all_data_carabids, aes(x = lying_dw_volume_tf, y = carabid_monspec_sr)) +
  geom_point(cex = 2, alpha = 0.7) +
  labs(x = "Lying deadwood volume (m3, log scale)", y = "Montane SR") +
  geom_line(aes(x = x*sd_dw + mean_dw, y = predicted), data = predicted_df5_dw, colour = "darkgreen", linewidth = 0.8) +
  geom_ribbon(aes(x = x*sd_dw + mean_dw, y = predicted, ymin = conf.low, ymax = conf.high), data = predicted_df5_dw, alpha = 0.1) +
  scale_x_continuous(breaks = c(log(10), log(20), log(30), log(40), log(50), log(100), log(150), log(200), log(300)), labels = exp) + #adjust X axis labels to invert log transformation
  theme_cowplot()
ggmonsr_dw

#for model 6: effect of altitude on montane specialist abundance

predicted_df6_alt <- ggpredict(m6, terms = "avg_alt_stdr [n = 500]")

ggmona_alt <- ggplot(all_data_carabids, aes(x = avg_alt, y = (carabid_monspec_n/trapdays)*mean(trapdays))) +
  geom_point(cex = 2, alpha = 0.7) +
  labs(x = "Altitude (m a.s.l.)", y = "Montane spp. AD") +
  geom_line(aes(x = x*sd_alt + mean_alt, y = predicted), data = predicted_df6_alt, colour = "darkgreen", linewidth = 0.8) +
  geom_ribbon(aes(x = x*sd_alt + mean_alt, y = predicted, ymin = conf.low, ymax = conf.high), data = predicted_df6_alt, alpha = 0.1) +
  theme_cowplot()
ggmona_alt

#for model 6: effect of canopy cover on montane specialist abundance

predicted_df6_cover <- ggpredict(m6, terms = "canopycover_stdr [n = 500]")

ggmona_can <- ggplot(all_data_carabids, aes(x = canopycover*100, y = (carabid_monspec_n/trapdays)*mean(trapdays))) +
  geom_point(cex = 2, alpha = 0.7) +
  labs(x = "% Canopy Cover", y = "Montane spp. AD") +
  geom_line(aes(x = (x*sd_cover + mean_cover)*100, y = predicted), data = predicted_df6_cover, colour = "darkgreen", linewidth = 0.8) +
  geom_ribbon(aes(x = (x*sd_cover + mean_cover)*100, y = predicted, ymin = conf.low, ymax = conf.high), data = predicted_df6_cover, alpha = 0.1) +
  theme_cowplot()
ggmona_can

# for model 7: effect of altitude on abundance of large carabids

predicted_df7_alt <- ggpredict(m7, terms = "avg_alt_stdr [n = 500]")

gglarge_alt <- ggplot(all_data_carabids, aes(x = avg_alt, y = (carabid_large_n/trapdays)*mean(trapdays))) +
  geom_point(cex = 2, alpha = 0.7) +
  labs(x = "Altitude (m a.s.l.)", y = "Large spp. AD") +
  geom_line(aes(x = x*sd_alt + mean_alt, y = predicted), data = predicted_df7_alt, colour = "darkgreen", linewidth = 0.8) +
  geom_ribbon(aes(x = x*sd_alt + mean_alt, y = predicted, ymin = conf.low, ymax = conf.high), data = predicted_df7_alt, alpha = 0.1) +
  theme_cowplot()
gglarge_alt

# for model 9: effect of pc_broadleaf on body size variance

predicted_df9_broad <- ggpredict(m9, terms = "pc_broadleaf_stdr [n = 500]")

ggcwv_broad <- ggplot(all_data_carabids, aes(x = pc_broadleaf*100, y = cwv)) +
  geom_point(cex = 2, alpha = 0.7) +
  labs(x = "% Broadleaf", y = "CWV body size") +
  geom_line(aes(x = (x*sd_broad + mean_broad)*100, y = predicted), data = predicted_df9_broad, colour = "darkgreen", linewidth = 0.8) +
  geom_ribbon(aes(x = (x*sd_broad + mean_broad)*100, y = predicted, ymin = conf.low, ymax = conf.high), data = predicted_df9_broad, alpha = 0.1) +
  theme_cowplot()
ggcwv_broad

# for model 9: effect of DBH mean on body size variance

predicted_df9_dbh <- ggpredict(m9, terms = "DBHMean_stdr [n = 500]")

ggcwv_dbh <- ggplot(all_data_carabids, aes(x = DBHMean/10, y = cwv)) +
  geom_point(cex = 2, alpha = 0.7) +
  labs(x = "Mean tree DBH (cm)", y = "CWV body size") +
  geom_line(aes(x = (x*sd_dbh + mean_dbh)/10, y = predicted), data = predicted_df9_dbh, colour = "darkgreen", linewidth = 0.8) +
  geom_ribbon(aes(x = (x*sd_dbh + mean_dbh)/10, y = predicted, ymin = conf.low, ymax = conf.high), data = predicted_df9_dbh, alpha = 0.1) +
  theme_cowplot()
ggcwv_dbh

#----------- 8. Create plot grids with all significant relationships -------------------------

# first create text graphic objects, with the names of the predictors, which can be used as "row names" in our plot grid
# \n can be used to create line breaks, just = "left" for aligning text and font size can also be specified with gp = gpar()
t1 <- textGrob("Altitude \n(m)", just = "left", gp = gpar(fontsize = 18))
t2 <- textGrob("Roe deer \n(rel. abd.)", just = "left", gp = gpar(fontsize = 18))
t3 <- textGrob("Broadleaf \n(%)", just = "left", gp = gpar(fontsize = 18))
t4 <- textGrob("Canopy cover \n(%)", just = "left", gp = gpar(fontsize = 18))
t5 <- textGrob("Mean DBH \n(cm)", just = "left", gp = gpar(fontsize = 18))
t6 <- textGrob("Lying DW \n(m3)", just = "left", gp = gpar(fontsize = 18))

# create a list including text graphical objects and ggplot objects, in the right order
relationship_plots <- list(t1, ggsr_alt, ggab_alt, ggbrachy_alt,
                       gglarge_alt, ggmonsr_alt, ggmona_alt,
                       t2, ggsr_deer, ggab_deer, ggbrachy_deer,
                       t3, ggab_broad, ggbrachy_broad, ggcwv_broad,
                       t4, ggmonsr_can, ggmona_can,
                       t5, ggsr_dbh, ggbrachysr_dbh, ggcwv_dbh,
                       t6, ggmonsr_dw)

# we need to remove the X axis titles and set all Y axis titles to the same font size
# these adjustments can be applied to all ggplots using map_if() from package purrr
relationship_plots <- relationship_plots %>% map_if(is.ggplot,
                                                    ~.x + labs(x = NULL) + theme(axis.title.y = element_text(size = 18)))

# create a plot grid using grid.arrange()
# argument grobs takes a list of graphical objects (ggplots, text grobs, etc.)
# layout matrix can be used to place the elements of the list (by their indices) on their intended location
# since the figure can't fit in a single page, we created one grid for half of the plots (indices 1 to 15), and another for the other half (indices 16 to 24)
fig3.1 <- grid.arrange(grobs = relationship_plots[1:15], nrow = 4, ncol = 4,
                     widths = c(1.5,2,2,2),
                     layout_matrix = rbind(c(1,2,3,4),
                                   c(NA,5,6,7),
                                   c(8,9,10,11),
                                   c(12,13,14,15))
                     )
fig3.2 <- grid.arrange(grobs = relationship_plots[16:24], nrow = 3, ncol = 4,
                       widths = c(1.5,2,2,2),
                       layout_matrix = rbind(c(16,17,18,NA),
                                             c(19,20,21,22),
                                             c(23,24,NA,NA))
                       )

# export plot grids with correct resolution using ggsave()
ggsave("fig3.1.png", fig3.1, width = 15, height = 15, dpi = 600, bg = "white")
ggsave("fig3.2.png", fig3.2, width = 15, height = 11, dpi = 600, bg = "white")

#----------- 9. NMDS analysis -----------------------------------

# package vegan requires a community dataframe, only containing species counts (and plot IDs only as row names)
# we will use the community_df created above, in section 3
View(community_df)
# convert community data to presence-absence
community_df[community_df > 1] <- 1

# run NMDS, with function metaMDS() from package vegan
set.seed(3)
nmds <- metaMDS(comm = community_df, distance = "bray", k = 2, try = 10000, tidy = TRUE)

# check output and stress diagnostics
nmds
stressplot(nmds)
plot(nmds)

# conduct null hypothesis testing to check strength of ordination (as per Dexter et al. 2018)
# null hypothesis is the absence of community structure
# data is simulated with absence of structure (using in this case a "quasiswap_count" algorithm, where row and column totals, as well as presences and absences, are maintained)
# stress value of observed data is tested against 1000 (recommended) stress values from null model permutations
stressTest <- oecosimu(comm = community_df, method = "quasiswap", 
                     nestfun = metaMDS, k = 2, distance = "bray",
                     nsimul = 1000,
                     parallel = 2,
                     statistic = "stress", 
                     alternative = "less",
                     trace = TRUE,
                     maxit = 1000, trymax = 200, sratmax = 0.9999999)
stressTest # check test output

# plot outcome of null-model testing
hist(as.vector(stressTest$oecosimu$simulated), xlim = c(0, max(stressTest$oecosimu$simulated) + .05), xlab = "Ecological null model stress value", ylab = "Frequency of stress value", main = "", breaks = 7)
abline(v = stressTest$oecosimu$statistic, col = "red", lty = 2)

# extract predictor data with plot names as row names for nmds post-hoc fitting
predictors_df <- all_data_carabids %>% dplyr::select(plot_id, avg_alt, DBHMean, lying_dw_volume, pc_broadleaf, roe_deer, canopycover, for_cover_100ha, sd_slope, northness) %>% as.data.frame()
rownames(predictors_df) <- predictors_df$plot_id # set plot IDs as row names
predictors_df$plot_id <- NULL # remove plot ID column

# do post-hoc fitting of predictors and permutation tests, using envfit() from package vegan
set.seed(5)
nmds_envfit <- envfit(nmds ~ avg_alt + for_cover_100ha + log(roe_deer) + canopycover + pc_broadleaf + DBHMean + log(lying_dw_volume), data = predictors_df, permutations = 10000)
nmds_envfit

# calculate Bonferroni-corrected significance threshold
bonferroni_threshold <- 0.05/length(nmds_envfit$vectors$pvals)
# plot NMDS with significant post-hoc fitted predictors
plot(nmds)
plot(nmds_envfit, p.max = bonferroni_threshold) # overlays vector on NMDS plot
# plot smoothed surfaces for elevation, using function ordisurf()
ordisurf_altitude_nmds <- ordisurf(scores(nmds), predictors_df$avg_alt, col = "black", add = TRUE, levels = c(500, 600, 700, 800, 900, 1000, 1100, 1200, 1300), lwd.cl = 0.5)
summary(ordisurf_altitude_nmds)

# plot NMDS with ggplot2

# first extract sites and species scores
site_scores_nmds <- as.data.frame(scores(nmds, display = "sites"))
species_scores_nmds <- as.data.frame(scores(nmds, display = "species"))
names(all_data_carabids)

# extract variable scores from envfit
variable_scores_nmds <- as.data.frame(scores(nmds_envfit, display = "vectors"))
variable_scores_nmds$pvals <- nmds_envfit$vectors$pvals # add p-values
variable_scores_nmds <- rownames_to_column(variable_scores_nmds, var = "predictor") # transfer variable names from row names to a column of their own
# rename altitude variable, to use as label when plotting
variable_scores_nmds[1,1] <- "Altitude"

# join altitude values to site scores
site_scores_nmds <- site_scores_nmds %>% rownames_to_column(var = "plot_id") %>%
  left_join(all_data_carabids[c(1,50)], by = "plot_id")

# join carabid trait data to species scores
species_scores_nmds <- species_scores_nmds %>% rownames_to_column(var = "species") %>%
  left_join(carabid_traits, by = "species")

# create vector of shortened species names (in the same sequence as species scores)
short_names <- c("Aova","Apap","Apal","Caut","Caur","Ccor","Cint","Cirr","Cnem","Cpro","Csyl","Cvio","Carc","Catt",
                 "Ccar","Dagi","Hlae","Hlat","Lruf","Mela","Mpic","Nbre","Nsal","Nbig","Nruf","Pnig","Paet","Pbur",
                 "Pobl","Pmad","Pmel","Pcri","Pnit","Tlae","Tnit","Pcup","Ppum","Pdil","Lass","Cmic","Lhof")

species_scores_nmds$species_short <- short_names #add short names to species scores data frame

# NMDS plot with sites
# to incorporate ordisurf contours, we use code from https://chrischizinski.github.io/rstats/ordisurf/)

# first add a z column to site scores, that can be filled with the values for contours
site_scores_nmds$z <- NA

# extract x, y and z that form contour lines, from the ordisurf object
ordisurf_altitude_nmds$grid
# get all combinations of x and y coordinates
ordisurf_contours <- expand.grid(x = ordisurf_altitude_nmds$grid$x, y = ordisurf_altitude_nmds$grid$y)
# join corresponding z values (altitudes) for every combination of x and y
ordisurf_contours <- cbind(ordisurf_contours, c(ordisurf_altitude_nmds$grid$z))
names(ordisurf_contours) <- c("x","y","z") # set variable names
head(ordisurf_contours)

# create ggplot
# use stat_contour to draw contour lines, based on the ordisurf_contours dataframe
# for the contour labels, use geom_text_contour() from package metR, it's really useful!
# in geom_text_contour(), the same data and aesthetics must be used as in stat_contour
# additionally, the "skip" argument ensures no contour is skipped in labelling, the "stroke" that the number is not crossed off by the contour line, and "label.placer" that the label angles are the most similar possible (but can also be used to adjust number of labels per contour, etc.)
nmds_plot_sites <- ggplot(site_scores_nmds, aes(x = NMDS1, y = NMDS2, colour = avg_alt, label = plot_id)) +
  geom_point(size = 3) + # first just add the points (colour is already specified in the plot aesthetics)
  geom_text(data = site_scores_nmds[21,], aes(colour = NULL), nudge_y = 0.05) + # then text labels for two specific outlying sites
  geom_segment(data = variable_scores_nmds[variable_scores_nmds$pvals < bonferroni_threshold,], aes(x = 0, xend = NMDS1, y= 0, yend = NMDS2), inherit.aes = FALSE, arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd = 1) + # add arrow (override plot aesthetics), but subset data for significant variables
  geom_hline(yintercept = 0, linetype = 2, colour = "grey10") + # reference line for y = 0
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10") + # reference line for x = 0
  stat_contour(data = ordisurf_contours, aes(x = x, y = y, z = z), inherit.aes = FALSE, colour = "black", alpha = 0.7, breaks = c(750, 800, 850, 900, 950, 1000, 1050)) + # add contours from ordisurf
  metR::geom_text_contour(data = ordisurf_contours, aes(x = x, y = y, z = z), inherit.aes = FALSE, stroke = 0.15, cex = 4, skip = 0, label.placer = label_placer_flattest(ref_angle = 45)) + # add contour labels
  scale_color_viridis_c() + # replace colour scale for viridis gradient, with better perceptuality
  labs(colour = "Altitude (m)", fill = NULL) + # specify title for altitude legend, and remove the title for polygon legend
  theme_cowplot() +
  theme(legend.title = element_text(size=12))
nmds_plot_sites # ignore warning message, it refers to the NAs on the contour lines dataframe

#plot species (including contour lines for altitude, but not including contour labels, so the plot is not overcrowded)
nmds_plot_species <- ggplot(species_scores_nmds, aes(x = NMDS1, y = NMDS2, colour = redlisted)) +
  geom_text_repel(aes(label = species_short), size = 5, max.overlaps = 30) +
  stat_contour(data = ordisurf_contours, aes(x = x, y = y, z = z), inherit.aes = FALSE, colour = "darkblue", alpha = 0.7, breaks = c(750, 800, 850, 900, 950, 1000, 1050)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey10") + #reference line for 0
  geom_vline(xintercept = 0, linetype = 2, colour = "grey10") + #reference line for 0
  theme_cowplot() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("red"), na.value = "black")
nmds_plot_species

# -------------- end of code ---------------------
