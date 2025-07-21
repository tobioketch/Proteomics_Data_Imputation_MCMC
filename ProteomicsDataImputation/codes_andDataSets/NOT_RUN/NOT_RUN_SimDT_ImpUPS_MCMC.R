# #########################################################################################
# 
# UPS1 DATA SET - MCMC IMPUTATION
# 
# ##########################################################################################


# script for mcmc_imputation
sampling_script <- paste0(here::here("codes_andDataSets", "NOT_RUN", "NOT_RUN_mcmc_scriptFILE.stan") )

## UPS1 Data
## subset simulated missing, Yeast data
name_missUPS1DT <- names(SimMissDT) %>%
  tidyselect::vars_select(., contains('UPS1'))

SimMissUPS1DT <- SimMissDT[name_missUPS1DT]


rstan_options(auto_write =TRUE)

# filter simulated missing data by amount missing
names_miss1 <- tidyselect::vars_select(
  names(SimMissUPS1DT),
  vroom::num_range("UPS1Mq_", seq(0.02, 0.08, 0.02)),
  vroom::num_range("UPS1Pg_", seq(0.02, 0.08, 0.02))
)

names_miss2 <- setdiff(names(SimMissUPS1DT), names_miss1)

# #######################################################################################
# 
# impute rstan_MCMC model -UPS1 data
# impute missing amount between 2 - 8 %
# 
# #######################################################################################

plan(multisession, workers = parallel::detectCores())

SimDTUPSModelMCMC..1 <- furrr::future_map(SimMissUPS1DT[names_miss1], function(dt) {
  impute_mcmc(dt, file = sampling_script)
})

plan(sequential)

# First Samples
SimDTUPSamplesMCMC..1 <- list_summarise_mcmc(SimDTUPSModelMCMC..1, 'miss_data')

# ######################################################################################
# 
# impute missing amount between 10% + (plus)
# 
# ######################################################################################

plan(multisession, workers = parallel::detectCores())
SimDTUPSModelMCMC..2 <- furrr::future_map(SimMissUPS1DT[names_miss2], function(dt) {
  impute_mcmc(dt, file = sampling_script)
})
plan(sequential)

# Second Samples
SimDTUPSamplesMCMC..2 <- list_summarise_mcmc( SimDTUPSModelMCMC..2, 'miss_data')

# Combining UPS1 Samples MCMC
SimDTUPSamplesMCMC <- c(SimDTUPSamplesMCMC..1, SimDTUPSamplesMCMC..2) 

save(
  SimDTUPSamplesMCMC,
  file = here::here("samples_intermediate_outputData", "SimDTUPSamplesMCMC.RData")
)

rm(
  SimDTUPSModelMCMC..1, SimDTUPSModelMCMC..2, SimDTUPSamplesMCMC..1, SimDTUPSamplesMCMC..2
)

# warnings()


