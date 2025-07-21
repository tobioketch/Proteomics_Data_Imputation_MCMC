# Load Data And Packages --------------------------------------------------

suppressPackageStartupMessages(
  source(file = here::here("codes_andDataSets", "NOT_RUN", "NOT_RUN_LoadingDT.R"))
)


# #####################################################################
# 
# MCMC IMPUTATION BEGINS HERE WITH YEAST DATA SET
# 
# #####################################################################


# mcmc impute simulated yeast  --------------------------------------------

# script for mcmc_imputation
sampling_script <- paste0(here::here("codes_andDataSets", "NOT_RUN", "NOT_RUN_mcmc_scriptFILE.stan") )


# subset simulated missing, Yeast data
name_missYeastDT <- names(SimMissDT) %>% 
  tidyselect::vars_select(., contains('Yeast') )

SimMissYeastDT <- SimMissDT[ name_missYeastDT]

# Sim Yeast MCMC Model

rstan_options(auto_write = TRUE)
plan(multisession, workers = parallel::detectCores())

SimDTYeastModelMCMC <- furrr::future_map(SimMissYeastDT, function(dt) {
  impute_mcmc(dt, file = sampling_script)
})

plan(sequential)

SimDTYeastSamplesMCMC <- list_summarise_mcmc(SimDTYeastModelMCMC, "miss_data")

# Save Yeast Imps MCMC 
save(
  SimDTYeastSamplesMCMC,
  file = here::here("samples_intermediate_outputData", "SimDTYeastSamplesMCMC.RData")
)

rm(SimDTYeastModelMCMC, SimDTYeastSamplesMCMC)

# warnings()


# ##########################################################################
# 
#  HERE WE IMPUTE UPS1 DATA SET BY CALLING A SCRIPT
# 
# ##########################################################################

# mcmc impute simulated UPS1 ----------------------------------------------

# Here we source the file for imputing UPS1 Dataset
source(file = here::here("codes_andDataSets", "NOT_RUN", "NOT_RUN_SimDT_ImpUPS_MCMC.R"))











