
# Compare Models on Simulated Missing DT ----------------------------------

# MCMC &
# other models [ mice, Qrilc, missForest ]  -------------------------------

# set parameters
ParamOtherModels <- within(list(), {
  Mice <- list(m = 5, print = FALSE, seed = 6075655L, remove.collinear = FALSE, diagnostic = 1L)
  MissForest <- list(parallelize = "forests", ntree = 100, verbose = FALSE)
  Qrilc <- list(seed = 36075655L)
})

# impute simulated/ test data
SimDTOtherModels <- impute_non_mcmc(SimMissDT, parList = ParamOtherModels)


# mcmc model --------------------------------------------------------------


# Load previously imputed data_SamplesByMCMC
# Loads Imputed Samples and Fitted Models for UPS1 and Yeast .RData

SimDTObjMCMC <- list.files(
  path = here::here("samples_intermediate_outputData"),
  pattern = "SamplesMCMC",
  full.names = TRUE
)

for (.objects in SimDTObjMCMC) {
  load(.objects)
}


# Process imputation results by models ------------------------------------


# mcmc samples/ imputed values
SimDTSamplMCMC <- c(SimDTYeastSamplesMCMC,SimDTUPSamplesMCMC)

SimCombinedModels <-  multi_tidy_merge(SimDTSamplMCMC, SimMissDT, SimDTOtherModels)


# #####################################
#  Qrilc Model
# #####################################

DtNames <- names(TestDTLogs)

QrilcNMS <- vars_select(names(SimCombinedModels), starts_with("Qrilc"))
QrilcDt <- SimCombinedModels[QrilcNMS] 

QrilcModel <- purrr::map(
  DtNames,
  \(dtnm) list_update_status (
    completed = QrilcDt,
    missing = SimMissDT,
    dtname = dtnm,
    non.miss = TestDTLogs,
    sim = TRUE
  )
)

names(QrilcModel) <- DtNames
QrilcModel <- purrr::list_flatten(QrilcModel) %>% dplyr::bind_rows(.id = "data_name") %>% separate(
  col = "data_name",
  into = c("data_name", "model", NA, NA),
  sep = "_"
)

# ######################################
# MissForest Model
# #######################################

MissForestNMS <- vars_select(names(SimCombinedModels), starts_with("MissForest"))
MissForestDt <- SimCombinedModels[MissForestNMS] 

MissForestModel <- purrr::map(
  DtNames,
  \(dtnm) list_update_status (
    completed = MissForestDt,
    missing = SimMissDT,
    dtname = dtnm,
    non.miss = TestDTLogs,
    sim = TRUE
  )
)

names(MissForestModel) <- DtNames
MissForestModel <- purrr::list_flatten(MissForestModel) %>% dplyr::bind_rows(.id = "data_name") %>% separate(
  col = "data_name",
  into = c("data_name", "model", NA, NA),
  sep = "_"
)

# ######################################
# MICE Model
# #######################################

MiceNMS <- vars_select(names(SimCombinedModels), starts_with("Mice"))
MiceDt <- SimCombinedModels[MiceNMS] 

MiceModel <- purrr::map(
  DtNames,
  \(dtnm) list_update_status (
    completed = MiceDt,
    missing = SimMissDT,
    dtname = dtnm,
    non.miss = TestDTLogs,
    sim = TRUE
  )
)

names(MiceModel) <- DtNames
MiceModel <- purrr::list_flatten(MiceModel) %>% dplyr::bind_rows(.id = "data_name") %>% separate(
  col = "data_name",
  into = c("data_name", "model", NA, NA),
  sep = "_"
)

# ######################################
# MCMC Model
# #######################################

MCMCNMS <- vars_select(names(SimCombinedModels), starts_with("MCMC"))
MCMCDt <- SimCombinedModels[MCMCNMS] 


MCMCModel <- purrr::map(
  DtNames,
  \(dtnm) list_update_status (
    completed = MCMCDt,
    missing = SimMissDT,
    dtname = dtnm,
    non.miss = TestDTLogs,
    sim = TRUE
  )
)

names(MCMCModel) <- DtNames
MCMCModel <- purrr::list_flatten(MCMCModel) %>% dplyr::bind_rows(.id = "data_name") %>% separate(
  col = "data_name",
  into = c("data_name", "model", NA, NA),
  sep = "_"
)


SimDTModelOutcomes <- dplyr::bind_rows(QrilcModel, MissForestModel, MiceModel, MCMCModel)

# exporting tibble of imputed results
save(
  SimDTModelOutcomes,
  file = here::here("samples_intermediate_outputData", "SimDTModelOutcomes.RData")
)



