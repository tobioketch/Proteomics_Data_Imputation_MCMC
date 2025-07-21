

# ******************************************************************************************
#  Preferred basic set ups to facilitate imputation process

#  main parameters for MCMC runs

n_chain <- 4         # number of MCMC_rstan chains
n_iter <- 2000        # Total iterations for MCMC_rstan
iter_warmUps <- 1000  # warm_up iterations 

# ******************************************************************************************
# Packages and Functions
suppressWarnings(suppressPackageStartupMessages(
  source(file = here::here("codes_andDataSets", "NOT_RUN", "NOT_RUN_Packages_andFunctions.R"))
) )


## Loading Data Files 
fileList <- list.files( path = here::here('codes_andDataSets','dataSets'), pattern = 'csv$' )

# Removing extensions from the file name
fileNames <- sub("\\.csv", "", fileList) 

# Loading & naming the list of loaded data data
dataFiles <- tibble(fileList, fileNames) %>%
  dplyr::mutate(originalDT = purrr::map(fileList, ~ read_csv(here::here("codes_andDataSets", "dataSets", .),
    na = c("0", "---"),
    col_types = cols(.default = col_guess())
  )))

dataFiles$originalDT = setNames(dataFiles$originalDT, fileNames)


# #######################################################################################################
# 
# SIMULATED MISSINGNESS (MAR)
# 
# ########################################################################################################


# Simulates Test Data Set -------------------------------------------------

# ******************************************************************************************
# 
# Amounts of missing
amount.miss <- c(0.02, 0.04, 0.08, 0.1, 0.15, 0.20, 0.25 )
# 
# ******************************************************************************************


# ##########################################################################################
# Test data and transformations
# Test data : Observed (complete ) data, excluding missing values
# ###########################################################################################
# set.seed(100)
TestDT <- dataFiles$originalDT %>%
  purrr::map(function(x)na.exclude(x) ) #%>% sample_n(20)) 



# Log transformation of test, observed data
# determining normalizing constant for log transformation 

norm.const = TestDT %>%
  purrr::map(~ norm_const_log(.[, sapply(., class) == "numeric"], na.rm = TRUE)) %>% reduce(sum)

TestDTLogs <- TestDT %>%
  purrr::map(~ purrr::lmap_if(.x, is.numeric, log_transform, norm.const))

TidyTestDTLogs <- TestDTLogs %>%
  purrr::map(~ is_numeric_gather(.x, key = "Samples", value = "log_values") %>%
               dplyr::mutate(Samples = readr::parse_number(Samples)) %>%  tibble::as_tibble())


# Missing data simulation

SimMissDT <- suppressMessages(purrr::map(amount.miss, function(x) {
  purrr::map(
    TestDTLogs, function(dt) {
      sim_miss_by_reps(dt, prop.miss = x, seed = 10208888) %>%
        dplyr::mutate(prop.miss = paste0(x)) 
    }
  )
}) )  #%>% unlist(., recursive = FALSE))

# renaming the simulated missing data 

names(SimMissDT) <- amount.miss
SimMissDT <- list_flatten(SimMissDT, name_spec = "{inner}_{outer}")













