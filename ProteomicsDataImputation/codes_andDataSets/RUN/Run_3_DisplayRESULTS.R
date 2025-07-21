# #################################################
# Loading Data & Packages and Functions
# ################################################

suppressPackageStartupMessages(
  source(file = here::here("codes_andDataSets", "NOT_RUN", "NOT_RUN_LoadingDT.R"))
)

# DATA DESCRIPTION --------------------------------------------------------


# A. overall missing proportions by original data sets 
purrr::map_df(dataFiles$originalDT,function(.data)
{
  dplyr::select(.data, where(is.numeric) ) %>% 
  naniar::prop_miss() 
})

# #######################################################################################
# B. missing proportions by data sets and experimental conditions
# #######################################################################################

( miss_props_originalDT <- list_gather_summarise(dataFiles$originalDT , what = "missing" ) )

# exporting table summary of total counts of features, and proportions of missing features
write_csv(miss_props_originalDT,
          here::here('samples_intermediate_outputData','table_miss_props_originalDT.csv' ) )



# observed transformed data
TidyTestDTLogs <- TidyTestDTLogs %>%
  dplyr::bind_rows(.id = "data.name") %>% 
  dplyr::mutate(data.groups = forcats::fct_collapse(.data[["data.name"]], Yeast = c("YeastMq","YeastPg"), UPS1 = c("UPS1Mq","UPS1Pg") ) )
  

# Distribution of the test data ( complete data with no missing values)
# both in original and log transformation scales

TidyTestDTLogs  %>%
  # dplyr::bind_rows(.id = "data.name") %>% 
  # dplyr::mutate(data.groups = case_when(data.name %in% c("UPS1Mq", "UPS1Pg") ~ "UPS1", .default = "Yeast" )) %>% 
  dplyr::rename(values = log_values) %>% 
  dplyr::group_by(data.groups) %>%
  do(LogScale = ggplot_DensityCurve( ., x = values, facet1 = data.name, facet2 = Samples, 
                                     title = paste0(unique(.data$data.groups)," Transformed Data") ),
      OrigScale = mutate( ., values = negate_log.transform(.data$values, norm.const ) ) %>%
       ggplot_DensityCurve( ., values, data.name, Samples, title = paste0(unique(.data$data.groups)," Data in Original Scale" ) )  ) %>%
  
  tibble::as_tibble() %>% dplyr::ungroup() %>%
  dplyr::mutate(extractedLogScalePlot = mapply( extract_plots, .$LogScale, paste0(.$data.groups,'LogDensity'), dir = "demo_plots",SIMPLIFY = FALSE ),
                extraxtedOrigScalePlot = mapply( extract_plots,.$OrigScale, paste0(.$data.groups,'OriginalDensity'), dir = "demo_plots", SIMPLIFY = FALSE )  )

(miss_props_CompleteData <- list_gather_summarise(TestDTLogs , what = "complete" ) )


# RESULT OUTPUTS ----------------------------------------------------------


# loading imputed data for analysis 
load(here::here('samples_intermediate_outputData','SimDTModelOutcomes.RData') )

SimOutcomess <- dplyr::mutate(
  SimDTModelOutcomes,
  prop.miss = as.numeric(prop.miss),
  data.groups = forcats::fct_collapse(
    data_name,
    Yeast = c("YeastMq", "YeastPg"),
    UPS1 = c("UPS1Mq", "UPS1Pg")
  )
)


# B. comparing distributions of observed vs imputed data [log scale] --------


# Tidy SimOutcomes DT
other.names <- names( dplyr::select(SimOutcomess, vroom::everything(), -c("status","imputed", "original") ) ) 
# here imputed & original are untidy numerical variables that we would like to tidy under a common numeric variable (log_values)
# and a key = status (imputed/ original -Observed)

# Filter imputed values only
# Drop the status variable. Since only imputed values are filtered.
# Tidy/ re-shaping data values under new binary status column, original & imputed
# the numerical columns (imputed and original) are coalesced (tidyed) under new column log_values

PlotingDT <- dplyr::select(filter(SimOutcomess, status == 'imputed'), -status) %>%
  tidyr::gather(key = status,
                value = log_values,
                -tidyselect::any_of(other.names))

PlotingDT$prop.miss <- paste0(PlotingDT$prop.miss * 100, " %")
# reordering levels of column prop.miss
PlotingDT$prop.miss <- factor(PlotingDT$prop.miss, levels = unique(PlotingDT$prop.miss))



# Density Plots -----------------------------------------------------------

PlotingDT %>%
  dplyr::filter(prop.miss %in% c("2 %", "4 %", "6 %", "15 %", "20 %", "25 %")) %>%
  dplyr::group_by(data.groups, model) %>%
  do (
    densityCurves = ggplot_CompareDensityCurves(
      .,
      xVar = log_values,
      groupVar = status,
      facetVar1 = data_name,
      facetVar2 = prop.miss,
      titleVar = model
    )
  ) %>%
  ungroup() %>%
  dplyr::mutate(
    extractdensityCurves = mapply(
      extract_plots,
      .$densityCurves,
      paste0(.$data.groups, .$model, 'DensityPlot'),
      dir = "demo_plots" ,
      SIMPLIFY = FALSE
    )
  )


# C. imputation error -----------------------------------------------------


 
  SimModelError <- SimOutcomess %>%
  dplyr::filter(status == "imputed") %>%
  dplyr::group_by(across(c(data.groups, data_name, prop.miss, model))) %>%
  dplyr::summarise(
    MAE = magrittr::add(boot_mean(abs.error, iter = 10000, seed = 15), 1),
    NRMSE = magrittr::add(boot_normalised_RMSE(original, imputed, iter = 10000, seed = 15), 1), .groups = "drop"
  )

# line plots for errors

dplyr::filter(SimModelError, model != "Qrilc") %>%
  dplyr::group_by(data.groups) %>%
  do(
    PlotNRMSE = ggplot_ErrorLines(., prop.miss, NRMSE, model, facetVar = data_name, titleVar = unique(.$data.groups)),
    PlotMAE = ggplot_ErrorLines(., prop.miss, MAE, model, facetVar = data_name, titleVar = unique(.$data.groups))
  ) %>%
  tibble::as_tibble() %>%
  dplyr::mutate(
    extractNRMSE = mapply(extract_plots, .$PlotNRMSE, paste0(.$data.groups, "NRMSE"), dir= "demo_plots",SIMPLIFY = FALSE),
    extractMAE = mapply(extract_plots, .$PlotMAE, paste0(.$data.groups, "MAE"), dir= "demo_plots", SIMPLIFY = FALSE)
  )


# Variance Component Graphs -----------------------------------------------

SimModelVar <- SimOutcomess %>% 
            dplyr::group_by(data_name, prop.miss, model ) %>% 
            do (imputedData = var(.data$imputed) ,  observedData = var(.data$original)  ) %>% 
  
            unnest(cols = c('imputedData', 'observedData') ) %>% 
  
            tidyr:: gather( key = Status, value = Variance, c('imputedData', 'observedData') ) %>% 
            dplyr::mutate( groupedData = as.factor( if_else(model != 'Qrilc', 'ComparableModels', 'QrilcModel') ),
                           nameFiles = data_name )

# I want to be consistent here.
# So, I will update the models variable with Observed as values from the status variable                      
SimModelVar[SimModelVar$Status== 'imputedData','Status'] = SimModelVar[SimModelVar$Status =='imputedData', 'model']

# variance plot for ; MCMC, missForest and mice models

SimModelVar %>%
  split(.$groupedData) %>% # this line ensures plotting Qrilc model separately. Bcs it's not comparable to other models
  purrr::map(~ dplyr::group_by(., data_name) %>%
    purrrlyr::by_slice(~ ggplot_VarLines(.x, xVar = prop.miss, yVar = Variance, groupVar = Status, titleVar = nameFiles), .to = "plotedCurves")) %>%
  purrr::map2_df(., names(.), ~ update_list(.x, groupedData = .y)) %>%
  dplyr::select(data_name, groupedData, vroom::everything()) %>%
  dplyr::mutate(extractedPlots = mapply(extract_plots, .$plotedCurves, paste0(.$data_name, .$groupedData, "VarPlots"), dir="demo_plots",SIMPLIFY = FALSE))



