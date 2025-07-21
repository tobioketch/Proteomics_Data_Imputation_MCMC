# installing and loading packages

# Start by :
# Install an appropriate version of Rtools
# Go HERE https://cran.r-project.org/bin/windows/Rtools/

# Installing BioConductor Packages ----------------------------------------

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if(!require('impute'))
  BiocManager::install("impute", force = TRUE)

if (!require('pcaMethods'))
  BiocManager::install("pcaMethods", force = TRUE)



# Loading Packages --------------------------------------------------------

if(!require('pacman'))
  install.packages('pacman')
suppressPackageStartupMessages(
  pacman::p_load(
    "magrittr",
    "MASS",
    "readr",
    "dplyr",
    "rlang",
    "tibble",
    "tidyr",
    "tidyselect",
    "purrr",
    "stringr",
    "glue",
    "ggplot2",
    "parallel",
    "devtools",
    "randomForest",
    "missForest",
    "ggpubr",
    "curl",
    "lattice",
    "mice",
    "rstan",
    "fabricatr",
    "doParallel",
    "naniar",
    "vroom",
    "gridExtra",
    "ggeasy",
    "furrr",
    "plan",
    "multidplyr",
    "fs",
    "purrrlyr",
    "tidyselect",
    "here",
    "corrr",
    "future",
    "imputeLCMD"
  )
)



# Project Functions -------------------------------------------------------

# Simulating normal data:  default column-wise univariate normal data
# simulate multivariate normal data if  mv.norm = TRUE

# sample.data : incomplete data set whose distribution we are learning
# seed : a numerical number to initialize a random sampling process

sim_rnorm<-function(sample.data,seed = NULL, mv.rnorm = FALSE)
{
  if(missing(seed))( seed = 5 ) 
  set.seed(seed)
  num.data = na.omit( dplyr::select_if(sample.data,is.numeric ) ) 
  N = nrow(sample.data)
  mean_vec = summarise_all( num.data, .funs = mean) %>% as.matrix()
  cov_mat = as.matrix( cov(num.data) )
  # multivariate simulation
  if(mv.rnorm)( sim.data = data.frame(MASS::mvrnorm( n = N , mean_vec, cov_mat, empirical = TRUE )) )
  else { 
    sim.data = setNames(vector("list", length = length(mean_vec)), names(num.data) )
    # univariate simulation
    for (i in seq_len( length(mean_vec) )){
      sim.data[[i]] = rnorm(n = N, mean = mean_vec[i], sd= sqrt(cov_mat[i,i]))
    }
    sim.data = data.frame(do.call("cbind", sim.data))
  }
  bind_cols(sim.data, sample.data[ , setdiff( names(sample.data), names(num.data) ) ] )[names(sample.data) ]
}

# Simulates complete data
# learning from the observed features in the missing data  
# miss.data   : a data frame with missing values

# out puts    : completed sample missing values are replaced with rnorm sample
#             : bootstrap sample - observed features are resample 
#             : rnorm.sample - independent univariate normal (mvrnorm if mv.norm = TRUE)

sim_completeData<-function(miss, seed = NULL, mv.rnorm = FALSE){
  complet = data.frame( miss )
  ind.miss = which(is.na(miss), arr.ind = TRUE )
  nMiss = length(ind.miss)
  nCases = nrow(miss)
  nSamples = ncol(miss )
  pMiss = nMiss/prod(nCases, nSamples )
  
  obs = na.omit(miss)
  if(missing(seed)) (seed = 5)
  set.seed(seed)    # randomizing sampling set
  bootsamp = sample_n( obs, size = nrow(obs) ) %>% 
    fabricatr::resample_data(., N = nCases )  # simulating complete test data
  rnorm.sample = sim_rnorm( bootsamp, seed, mv.rnorm )
  for ( i in 1:nrow(ind.miss) ){
    complet[ind.miss[i,1], ind.miss[i,2]] = rnorm.sample[ind.miss[i,1], ind.miss[i,2]] 
  }
  return(list( complet = complet , prop.miss = pMiss, 
               bootstrap = bootsamp , normal = rnorm.sample ) )
}

# Takes data frame with replicates of samples as column major variables and
# gathers (stretch out ) the replicates columns (numeric columns) into 
# a single column of Reps (key), and values; 

# Additional columns created to identify corresponding samples (conditions);
# and featureIds.Samples

# featureIds.Samples can be used to compute row/ feature summaries by samples ( row/ feature means by samples) 

# I. Integrated Option

pivot_long_reps <- function(.data, cols_vary = "fastest") {
  temp <- .data
  non_num_nms <- names(dplyr::select(temp, where(Negate(is.numeric))))
  temp <- temp %>%
    tibble::rowid_to_column(., var = "FeatureId") %>%
    pivot_longer(
      cols = -c(FeatureId, tidyselect::all_of(non_num_nms)),
      names_to = "Reps", values_to = "Value", cols_vary = cols_vary
    )
  # Sample/ Conditions/ Factors
  mutate(temp,
         Condition = as.character(readr::parse_number(Reps)),
         SampFeatureId = paste(Condition, FeatureId, sep = "_")
  ) %>%
    dplyr::select(tidyselect::all_of(non_num_nms), SampFeatureId, Reps, Condition, Value)
}

# 2. Improved - Not integrated yet

pivot_long_reps_v2 <- function(dt, into = c(NA, "condition", NA), sep = "_", ...) {
  non_numeric_nms <- dplyr::select_if(dt, Negate(is.numeric)) %>% names()
  new_vars <- into[which(!is.na(into))]
  res <- dt %>%
    tibble::rowid_to_column(., var = "row_id") %>%
    dplyr::mutate(row_id = as.character(row_id)) %>%
    is_numeric_gather(., key = "samples", value = "values") %>%
    tidyr::separate(samples, into = into, remove = FALSE, convert = TRUE, sep = sep, ...)
  repl_grp <- purrr::map(new_vars, ~ dplyr::pull(res, .x))
  repl_grp <- do.call("paste", c(repl_grp, sep = "_"))
  
  res %>%
    dplyr::mutate(repl_grp = paste(repl_grp, .data$row_id, sep = "_")) %>%
    dplyr::select(
      tidyselect::all_of(non_numeric_nms),
      vroom::everything(), -row_id
    )
  res
}




# Takes data frame with replicates of samples as column major variables (wide data frame)  and 
# uses numeric variables (replicates of samples/ conditions) to compute 
# rowmeans (means of features) grouped by samples/ trt conditions

feature_averages <- function(data, digits = 3L) {
  nums <- dplyr::select(data, where(is.numeric))
  longDt <- pivot_long_reps(data, cols_vary = "slowest")
  avgsDt <- dplyr::group_by(longDt, SampFeatureId) %>%
    mutate(rowAvgs = round(mean(Value, na.rm = TRUE), digits = digits)) %>%
    dplyr::ungroup()
  avgsDt <- matrix(avgsDt$rowAvgs, nrow(nums), ncol(nums)) %>% as.data.frame(.)
  colnames(avgsDt) <- colnames(nums)
  avgsDt
}


# Takes data with replicates as column major variables (wide data frame ) and 
# coalesces replicates into samples/ conditions major variables
# .data : a data frame with replicates as variables
# Id :    identifier variable
# drop.Id : default is TRUE ( Identifier variable not returned )
# as_experiment_vars
as_experiment_vars<-function(.data, Id = NULL, drop.Id = TRUE){
  nums = dplyr::select(.data, where(is.numeric) )
  grps = stringr::str_sub(names(nums), start = 1, end = -3 ) %>% unique()
  if(!missing(Id)){
    temp = dplyr::bind_cols( dplyr::select(.data, Id), nums ) %>% 
      tidyr::gather(.k, .v, -Id )
  }
  if(missing(Id)){
    temp = nums %>% 
      tibble::rowid_to_column(., "Id") %>% tidyr::gather(., .k, .v, -Id)
  }
  temp = temp %>% 
    dplyr::mutate(.k = str_sub(.$.k, start = 1, end = -3) ) %>%  split(.$.k)
  # re-arranging the list by sample groups, select and bind cols of values
  Out = purrr::map(temp[grps], ~dplyr::select(., .v) ) %>% reduce(cbind)
  names(Out) = grps 
  if(drop.Id)return(Out)
  if(!drop.Id)return(dplyr::bind_cols(dplyr::select(temp[[1]], Id), Out ) )
}

# inserts /or creates MAR pattern on a complete data set
insert_missingValues <- function(complete.data, frac.miss, seed = NULL ){
  num.data = dplyr::select_if( complete.data , is.numeric)
  N = nrow(num.data);  K = ncol(num.data)
  
  adjust_miss_allRows<-function(miss_data)
  {
    rowise.nMiss = rowSums( is.na(miss_data) )
    NA.allrows = which(rowise.nMiss > subtract( K, 1)  )
    adjust.NAcols = sample(1:K, length(NA.allrows), replace = TRUE)
    for( i in 0:length(adjust.NAcols) ){
      miss_data[NA.allrows[i], adjust.NAcols[i] ] = num.data[NA.allrows[i], adjust.NAcols[i] ]
    }
    bind_cols( dplyr::select_if( complete.data,Negate(is.numeric) ), miss_data )[names(complete.data)]
  }
  if(missing(seed))(seed = 5)
  generate_miss<- function(num.data){
    temp_num_1 = num.data
    for( k in 1:K ){
      set.seed( add(seed,k) )
      samp_indx <- which(temp_num_1[[k]] < quantile(temp_num_1[[k]], probs = 0.5) )
      temp_num_1[sample(samp_indx, ceiling( frac.miss*N ) ), k] <- NA
    }
    temp_num_1 = adjust_miss_allRows(temp_num_1) 
    return(list(miss.data = temp_num_1,
                p_miss = naniar::prop_miss(dplyr::select_if(temp_num_1, is.numeric) ) ) )
  }
  return(generate_miss(num.data) )
}

# simulates missing by replicates
sim_miss_by_reps <- function(df, prop.miss, seed) {
  # df : dataframe without missing values
  # prop.miss: Proportions of missing as in decimals (0.2 for 20%)
  # seed : for randomness
  PivotLong <- df %>% pivot_long_reps()
  unique_nms <- unique(PivotLong$Condition)
  SplitPvtLong <- PivotLong %>% base::split(.$Condition)
  SplitPvtLong <- SplitPvtLong[unique_nms]
  nums <- SplitPvtLong %>% map( ~ select(.x, where(is.numeric)))
  non_nums <- SplitPvtLong %>% map( ~ select(.x, where(Negate(is.numeric))))
  nums_missing <- bind_cols(nums) %>% suppressMessages() %>%
    insert_missingValues(frac.miss = prop.miss, seed = seed) %>%
    magrittr::extract2(1)
  Value <- data.frame(Value = matrix(data.matrix(nums_missing), ncol = 1L))
  sim_miss <- bind_cols(non_nums %>% bind_rows(), Value)
  sim_miss <- select(sim_miss, identifier, Reps, Value) %>%
    pivot_wider(id_cols = identifier,
                names_from = Reps,
                values_from = Value) %>%
    as.data.frame()
  
  return(sim_miss)
  
}


# updates the missing row means by simulating row means from the observed features 
# the means are computed sample by sample

# meansdata : data to be updated
# missdata   : data from which the meansdata has been generated
# mv.rnorm  : TRUE/ FALSE (implies multivariate/ univariate normal dependency between the samples)
adjust_meansMiss<-function(meansdata, missdata, seed = NULL, mv.rnorm = FALSE)
{ 
  if(missing(seed))(seed = 2020) # try 5 small seed 
  means_miss_Idx = which(is.na(meansdata), arr.ind = TRUE )
  # simulating univariate / multivariate normal data if mv.rnorm = TRUE/ FALSE respectively
  sim.means = sim_rnorm(missdata, seed = seed, mv.rnorm ) %>% feature_averages() 
  for ( i in 1:nrow(means_miss_Idx) ){ # updating NaN with simulated values
    meansdata[means_miss_Idx[i,1], means_miss_Idx[i,2]] = sim.means[means_miss_Idx[i,1],means_miss_Idx[i,2]]
  }
  meansdata
}

# MCMC model
impute_mcmc <- function(.data, chains = 4L, iterate = 2000L,
                        warmups = base::floor(iterate/2), factor.vars = NULL, seed = 290, model_name = NULL, ...) {
  if (!missing(factor.vars)) {
    modelname <- paste0(unname(unique(.data[factor.vars])), collapse = "_")
  }
  if (missing(factor.vars)) {
    modelname <- model_name
  }
  stan_data <- get_stanmodel_dt(.data, regroupReps = TRUE, mv.rnorm = FALSE)
  rstan::stan(
    data = stan_data,
    seed = seed,
    cores = parallel::detectCores(),
    control = list(adapt_delta = 0.90),
    verbose = FALSE,
    chain = chains,
    iter = iterate,
    warmup = warmups,
    save_warmup = TRUE,
    model_name = modelname,
    ...
  )
}

# Adjust means fn                        ... mv.rnorm TRUE/ FALSE
get_stanmodel_dt<-function(missdata, regroupReps = FALSE, init.miss = NULL,...)
{ 
  num.data = data.matrix(dplyr::select_if(missdata, is.numeric) )
  if(regroupReps)(num.data = data.matrix(as_experiment_vars(missdata)) )
  if(!regroupReps)(num.data = num.data)
  if(missing(init.miss))(init.miss = 40)
  
  # initializing missing data
  data_matrix = num.data
  data_matrix[is.na(data_matrix)] = init.miss
  # indicator for indexes with missing values
  miss_indicator = num.data[]
  miss_indicator[ ] = foreach( i = 1:dim(miss_indicator)[1], .combine = rbind ) %do%
    if_else(is.na(miss_indicator[i, ]), 0, 1 )
  
  # reshaping data to compute row means
  means.data = feature_averages(missdata)
  if(anyNA(means.data) ){ # Updating NaN in row means data ( priors)
    means.data = adjust_meansMiss(means.data, missdata,...) 
  }
  if(!regroupReps)(means.data = means.data)
  if(regroupReps)( means.data = data.matrix(as_experiment_vars(means.data) ) )
  return( within( list(),{
    K = ncol(num.data )
    N = nrow(num.data )
    poss_obs = which(!is.na(num.data), arr.ind = TRUE)
    poss_miss = which(is.na(num.data), arr.ind = TRUE)
    Nmiss = nrow(poss_miss)
    Nobs = nrow(poss_obs)
    miss_priors = data.matrix(means.data)[which(is.na(num.data) )]
    mu_priors = as.numeric( round( colMeans(num.data, na.rm =  TRUE ), 4 )  )
    data_matrix = data_matrix
    miss_indicator = miss_indicator
  }
  ))
}

# A function that takes in x, a vector of values 
# and returns a value member of x with minimum deviation from the overall mean
min.dev <-  function(x) ( x[ which.min( abs( x -mean(x, na.rm =TRUE ) ) ) ] )

# #
# # The function is clumsy, and can be replaced by simple rstan functions like summary$summary['mean'] if all we want is mean, the theoretical cornerstone of MCMC
# 
# summarise_mcmc1 <- function( rstan.obj, par.obj, fun = 'mean'){
#   purrr::map(rstan.obj, ~ as.data.frame(.x, pars = par.obj ) %>%
#                dplyr::summarise_all(., .funs = get(fun) ) %>% 
#                tidyr::gather(., index, value) ) 
# }
# 

summarise_mcmc <- function(obj,...){
  # summarise Rstan Samples
  # obj, rstan object
  # ..., additional parameters passed on to rstan summary function
  # purrr::map(obj, function(obj){
    rstan::summary(obj, ...)$summary %>%
      as.data.frame() %>%
      dplyr::mutate(param = rownames(.) ) %>%
      dplyr::select( param, mean) %>% tibble::as_tibble()
  # } )
}

list_summarise_mcmc <- function(obj,...){
  # summarise Rstan Samples
  # parallel summary in the sense of list operations
  # obj, rstan object
  # ..., additional parameters passed on to rstan summary function
  purrr::map(obj, \(x){
  rstan::summary(x, ...)$summary %>%
    as.data.frame() %>%
    dplyr::mutate(param = rownames(.) ) %>%
    dplyr::select( param, mean) %>% tibble::as_tibble()
  } )
}



list_update_rstansamples <- function(completed, missing) {
  # A function that take a list of imputed (completed) data set  
  # and the corresponding data set with missing cases
  # updates the status column of the imputed data as imputed or observed.
  purrr::map2(completed, missing, function(.x, .y) {
    .y[is.na(.y[["imputed"]]), "imputed"] <- .x[["mean"]]
    return(.y)
  })
}

# A wrapper of mice model
# produce m copies of imputed data and 
#   returns the average of the m completed data sets, and a list of m imputed sets.

impute_MICE <-function(data,m,...) {
  init.data = data
  num.data = tibble::as_tibble(data) %>% dplyr::select_if(., is.numeric)
  mice_model.data = mice::mice(num.data,m,...)
  list.complete= list()
  for(i in 1:m){
    complete.mice_impute = mice::complete(mice_model.data, i) 
    list.complete[[length(list.complete) + 1 ]] = complete.mice_impute
  }
  mice.imputed = bind_cols( Reduce("+", list.complete )/ length(list.complete), dplyr::select_if(init.data, Negate(is.numeric)) ) %>%
    dplyr::select(names(init.data)) %>% tibble::as_tibble()
  return(list(mice.imputed = mice.imputed, mice_model.data = mice_model.data ) )
}

# Function to extend QRILC imputation model
impute_qrilcLCMD <- function(data, seed = TRUE) {
  temp <- data
  # Qrilc model requires dataframe structures
  num <- dplyr::select(data, where(is.numeric)) |> as.data.frame()
  # initializing random number generator
  set.seed(seed)
  mdl <- imputeLCMD::impute.QRILC(num) %>% extract2(1)
  Out <- bind_cols(mdl, dplyr::select(temp, where(Negate(is.numeric))))[names(temp)]
  Out
}

# use number of numeric and factor variables in a list 
# to determine number of clusters to use in the imputation

clusters_missForest <- function(df){
  p <- select(df, where(is.numeric) | where(is.factor) ) |> names() |> length()
  min ( p, parallel::detectCores() )
}

# Fn to extend missForest imputation model
impute_missForest <- function(data, ...) {
  temp <- as.data.frame(data)
  nms <- names(temp)
  nms.imp <- names(temp[, sapply(temp, class) %in% c("numeric", "factor")])
  nms.char <- base::setdiff(nms, nms.imp)
  impDt <- dplyr::select(temp, tidyr::all_of(nms.imp))
  non.impDt <- dplyr::select(temp, tidyr::all_of(nms.char))
  Out <- missForest::missForest(impDt, ...)
  dplyr::bind_cols(Out$ximp, non.impDt)[nms]
}


# use missing cases to 
# declares the status of rows of the imputed data set 
# as either observed/ or imputed

# assumes long format with single numeric column

# new proposed name :  update_status_long ()
update_status <- function(completed, missing){
  # Uses indexes of the missing cases in missing data to designated values of a newly created column
  # status as observed or missing
  completed = dplyr::mutate(completed, status = 'observed' )
  ind.miss = which( is.na(dplyr::select_if(missing, is.numeric) ))
  completed[ ind.miss, 'status'] = 'imputed'
  return( completed)
}

# derives numerical constant to facilitate log transformation
# x :        the input data
norm_const_log<-function(x,...){
  min.value = min(x,...)
  if(min.value > 1 )(num.const_x = 0 )
  if(min.value <= 1 )(num.const_x = add( abs(min.value), 1 ) %>% round(., 1) )
  return( num.const_x )
}

#
log_transform<-function(x, norm.const = NULL){
  if(missing(norm.const))( norm.const = 0 )
  log( magrittr::add( x, norm.const ) )
}

#
negate_log.transform<-function(x, norm.const = NULL){
  if(missing(norm.const))( norm.const = 0 )
  magrittr::raise_to_power( exp(1), x) %>% 
    magrittr::subtract(., norm.const)
}

#
# normalised_RMSE<-function(obs.x, est.x){
#   if(require(magrittr)){
#     sq.Err = raise_to_power(subtract(obs.x, est.x), 2)
#   } 
#   sqrt( mean(sq.Err) ) / sd(obs.x) 
# }

normalised_RMSE <- function(obs.x, est.x) {
  sq.Err = (obs.x - est.x) ^ 2
  sqrt(mean(sq.Err)) / sd(obs.x)
}


# # data specific wrapper function for correlation statistic
# 
# correlation <- function(.data) {
#     tempData <- .data
#     # split data by samples/ conditions & re-order the list by sample names
#     splitData <- tempData %>% split(.$samples)
#     splitData <- splitData[unique(tempData$samples)] 
#     # create rectangular with samples/ conditions as column major variables
#     sampleValues <- splitData %>% 
#                     purrr::map(~dplyr::select(.x, .data$log_values)) %>% purrr::reduce(cbind) 
#     names(sampleValues) <- unique(tempData$samples)
#     # calculates pearson correlation between samples 
#     corrData <- corrr::correlate(sampleValues, quiet = TRUE) %>% 
#     
#     corrr::stretch() %>%  dplyr::filter( !is.na(r))
#     # re-group the statistic summary with other variables of interest
#     dplyr::bind_cols(dplyr::select_if(splitData[[1]], 
#                                       Negate(is.numeric) )[1:dim(corrData)[1], ], corrData)
#   
# }


# plotting functions

# Customized theme for plots
customize_themePlot <- function() {
  theme(plot.title=element_text(hjust = 0.5, size = 10, face="plain"),
        legend.text = element_text(size = 7, hjust = 0.7),
        legend.title = element_blank(),
        legend.position = "bottom",
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 7)
  )
}



# plots distribution of data into two facet: 
#     used for plotting density of original into data name against conditions (samples)

ggplot_DensityCurve <-function(data,x,facet1,facet2,title){
  require(rlang)
  args = match.call()
  ggplot(data, aes({{x}}) ) + 
    geom_density(linewidth = 0.8 ) + 
    facet_grid( eval(args$facet1) ~ eval(args$facet2) , scales = "free") + 
    theme_classic() +
    ggtitle({{title}} ) +  
    customize_themePlot()
}

# plots the distribution of data into three facets and one grouping var ( color option ) ; 
#       used for plotting data density into pct missing against data name and conditions (samples)

ggplot_CompareDensityCurves<-function(data, xVar, groupVar, facetVar1, facetVar2, titleVar, xlab=NULL, ylab=NULL )
{
  require(rlang)
  args = match.call()
  .title = unique(dplyr::select(data, {{titleVar}} )  )
  if(missing(xlab)) xlab = "Log values"
  if(missing(ylab)) ylab = "Density"
  
  data %>% 
    ggplot(data = . ) +  aes({{xVar}} , colour = {{groupVar}})  + 
    geom_density(linewidth = 0.5) + 
    facet_grid( eval(args$facetVar1) ~ eval(args$facetVar2) , scales = "free") + 
    theme_classic() +
    labs(title = .title, y = ylab, x = xlab ) + 
    customize_themePlot()
}

#
ggplot_CompareBoxPlots<-function(data, xVar, yVar, facetVar1, titleVar)
{
  require(rlang)
  .title = unique(dplyr::select(data, {{titleVar}} )  )
  args = match.call()
  data %>% 
    ggplot(data = . ) +  aes({{xVar}}, {{yVar}} , fill = {{xVar}} )  + 
    geom_boxplot( size = 0.8 , notch = FALSE ) + 
    scale_x_discrete(guide = guide_axis(n.dodge = 2) )+
    scale_y_continuous( labels = scales::number_format(accuracy = 0.01) ) +
    facet_grid( . ~ eval(args$facetVar1) , scales = "free") + 
    theme_classic() +
    # theme_cleveland() + 
    # theme(plot.title=element_text(hjust = 0.5, size = 13, face="bold") ) +
    # #axis.text.x = element_text(hjust = 1 ) ) + 
    labs(title = titleVar) + #, y = ylab, x = xlab ) + 
    customize_themePlot() +
    theme(axis.text.x = element_text(hjust = 1 ) )
  # + 
  #   ggeasy::easy_remove_legend_title() + ggeasy::easy_legend_at(to = "bottom")
}

# takes two coordinates x and y to plot lines into two facets and one grouping variables (color option)
#         used for comparing MAE and NRMSE for models (grouping var) for each data vs sample conditions
ggplot_ErrorLines<-function(data, xVar, yVar, groupVar,facetVar, xlab=NULL, ylab=NULL, titleVar = NULL)
{ 
  require(rlang)
  .title = names(dplyr::select(data, {{yVar}}) )
  if(!missing(titleVar))(.title = paste0(.title," ",titleVar) )
  if(missing(xlab)) xlab = "Proportion  missing"
  if(missing(ylab)) ylab = "Imputation error"
  args = match.call()
  ggplot(data, aes(x = {{xVar}}, y = {{yVar}}, colour = {{groupVar}} ) )  + 
    geom_point( aes(shape = {{groupVar}} ) , size = 3 )  + 
    geom_line(linewidth = 0.5 ) + 
    scale_x_continuous( labels = scales::percent_format() ) +
    scale_y_continuous( trans = "log", labels = scales::number_format(accuracy = 0.01) ) + 
    facet_grid(eval(args$facetVar) ~ . , scales = "fixed" ) + 
    theme_classic() +
    customize_themePlot()+
    labs(title = .title, y = ylab, x = xlab ) 
}

# extract jpeg graphs from plotted data
# plot.object              : R object created as a result of plotting data
# names.graph             : name(s) given to the plot to be loaded
# dir                     : the folder/ directory to which a plot is saved. 

#                         ".../paper_writeup_andFiles/plots/" is the the default directory 

extract_plots<-function( plot.object, names.graph, 
                         nrow = NULL, height = NULL, width = NULL, 
                         dir = fs::path('paper_writeup_andFiles', 'plots') )
{
  if(missing(nrow))( nrow = 1)
  if(missing(height)) (height = 1050 )
  if(missing(width)) (width = 1850 )
  grid.arrange(plot.object, nrow = nrow )
  # res (resolution) corresponds to dpi. Its for JPEG and PNG devices
  # quality -compression for JPEG devices
  dev.copy(jpeg, here::here(dir, paste0(names.graph,'.jpeg') ), 
           res = 300, quality = 100, height = height, width = width ) 
  dev.off()
}

#
ggplot_VarLines<-function(data, xVar, yVar, groupVar, titleVar, xlab=NULL, ylab=NULL )
{ 
  require(rlang)
  .title = unique(dplyr::select(data, {{titleVar}} ) )
  if(missing(xlab)) xlab = "Proportion  missing"
  if(missing(ylab)) ylab = "Variance component "
  args = match.call()
  ggplot(data, aes(x = {{xVar}}, y = {{yVar}} , group = {{groupVar}}, colour = {{groupVar}}  ) )  + 
    geom_line( linewidth = 0.8 )  + 
    geom_point( aes(shape = {{groupVar}}) , size = 3 )  + 
    scale_x_continuous( labels = scales::percent_format() ) +
    scale_y_continuous(trans = "log", labels = scales::number_format(accuracy = 0.01 ) ) +
    theme_classic() +
    customize_themePlot() +
    labs( title = .title, y = ylab, x = xlab ) 
}



# diagnostic and inference plots2
# trace plots for MCMC samples

# stan_obj, stan_fit (model object)
# n_chains, number of chains used to generate the stan_obj
# title, desired plot title
# ..., other arguments (pars) passed to stan_trace, and 
#       optional plot arguments (linetype, size, alpha) passed to geom_path

TracePlot <- function(stan_obj,  n_chains, title, add_burn_in = TRUE, ...) 
{
  suppressMessages(
    rstan::stan_trace(stan_obj, inc_warmup = add_burn_in, ...) +
      ggtitle({{ title }}) +
      scale_color_manual(values = rep("black", n_chains)) +
      theme_classic() +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 13, face = "bold")
      )
  )
}


# sampler parameter and divergence

get_sampler_params <- function(fitObj) {
  # lists by number chains used
  # combine chains together into a tbl
  sampler_param <-
    rstan::get_sampler_params(fitObj, inc_warmup = FALSE) %>%
    do.call("rbind", .)
  
  logPosterior <-
    rstan::get_logposterior(fitObj, inc_warmup = FALSE) %>%
    unlist(recursive = TRUE)
  
  sampler_param <-
    dplyr::bind_cols(sampler_param,
                     logPosterior = logPosterior,
                     Divergent = dplyr::if_else(get_divergent_iterations(fitObj), "Divergent", "Not Divergent")
    )
  
  sampler_param
}

# plots Mean Metropolis acceptance rate Vs divergence (Divergent / Not Divergent)
# samplerObj, generated by get_sampler_params
# 
sampler_divergence <- function(samplerObj) {
  ggplot2::ggplot(samplerObj, aes(x = Divergent, y = accept_stat__)) +
    geom_violin() +
    theme_classic() +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 13,
        face = "bold"
      ),
      axis.title = element_text(face = "bold")
    ) +
    labs( title = "Divergence",
          y = ("Mean Metropolis Accepatance Rate"),
          x = NULL
    )
}

# summarizes a list of data 
# Provides a basic counts of observed/missing data grouped by file names and sample/conditions
list_gather_summarise <- function(.data, what = c("missing", "complete")) {
  Dt <- .data %>%
    map2(., names(.), ~ update_list(.x, "dataName" = .y)) %>%
    map_df(., \(x) {
      x %>%
        dplyr::select(names(dplyr::select_if(., is.numeric)), dataName) %>%
        tidyr::gather(key = condition, value = N, -"dataName") %>%
        mutate(condition = gsub("[A-z]_||[A-z]||_[1-9]", "", .data$condition)) %>%
        dplyr::group_by(dataName, condition) 
    }) 
  if (what == "missing") {
    return(
      Dt %>% 
      dplyr::summarise(nFeatures = n(), nMiss = naniar::n_miss(.data$N),
                pct_miss = naniar::pct_miss(.data$N) )
    )
  }
  if (what == "complete") {
    
    return(Dt %>% 
              dplyr::summarise(nFeatures = n())  )
  }
}


# Create Cluster For Parallel MissForest

set_parallelMissForest <- function(df, free.cores = 0L) {
  ncl = magrittr::subtract(clusters_missForest(df), free.cores)
  if (.Platform$OS.type == "windows") {
    return(parallel::makePSOCKcluster(ncl) )
  }
  parallel::makeForkCluster(ncl) #for running in linux based OS
}

# Imputes Parallel MissForest

par_imputeMissForest <- function(data, free.cores = 0L, ...) {
  if (is.data.frame(data)) {
    dArgs <- list(...)
    cl <- set_parallelMissForest(df, free.cores = free.cores)
    doParallel::registerDoParallel(cl)
    Out <- exec("impute_missForest", data, !!!dArgs)
    parallel::stopCluster(cl) |> suppressMessages()
    cl <- c()
    return(Out)
  }
  if (!is.list(data)) stop("use a list or dataframe")
  indx <- map(data, \(x)NCOL(x)) |> list_c() |> which.min()
  cl <- set_parallelMissForest(data[[indx]], free.cores = free.cores)
  doParallel::registerDoParallel(cl)
  Out <- furrr::future_map(data, \(x)impute_missForest(x, ...),
                           .options = furrr_options(seed = TRUE), .progress = TRUE
  )
  parallel::stopCluster(cl) |> suppressMessages()
  cl <- c()
  return(Out)
}



# Impute Parallel MICE

par_imputeMice <- function(ldf, diagnostic = 2L, ...) {
  future::plan(multisession, workers = parallel::detectCores())
  Out <- furrr::future_map(ldf, \(x)impute_MICE(x, ...) %>% magrittr::extract2(diagnostic),
                           .options = furrr_options(seed = TRUE), .progress = TRUE
  )
  on.exit(future::plan(sequential), add = TRUE) |> suppressMessages()
  Out
}

# Impute Parallel Qrilc
par_imputeQrilc <- function(ldf, seed = TRUE,...) {
  future::plan(multisession, workers = parallel::detectCores())
  Out <- furrr::future_map(ldf, \(x) impute_qrilcLCMD(x, seed = seed, ...),
                           .options = furrr_options(seed = TRUE), .progress = TRUE
  )
  on.exit(future::plan(sequential), add = TRUE) |> suppressMessages()
  Out
}

# parallel imputation for mice, Qrilc, and missForest models
# lst- a list of missing data
# .paramList - a list of parameters for each of the models; name of each list of parameter must correspond to the model for which it is specified

# seed - passed to the corresponding model function, especially an an option in parallel execution

impute_non_mcmc <- function(ldf, parList) {
  map2(parList, names(parList), \(p, m){
    exec(paste0("par_impute", m), ldf, !!!p)
  }) |> list_flatten()
}

#list_ByNamesBindRows

# Bind list by rows, adding a column of name of each list.
# names_to = "models", the default column name

list_bind_rows <- function(.lst, names_to = "models") {
  df <- purrr::map(.lst, as.data.frame) %>% do.call("rbind", .)
  df[names_to] <- gsub("[0-9.]", "", rownames(df))
  tibble::as_tibble(df)
}


is_numeric_pivot_long <- function(dt, ...) {
  is.num.cols <- sapply(dt, class) == "numeric"
  nms.nums <- names(dt[, is.num.cols])
  Out <- tidyr::pivot_longer(
    data = dt,
    cols = tidyselect::all_of(nms.nums), ...)
  Out
}

# takes a wide data and gathers numeric columns into long format 
# ... , arguments passed to dplyr::gather

is_numeric_gather <- function(dt, ...) {
  nms <- names(dt)
  nums.cols <- sapply(dt, class) == "numeric"
  nms.nums <- names(dt[, nums.cols])
  non.nums.nms <- dplyr::setdiff(nms, nms.nums)
  tidyr::gather(dt, -tidyselect::all_of(non.nums.nms), ...)
}


split_by <- function(dt, ...) {
  split(dt, list(...) )
}


is_numeric_join <- function(dt1, dt2, use_var = 2L, join_fun = "left_join", ...) {
  FUN <- get(join_fun)
  nms <- names(dt2)
  non.numcols <- sapply(dt2, class) != "numeric"
  nms.non.num <- function()names(dt2[ , non.numcols])
  # FUN(dt1, dt2, by = tidyselect::all_of(nms.non.num))
  FUN(dt1, dt2, by = nms.non.num())
}



# combine all model results into a list of lists and  aligns mcmc results with test/ simulated data 

# other_models , a list of non_mcmc_imputed (Qrilc, mice, and missForest)
# mcm_sampl, a list of mcmc_samples
# sim_miss, a simulated missing 

# combine_imputed_models1 <- function(other_models, mcmc_sampl, sim_miss) {
#   combined <- c(other_models, MCMC = list(sim_miss)) %>%
#     purrr::map(~ purrr::map(., ~ is_numeric_gather(., key = "samples", value = "imputed")))
#   
#   # align mcmc_sampl with sim_miss
#   names_by_sim <- names(combined$MCMC)
#   mcmc_sampl <- mcmc_sampl[names_by_sim]
#   # add mcmc_sampl
#   combined[["MCMC"]] <-
#     list_update_rstansamples(imp.data = mcmc_sampl, miss.data = combined[["MCMC"]])
#   
#   combined
# }



update_rstansamples <- function(completed, missing) {
  # takes a data frame of imputed (completed) data set  
  # and the corresponding data set with missing cases
  # adds imputed mcmc values to the missing indexes
  miss_idx <- is.na(missing[['imputed']])
  missing[miss_idx, 'imputed'] = completed[['mean']]
  missing
}

merge_mcmc <- function(mcmc, miss, miss.tidy = FALSE) {
  # Adds MCMC imputed samples to missing cases in the original dataset with missing values
  # mcmc, imputed mcmc values
  # miss, original dataset with missing values
  # miss.tidy = FALSE (default) converts miss to a long tidy format
  if (miss.tidy) {
    miss <- miss
  } else {
    miss <- is_numeric_gather(miss, key = "samples", value = "imputed")
  }
  Out <- update_rstansamples(completed = mcmc, missing = miss)
  Out
}

# Use either 
#' merge_mcmc / multi_merge_mcmc to combine mcmc samples
#' Use merge models to combine both mcmc and non mcmc imps
#' No function to merge non mcmc samples only.
#' You may skip merging non mcmc (only) and use either process_models and list_process_models

# adds raw rstan samples (miss_data parameter values ) to missing cases 
multi_merge_mcmc <- function(mcmc, miss, miss.tidy = FALSE, order.names = TRUE) {
  # A convenient function to merge MCMC imputed Missing values to the original list of Data Sets With Missing Values
  # mcmc - a list of imputed mcmc samples
  # miss - list of data sets with missing values
  # miss.tidy = TRUE/FALSE transforms missing data sets to a long tidy format (default)
  # order.names = TRUE/FALSE Orders a list of mcmc.sampl by names of miss
  if (miss.tidy) {
    miss <- miss
  } else {
    miss <- purrr::map(miss, \(x) is_numeric_gather(x, key = "samples", value = "imputed"))
  }
  if (order.names) {
    lstnms <- names(miss) # align mcmc_sampl with list names for miss
    mcmc <- mcmc[lstnms]
    Out <- list_update_rstansamples(completed = mcmc, missing = miss)
  } else {
    Out <- list_update_rstansamples(completed = mcmc, missing = miss)
  }
  Out
}


# other_imp_model=NULL (default), otherwise provide two level list (a list of list of models)
# Each list of model either has a single/ several list(s) of imputed data
# mcm_sampl , imp from mcmc model
# miss,   the missing data set that's to be imputed

multi_tidy_merge <- function(mcmc, miss, nonmcmc = NULL, ...) {
  # Inserts MCMC samples into missing cases of miss dataset
  # then concatenates completed MCMC dataset with completed datasets from non MCMC models if
  # such completed datasets by such models are provided as nonmcmc.sampl
  mcmc <- multi_merge_mcmc(mcmc = mcmc, miss = miss, ...)
  if (missing(nonmcmc)) {
    Out <- mcmc
  } else {
    nonmcmc <- nonmcmc %>%
      purrr::map(\(x) is_numeric_gather(x, key = "samples", value = "imputed"))
    mcmc <- stats::setNames(mcmc, paste("MCMC", base::names(mcmc), sep = "_"))
    Out <- c(nonmcmc, mcmc)
  }
  Out
}

# takes data frames of the respective inputs

summarize_imputed <- function(miss, non_miss, combined_imp, summary_digits = 2L, mcmc = TRUE) {
  # same long formats as combined_imp
  # add/ join original values from non_miss_test dt
  # update status, imputed/ observed

  is_nums_gather__1 <- \(dt) { # value = "original"
    is_numeric_gather(dt, key = "samples", value = "original")
  }
  is_nums_gather__2 <- \(dt) { # value ="imputed"
    is_numeric_gather(dt, key = "samples", value = "imputed")
  }
  long_sim_miss_ <- is_nums_gather__2(miss)
  if (!mcmc) {
    combined_imp <- is_nums_gather__2(combined_imp)
  }
  non_missLong <- is_nums_gather__1(non_miss)
  combine_updated <- update_status(imputed_data = combined_imp, missdata = long_sim_miss_) %>%
    is_numeric_join(dt..1 = ., dt..2 = non_missLong) %>%
    dplyr::mutate(abs.error = abs(subtract(original, imputed)) %>% round(summary_digits))
  
  combine_updated
}

# list_update_status <- function(completed, missing, dtname){
#   # completed: imputed list data
#   # missing: tidy missing list data
#   # dtname: data name from the missing list data
#   
#   nms <- tidyselect::vars_select(names(completed), vroom::contains(dtname))
#   Out <- BiocGenerics::lapply(completed[nms], FUN = update_status, missing=OrigMissTidy[[dtname]])
#   Out
# }

list_update_status <- function(completed, missing, dtname, non.miss = NULL,  miss.tidy = FALSE, nonmiss.tidy = FALSE,sim= FALSE ,digits = 2L){
  # completed: imputed list data
  # missing: tidy missing list data
  # dtname: data name from the missing list data
  
  if (miss.tidy) {
    miss <- missing
  } else {
    miss <- BiocGenerics::lapply(missing, is_numeric_gather, key = "samples", value = "imputed")
    # missing <- purrr::map(miss, \(x) is_numeric_gather(x, key = "samples", value = "imputed"))
  }
  
  nms <- tidyselect::vars_select(names(completed), vroom::contains(dtname))
  DtnmMiss <- vars_select(names(miss), vroom::contains(dtname) )
  DtnmComp <- vars_select(names(completed), vroom::contains(DtnmMiss))
  
  if(sim){
    
    Out= purrr::map2(completed[DtnmComp], miss[DtnmMiss], \(c, m)update_status(completed = c, missing = m))
  } else {
    Out <- BiocGenerics::lapply(completed[nms], FUN = update_status, missing=miss[[dtname]])
  }
  
  if(!missing(non.miss)){
    
    if(nonmiss.tidy){
      nonmiss <- non.miss[dtname]
    } else {
      nonmiss <- BiocGenerics::lapply(non.miss[dtname], FUN = is_numeric_gather, key = "samples", value = "original")
      # non.miss <- purrr::map(non.miss, \(x)is_numeric_gather(x, key = "samples", value = "original") )
    }
    
    join_nonmiss_data <- function(xcompleted, xnonmiss){
      
      res <- is_numeric_join(dt1 = xcompleted, dt2 = xnonmiss) |>
        dplyr::mutate(abs.error = abs(subtract(original, imputed)) |>  round(digits)) |>
        tibble::as_tibble()
      res
      
    }
    
    Out <- BiocGenerics::lapply(Out[DtnmComp], FUN = join_nonmiss_data, xnonmiss=nonmiss[[dtname]])
    # Out <- purrr::map2(Out[DtnmComp], miss[DtnmMiss], \(c, m)join_nonmiss_data(xcompleted = c, xnonmiss = m)
  }
  
  Out
  
}






# adds original values to the imputed data
# adds status column /imputed/observed

# sim_miss, wide format simulated/ original missing data
# non_miss_test, wide format non missing test data 
# combined_imp, long format imputed data

multi_process_mcmc <- function(miss, non_miss, combined_imp, add_summary = TRUE, print_digits = 2L) {
  long_miss <- purrr::map(miss, \(x) {
    is_numeric_gather(x, key = "samples", value = "imputed")
  })
  long_non_miss <- purrr::map_dfr(non_miss, \(dt) {
    is_numeric_gather(dt, key = "samples", value = "original")
  })
  comb_updated <- purrr::map2(combined_imp, long_miss, \(x, y) {
    update_status(imputed_data = x, missdata = y)
  }) %>%
    purrr::map_dfr(\(dt) {
      is_numeric_join(dt..1 = dt, dt..2 = long_non_miss)
    })
  if (add_summary) {
    comb_updated <- dplyr::mutate(comb_updated, abs.error = abs(subtract(original, imputed)) %>% round(print_digits))
  }
  comb_updated
}

any_missing <- function(x) {
  # x,  a vector of values
  # returns TRUE if a a vector has a missing value
  apply(x, 1 , anyNA)
}

remove_completely_miss <- function(dt, missing.rows = FALSE) {
  # uses row means data to remove rows of replicates that  are completely missing
  # dt -untidy data frame with missing values
  # missing.rows = TRUE / FALSE, returns rows of features, in each experimental condition, that has at least a data point
  row_mns <- feature_averages(dt) # feature averages by experimental conditions
  rmiss <- any_missing(row_mns)
  if (missing.rows)
    return(dt[rmiss,])
  dt[!rmiss,]
}


# sums up replicate values (by features) into a single (total) quantity
#  rep_groups, conditions
#  value, the variable holding numeric quantity

get_replicate_sums <- function(.data, rep_groups, value) {
  temp_df <- .data %>%
    dplyr::group_by(.data[[rep_groups]]) %>%
    dplyr::mutate(Totals = sum(.data[[value]], na.rm = TRUE)) %>%
    dplyr::ungroup()
  indx_repl_sum <- !base::duplicated(temp_df$Totals)
  temp_df <- temp_df[indx_repl_sum,]
  temp_df
}

# Can be used to automate finding rep_groups in get_replicate_sums function

get_numeric_names <- function(dt) {
  num_nms <- names(dt[, sapply(dt, class) == "numeric"])
  unique(readr::parse_number(num_nms))
}

# Bootstrap test for constant variance :fligner.test for constant variance
# returns average p value

boot_var_test <- function(dt, values, category, iter = 1000) {
  tot_pvalue = 0
  B = 1
  repeat{
    temp_dt <- dplyr::slice_sample(dt, prop = 1, replace = TRUE)
    p_value <- fligner.test(temp_dt[[values]] ~ temp_dt[[category]], data = temp_dt)$p.value
    tot_pvalue <- p_value + tot_pvalue
    B <- B + 1
    if (B > iter) break
  }
  avg_pval <- tot_pvalue / iter
  avg_pval
}

# Bootstrap test for constant variance
# Returns Bootstrap quantile intervals for the average p values 

var_test_intervals <- function(dt, values, category, iter=50, N =1000) {
  n = 1
  avg_pvalue <- vector(length = N)
  
  while (n <= N) {
    tot_pvalue = 0
    B = 1
    repeat {
      temp_dt <- dplyr::slice_sample(dt, prop = 1, replace = TRUE)
      p_value <- fligner.test(temp_dt[[values]] ~ temp_dt[[category]], data = temp_dt)$p.value
      tot_pvalue <- p_value + tot_pvalue
      B <- B + 1
      if (B > iter) break
    }
    avg_pvalue[n] <- tot_pvalue / iter
    n = n + 1
  }
  quantile(avg_pvalue, probs = c( 0.025, 0.975) )
}

# returns Bootstrap pooled mean and variances
# group_var, { status (imputed/ observed) and Conditions( 25, 50 , 100 ) } two grouping variables to over which the summaries are pooled

boot_pooled_summary <- function(data, group_vars, values, iter ) {
  stopifnot(length(group_vars) == 2)
  grp_1_nm <- group_vars[1]
  grp_1_val <- unique(data[[grp_1_nm]])
  grp_1_len <- length(grp_1_val)
  mat_res <- array(rep(0, grp_1_len * 3), dim = c(grp_1_len, 3), dimnames = list(grp_1_val, c("Means", "Var_p", "Sd_p")))
  B <- 0
  while (B < iter) {
    temp_res <- dplyr::slice_sample(data, prop = 1, replace = TRUE) %>%
      dplyr::group_by(across(tidyselect::all_of(group_vars))) %>%
      dplyr::summarise(Means = mean(.data[[values]]), Var = var(.data[[values]]), Sd = sqrt(Var), n_adj = n() - 1, SS = n_adj * Var, .groups = "drop_last") %>%
      dplyr::summarise(Means = mean(Means), Var_p = sum(SS) / sum(n_adj), Sd_p = sqrt(Var_p), .groups = "drop")
    
    temp_res_num  <- dplyr::select(temp_res, where(is.numeric))
    mat_res <- mat_res + as.matrix(temp_res_num )
    B <- B + 1
  }
  mat_res <- as.data.frame(mat_res / iter)
  rownames(mat_res) <- temp_res[[grp_1_nm]]
  mat_res <- tibble::rownames_to_column(mat_res, var = grp_1_nm)
  mat_res 
}

# returns bootstrap mean of x
boot_mean <- function(x, iter = 10000, seed = NULL, ...){
  B= 0
  temp_res = 0
  while(B < iter){
    if(!missing(seed))set.seed(seed + B)
    temp_x <- sample(x, size = length(x), replace = TRUE)
    temp_res = temp_res + mean(temp_x, ...)
    B = B + 1
  }
  temp_res/ iter
}

# returns bootstrap Normalized RMSE
boot_normalised_RMSE <- function(obs.x, est.x, iter = 1, seed = NULL) {
  B <- 0
  res <- 0
  temp_df <- cbind.data.frame(obs.x, est.x)
  require(dplyr)
  repeat{
    B <- B + 1
    if (B > iter) break
    if (!missing(seed)) set.seed(seed + B)
    df <- dplyr::slice_sample(temp_df, prop = 1, replace = TRUE)
    sq.Err <- (df[['obs.x']] - df[['est.x']])^2
    temp_res <- sqrt(mean(sq.Err)) / sd(df[['obs.x']])
    res <- res + temp_res
  }
  res / iter
}


boot_pooled_summary_2 <- function(data, group_var, values, iter, seed = NULL) {
  mat_res <- array(rep(0, 3), dim = c(1, 3), dimnames = list(NULL, c("Means", "Var_p", "Sd_p") ) )
  B <- count <- 0
  while (B < iter) {
    if(!missing(seed))set.seed(seed + B)
    temp_df <- dplyr::group_by(data, across(tidyselect::all_of(group_var)))  %>%
      dplyr::slice_sample(prop = 1, replace = TRUE)
    temp_res <- temp_df %>%
      dplyr::summarise(Means = mean(.data[[values]]), Var = var(.data[[values]]), Sd = sqrt(Var), n_adj = n() - 1, SS = n_adj * Var, .groups = "drop") %>%
      dplyr::summarise(Means = mean(Means), Var_p = sum(SS) / sum(n_adj), Sd_p = sqrt(Var_p))
    
    mat_res <- mat_res + as.matrix(temp_res)
    B <- B + 1
  }
  as.data.frame(mat_res / iter)
}






